"""Context computed from the sequence alone, rendered for Jev and as a matrix for the controls.

Level 1: physico-chemical properties and amino-acid composition.
Level 2: level 1 plus matches to a fixed list of textbook sequence motifs.

The motif list was fixed from general biochemistry before any context run and is not tuned on the
benchmarks. Descriptions give the structural or chemical role a motif is known for; none names an
EC class. Short patterns also match by chance, and the prompt says so.
"""
import re

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

AA = "ACDEFGHIKLMNPQRSTVWY"

# (key, regex, description shown to the model)
MOTIFS = [
    ("rossmann", r"G.G..G", "GxGxxG glycine-rich loop, typical of NAD(P)/FAD-binding Rossmann folds"),
    ("p_loop", r"[AG].{4}GK[ST]", "Walker A / P-loop [AG]x4GK[ST], binds the phosphates of ATP or GTP"),
    ("gxsxg", r"G.S.G", "GxSxG 'nucleophile elbow' carrying a catalytic serine (alpha/beta-hydrolase fold)"),
    ("hexxh", r"HE..H", "HExxH zinc-binding motif, zinc-dependent catalytic site"),
    ("kinase_loop", r"H[RL]D[LIVM][AK].{2}N", "HRDLKxxN catalytic loop of protein-kinase-like folds"),
    ("aars_high", r"H[IMLV]G[HN]", "HIGH motif of class I aminoacyl-tRNA synthetase-like nucleotidylyl folds"),
    ("aars_kmsks", r"K[MIL][SA][KR]S", "KMSKS motif of class I aminoacyl-tRNA synthetase-like folds"),
    ("p450_heme", r"[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]", "cysteine heme-iron ligand loop of cytochrome P450s"),
    ("sdr_yxxxk", r"Y[PSTAGNCV][STAGNQCIVM][STAGC]K", "YxxxK catalytic tyrosine/lysine pair of short-chain dehydrogenases/reductases"),
    ("thioredoxin", r"WC[GAP][PH]C", "WCGPC redox-active disulfide of thioredoxin-like folds"),
    ("dead_box", r"DE[AH][DH]", "DEAD/DEAH box of RNA/DNA helicases (ATP hydrolysis)"),
    ("amp_binding", r"[STG][STAG]G[ST][STEI][SG].[PASLIVM][KR]", "AMP-binding motif of adenylate-forming enzymes (acyl-AMP intermediates)"),
    ("radical_sam", r"C...C..C", "CxxxCxxC cysteine triad that binds a [4Fe-4S] cluster in radical SAM enzymes"),
    ("plp_lysine", r"[ST].{2}K[ST][LIVMFYW]", "PLP-binding lysine context of pyridoxal-phosphate enzymes"),
]


def physchem(seq: str) -> dict:
    pa = ProteinAnalysis(seq)
    return {
        "length": len(seq),
        "mw_kda": pa.molecular_weight() / 1000,
        "pi": pa.isoelectric_point(),
        "charge_ph7": pa.charge_at_pH(7.0),
        "gravy": pa.gravy(),
        "aromaticity": pa.aromaticity(),
        **{f"aa_{a}": seq.count(a) / len(seq) for a in AA},
    }


def motif_hits(seq: str) -> dict:
    """key -> list of 1-based start positions (non-overlapping regex matches)."""
    return {k: [m.start() + 1 for m in re.finditer(rx, seq)] for k, rx, _ in MOTIFS}


def matrix(seqs, level: int) -> tuple[np.ndarray, list]:
    """Feature matrix for the logistic-regression control at a context level."""
    rows, names = [], None
    for s in seqs:
        f = physchem(s)
        f["log_length"] = float(np.log(f.pop("length")))
        if level >= 2:
            f.update({f"motif_{k}": float(len(v) > 0) for k, v in motif_hits(s).items()})
        names = names or list(f)
        rows.append([f[n] for n in names])
    return np.array(rows), names


def render(seq: str, level: int) -> str:
    """Text context appended after the sequence in Jev's state."""
    if level == 0:
        return ""
    f = physchem(seq)
    comp = ", ".join(f"{a} {100 * f['aa_' + a]:.1f}" for a in AA)
    out = [
        "",
        "Properties computed from the sequence:",
        f"- length {f['length']} residues, molecular weight {f['mw_kda']:.1f} kDa",
        f"- isoelectric point {f['pi']:.2f}, net charge at pH 7 {f['charge_ph7']:+.1f}",
        f"- mean hydropathy (GRAVY) {f['gravy']:+.3f}, aromatic residues {100 * f['aromaticity']:.1f}%",
        f"- amino-acid composition (%): {comp}",
    ]
    if level >= 2:
        hits = motif_hits(seq)
        found = [(k, desc, hits[k]) for k, _, desc in MOTIFS if hits[k]]
        out.append("")
        out.append(f"Sequence motif scan ({len(MOTIFS)} textbook patterns; short patterns also occur by "
                   f"chance, so a single match is weak evidence):")
        if found:
            for _, desc, pos in found:
                where = ", ".join(map(str, pos[:5])) + (" ..." if len(pos) > 5 else "")
                out.append(f"- {desc}: {len(pos)} match{'es' if len(pos) > 1 else ''} (starting at residue {where})")
        else:
            out.append("- none of the screened motifs matched")
        absent = [desc.split(",")[0].split(" (")[0] for k, _, desc in MOTIFS if not hits[k]]
        if found and absent:
            out.append(f"- not found: {'; '.join(absent)}")
    return "\n".join(out)
