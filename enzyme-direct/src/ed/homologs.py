"""Homolog evidence for Jev, following ../enzyme-evidence (re-implemented here, not imported).

MMseqs2 search of each benchmark protein against CARE's training set (Swiss-Prot enzymes) with the
benchmark proteins removed; the top-k hits become one candidate per exact EC with its evidence.
Rendering matches enzyme-evidence's prompt: candidates shuffled per protein with opaque letters,
no accession, entry name or sequence, plus an 'N' (none) option.
"""
import hashlib
import random
import string
import subprocess
import urllib.request

import pandas as pd

from . import config, data

FMT = "query,target,fident,qcov,tcov,qlen,tlen,evalue,bits"
NONE_LABEL = "N"
INSTRUCTIONS = (
    "Which candidate enzyme function (EC number) is best supported by the homologs retrieved "
    "for this protein? Support is stronger with higher sequence identity, higher alignment "
    "coverage, a better (lower) hit rank, and more agreeing homologs. Choose N if no candidate "
    "is convincingly supported."
)


def reference(cfg, exclude: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    """CARE training proteins minus `exclude` (by accession and sequence); Entry -> [ECs]."""
    rel = cfg["controls"]["nn_reference"]
    path = data.fetch_care(cfg) / rel
    if not path.exists():
        urllib.request.urlretrieve(f"{cfg['care']['raw_base']}/{cfg['care']['commit']}/{rel}", path)  # nosec
    ref = pd.read_csv(path, usecols=["Entry", "Sequence", "EC number"])
    ref = ref[~ref.Entry.isin(exclude.Entry) & ~ref.Sequence.isin(exclude.Sequence)]
    ecs = ref.groupby("Entry")["EC number"].agg(lambda s: sorted(set(s)))
    return ref.drop_duplicates("Entry")[["Entry", "Sequence"]], ecs


def search(cfg, queries: pd.DataFrame, ref: pd.DataFrame, tag: str) -> pd.DataFrame:
    h = cfg["homologs"]
    work = config.DATA / "homologs" / tag
    work.mkdir(parents=True, exist_ok=True)
    out = work / "hits.tsv"
    if not out.exists():
        for df, fa in ((ref, "ref.fasta"), (queries, "q.fasta")):
            with open(work / fa, "w") as f:
                for e, s in zip(df.Entry, df.Sequence):
                    f.write(f">{e}\n{s}\n")
        subprocess.run([str(data.mmseqs_bin(cfg)), "easy-search", work / "q.fasta", work / "ref.fasta",
                        out, work / "tmp", "-s", str(h["sensitivity"]), "--max-seqs", str(h["max_seqs"]),
                        "-e", str(h["evalue"]), "--format-output", FMT, "--threads", "4"],
                       check=True, stdout=subprocess.DEVNULL)
    hits = pd.read_csv(out, sep="\t", header=None, names=FMT.split(","))
    seq = dict(zip(queries.Entry, queries.Sequence))
    tseq = dict(zip(ref.Entry, ref.Sequence))
    keep = [q != t and seq[q] != tseq[t] for q, t in zip(hits["query"], hits.target)]
    return hits[keep]


def candidates(cfg, hits: pd.DataFrame, ref_ecs: pd.Series) -> pd.DataFrame:
    """One row per (query, exact EC) from the top-k hits, ordered NN-first."""
    k, max_c = cfg["homologs"]["top_k_hits"], cfg["homologs"]["max_candidates"]
    h = hits.sort_values(["query", "bits"], ascending=[True, False]).groupby("query").head(k).copy()
    h["rank"] = h.groupby("query").cumcount() + 1
    h["t_ecs"] = h.target.map(ref_ecs)
    h["w"] = h.bits / h.t_ecs.str.len()      # a multi-EC neighbour splits its vote
    per_q = h.groupby("query").agg(n_hits=("rank", "size"), tot_bits=("bits", "sum"))
    x = h.explode("t_ecs").rename(columns={"t_ecs": "ec"})
    first = x.sort_values("rank").groupby(["query", "ec"]).first()
    c = x.groupby(["query", "ec"]).agg(n_support=("rank", "size"), w_sum=("w", "sum"))
    c["best_rank"], c["best_fid"], c["best_qcov"], c["qlen"] = (
        first["rank"], first.fident, first.qcov, first.qlen)
    c = c.reset_index().join(per_q, on="query")
    c["bits_share"] = c.w_sum / c.tot_bits
    c = c.sort_values(["query", "best_rank", "bits_share"], ascending=[True, True, False])
    return c.groupby("query").head(max_c).reset_index(drop=True)


def render(query: str, cands: pd.DataFrame, names: dict) -> dict:
    rows = cands.to_dict("records")
    random.Random(int(hashlib.sha256(query.encode()).hexdigest()[:16], 16)).shuffle(rows)
    labels = list(string.ascii_uppercase.replace(NONE_LABEL, ""))[: len(rows)]
    n_hits = int(rows[0]["n_hits"])
    evidence, criteria, label_to_ec = [], {}, {}
    for lab, r in zip(labels, rows):
        name = names.get(r["ec"], "unnamed enzyme")
        evidence.append(
            f"Candidate {lab}: EC {r['ec']} ({name}). Best homolog: {r['best_fid'] * 100:.0f}% "
            f"identity, {r['best_qcov'] * 100:.0f}% query coverage, hit rank {int(r['best_rank'])} "
            f"of {n_hits}. Supported by {int(r['n_support'])} of {n_hits} homologs "
            f"({r['bits_share'] * 100:.0f}% of total alignment score).")
        criteria[lab] = f"EC {r['ec']} ({name})"
        label_to_ec[lab] = r["ec"]
    criteria[NONE_LABEL] = "none of the listed candidates is convincingly supported"
    state = (f"Query protein: {int(rows[0]['qlen'])} amino acids. {n_hits} homologs retrieved from "
             f"a reference set of annotated enzymes.\n" + "\n".join(evidence))
    return {"state": state, "questions": {"ec": {"type": "choice", "instructions": INSTRUCTIONS,
                                                 "criteria": criteria}},
            "label_to_ec": label_to_ec}
