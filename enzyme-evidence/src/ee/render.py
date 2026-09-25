"""Candidate table -> compact text evidence for decision models (Jev now, Laya in phase 2).

Deliberately excluded: the query accession / entry name (entry names such as LDH_BACSU encode
function, and hosted models may have memorised accession->EC), sequences, neighbour
accessions, and e-values in scientific notation. Candidates are shuffled per query with a
deterministic seed and given opaque letters, so position and label wording carry no signal.
"""
import hashlib
import random
import string

import pandas as pd

NONE_LABEL = "N"
INSTRUCTIONS = (
    "Which candidate enzyme function (EC number) is best supported by the homologs retrieved "
    "for this protein? Support is stronger with higher sequence identity, higher alignment "
    "coverage, a better (lower) hit rank, and more agreeing homologs. Choose N if no candidate "
    "is convincingly supported."
)


def _rng(query: str) -> random.Random:
    return random.Random(int(hashlib.sha256(query.encode()).hexdigest()[:16], 16))


def render(query: str, cands: pd.DataFrame, ec_names: dict) -> dict:
    """Returns {'state', 'questions', 'label_to_ec'} for one query's candidates."""
    rows = cands.to_dict("records")
    _rng(query).shuffle(rows)
    labels = list(string.ascii_uppercase.replace(NONE_LABEL, ""))[: len(rows)]
    n_hits = int(rows[0]["n_hits"]) if rows else 0
    evidence, criteria, label_to_ec = [], {}, {}
    for lab, r in zip(labels, rows):
        name = ec_names.get(r["ec"], "unnamed enzyme")
        evidence.append(
            f"Candidate {lab}: EC {r['ec']} ({name}). Best homolog: {r['best_fid'] * 100:.0f}% "
            f"identity, {r['best_qcov'] * 100:.0f}% query coverage, hit rank {int(r['best_rank'])} "
            f"of {n_hits}. Supported by {int(r['n_support'])} of {n_hits} homologs "
            f"({r['bits_share'] * 100:.0f}% of total alignment score)."
        )
        criteria[lab] = f"EC {r['ec']} ({name})"
        label_to_ec[lab] = r["ec"]
    criteria[NONE_LABEL] = "none of the listed candidates is convincingly supported"
    qlen = int(rows[0]["qlen"]) if rows else 0
    state = (f"Query protein: {qlen} amino acids. {n_hits} homologs retrieved from a reference "
             f"set of annotated enzymes.\n" + "\n".join(evidence))
    questions = {"ec": {"type": "choice", "instructions": INSTRUCTIONS, "criteria": criteria}}
    return {"state": state, "questions": questions, "label_to_ec": label_to_ec}
