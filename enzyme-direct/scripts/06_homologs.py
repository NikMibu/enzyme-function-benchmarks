"""The enzyme-evidence setup on the ec1 benchmark: Jev chooses among homolog-derived EC candidates.

    python scripts/06_homologs.py --pilot 3     # plumbing check on pool proteins outside the benchmarks
    python scripts/06_homologs.py               # full run (~200 Jev calls) + evaluation

Methods scored on the same 210 proteins, at level 1 (EC class) and level 4 (exact EC):
    nn-homologs     EC of the best hit (first candidate)
    vote-homologs   candidate with the largest share of alignment score among the top-k hits
    jev-homologs    Jev's choice among the candidates ('none' mass is not given to any candidate)
    jev-hybrid      jev-homologs, and Jev's level-2 (motif) answer for proteins without a hit
Proteins without a hit get no prediction from a hit-based method and are scored wrong.
"""
import _common  # noqa: F401
import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats

from ed import config, data, homologs, metrics, models, prompts

ap = argparse.ArgumentParser()
ap.add_argument("--pilot", type=int)
args = ap.parse_args()
cfg = config.load()
classes = list(prompts.EC1_CLASSES)
ec1 = pd.read_csv(config.BENCH / "ec1.tsv", sep="\t", dtype=str)
ec4 = pd.read_csv(config.BENCH / "ec4.tsv", sep="\t", dtype=str)
names, _ = data.ec_names(cfg)

if args.pilot:
    pool = pd.read_csv(config.DATA / "raw" / "uniprot_pool.tsv", sep="\t", dtype=str)
    queries = pool[~pool.Entry.isin(pd.concat([ec1, ec4]).Entry)].sample(args.pilot, random_state=0)
    tag, cache = "pilot", config.DATA / "pilot" / "ec1_jev-homologs.jsonl"
else:
    queries, tag, cache = ec1, "ec1", config.RESULTS / "ec1" / "jev-homologs_responses.jsonl"

ref, ref_ecs = homologs.reference(cfg, pd.concat([ec1, ec4, queries]))
hits = homologs.search(cfg, queries, ref, tag)
cands = homologs.candidates(cfg, hits, ref_ecs)
rendered = {q: homologs.render(q, g, names) for q, g in cands.groupby("query")}
items = [(q, r["state"], r["questions"]) for q, r in rendered.items()]
resp = models.run("jev", cfg, items, cache)
if args.pilot:
    for q, r in rendered.items():
        print(r["state"][:600], "\n->", json.dumps(resp[q]["answers"])[:300], "\n")
    raise SystemExit

# ---------------------------------------------------------------- per-protein predictions
out = config.RESULTS / "ec1"
by_q = {q: g for q, g in cands.groupby("query")}
ctx2 = pd.read_csv(out / "jev-ctx2_predictions.csv", dtype=str).set_index("entry")


def jev_probs(q):
    r = rendered[q]
    p = resp[q]["answers"]["ec"]["probabilities"]
    ec_p = {ec: float(p.get(lab, 0.0)) for lab, ec in r["label_to_ec"].items()}
    return ec_p, float(p.get(homologs.NONE_LABEL, 0.0))


def pick(method, q):
    """-> (exact EC or None, {class: prob})"""
    if q not in by_q:
        if method == "jev-hybrid":
            row = ctx2.loc[q]
            return None, {c: float(row[f"p_{c}"]) for c in classes}
        return None, None
    g = by_q[q]
    if method == "nn-homologs":
        ec = g.ec.iloc[0]
        return ec, {ec.split(".")[0]: 1.0}
    if method == "vote-homologs":
        ec = g.sort_values("bits_share", ascending=False).ec.iloc[0]
        return ec, {ec.split(".")[0]: 1.0}
    ec_p, _ = jev_probs(q)
    ec = max(ec_p, key=ec_p.get)
    cls = {}
    for e, v in ec_p.items():
        cls[e.split(".")[0]] = cls.get(e.split(".")[0], 0.0) + v
    return ec, cls


METHODS = ["nn-homologs", "vote-homologs", "jev-homologs", "jev-hybrid"]
summary, tables = {}, {}
for m in METHODS:
    recs = []
    for q, true_ec, true_c in zip(ec1.Entry, ec1.ec, ec1.ec1):
        ec, cls = pick(m, q)
        P = np.array([(cls or {}).get(c, 0.0) for c in classes])
        P = P if P.sum() > 0 else np.full(len(classes), 1 / len(classes))
        pred_c = classes[int(P.argmax())] if cls else "none"
        rec = {"entry": q, "true": true_c, "pred": pred_c, "true_ec": true_ec, "pred_ec": ec or "none",
               "has_hit": q in by_q, **{f"p_{c}": P[i] for i, c in enumerate(classes)}}
        if m.startswith("jev") and q in by_q:
            rec["p_none_option"] = jev_probs(q)[1]
        recs.append(rec)
    pred = pd.DataFrame(recs)
    pred.to_csv(out / f"{m}_predictions.csv", index=False)
    s = metrics.summarize(pred, classes, cfg["gate"])
    k4 = int((pred.pred_ec == pred.true_ec).sum())
    lo, hi = metrics.wilson(k4, len(pred))
    s["model"] = m
    s["level4_accuracy"] = round(100 * k4 / len(pred), 1)
    s["level4_ci95"] = [round(100 * lo, 1), round(100 * hi, 1)]
    json.dump(s, open(out / f"{m}_summary.json", "w"), indent=2)
    metrics.confusion(pred, classes).to_csv(out / f"{m}_confusion.csv")
    summary[m] = s
    tables[m] = pred

# oracle: the true EC / class is among the candidates
has = ec1.Entry.isin(by_q)
orc4 = np.mean([q in by_q and t in set(by_q[q].ec) for q, t in zip(ec1.Entry, ec1.ec)])
orc1 = np.mean([q in by_q and c in {e.split(".")[0] for e in by_q[q].ec} for q, c in zip(ec1.Entry, ec1.ec1)])


def paired(a, b, col_true, col_pred):
    ca = (tables[a][col_true] == tables[a][col_pred]).to_numpy()
    cb = (tables[b][col_true] == tables[b][col_pred]).to_numpy()
    w, l = int((ca & ~cb).sum()), int((~ca & cb).sum())
    p = stats.binomtest(w, w + l, 0.5).pvalue if w + l else 1.0
    return {"a": a, "b": b, "delta": round(100 * (ca.mean() - cb.mean()), 1),
            "a_only_correct": w, "b_only_correct": l, "mcnemar_p": float(f"{p:.3g}")}


tests = {lvl: [paired(a, b, t, pcol) for a, b in [("jev-homologs", "nn-homologs"),
                                                  ("jev-homologs", "vote-homologs"),
                                                  ("jev-hybrid", "nn-homologs")]]
         for lvl, t, pcol in [("level1", "true", "pred"), ("level4", "true_ec", "pred_ec")]}
t1 = tests["level1"][0]
report = {
    "n": len(ec1), "with_hit": int(has.sum()),
    "median_candidates": float(cands.groupby("query").size().median()),
    "median_top_identity": round(float(hits.sort_values("bits", ascending=False)
                                       .drop_duplicates("query").fident.median()), 3),
    "unnamed_candidate_ecs": int(sum(ec not in names for ec in cands.ec.unique())),
    "oracle_level1": round(100 * orc1, 1), "oracle_level4": round(100 * orc4, 1),
    "accuracy": {m: {"level1": summary[m]["accuracy"], "level1_ci95": summary[m]["accuracy_ci95"],
                     "level4": summary[m]["level4_accuracy"], "level4_ci95": summary[m]["level4_ci95"],
                     "log_loss_level1": summary[m]["log_loss"]} for m in METHODS},
    "paired_tests": tests,
    "pre_registered_reading": {"jev_adds_to_homolog_evidence": bool(t1["mcnemar_p"] < 0.05 and t1["delta"] > 0)},
}
json.dump(report, open(out / "homologs_report.json", "w"), indent=2)
print(json.dumps(report, indent=1))
