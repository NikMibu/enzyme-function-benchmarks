"""ESM-2 as the sequence-only reference on ec1, and the hard subset where homology search is weak.

    python scripts/07_esm.py        # ~1 h on 4 CPU cores (embeddings are cached)

See the `esm` section of configs/default.yaml for the methods and the pre-registered reading.
"""
import _common  # noqa: F401
import json

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ed import config, data, features, homologs, metrics, prompts
from ed.esm import Embedder

cfg = config.load()
E = cfg["esm"]
classes = list(prompts.EC1_CLASSES)
out = config.RESULTS / "ec1"
ec1 = pd.read_csv(config.BENCH / "ec1.tsv", sep="\t", dtype=str)
ec4 = pd.read_csv(config.BENCH / "ec4.tsv", sep="\t", dtype=str)
log = {}

# ---------------------------------------------------------------- training set: no ec1 homologs
ref, ref_ecs = homologs.reference(cfg, pd.concat([ec1, ec4]))
hits = homologs.search(cfg, ec1, ref, "ec1")
homolog_targets = set(hits.target)
u = cfg["uniprot"]
ref = ref.assign(ecs=ref.Entry.map(ref_ecs))
ok = (ref.ecs.str.len().eq(1)
      & ref.ecs.str[0].str.fullmatch(r"\d+\.\d+\.\d+\.\d+").astype(bool)
      & ref.Sequence.str.len().between(u["min_len"], u["max_len"])
      & ref.Sequence.str.fullmatch(data.STANDARD_AA.pattern).astype(bool)
      & ~ref.Entry.isin(homolog_targets))
cand = ref[ok].assign(ec1=lambda d: d.ecs.str[0].str.split(".").str[0])
cand = cand[cand.ec1.isin(classes)]
log["train_pool_without_ec1_homologs"] = cand.ec1.value_counts().sort_index().to_dict()
train = (cand.groupby("ec1").sample(E["per_class"], random_state=cfg["seed"])
         .sort_values(["ec1", "Entry"]).reset_index(drop=True))
log["train_n"] = len(train)
log["excluded_homologs_of_ec1"] = len(homolog_targets)

# ---------------------------------------------------------------- embeddings
emb = Embedder(E["model"])
Xte = emb(ec1.Entry, ec1.Sequence)
Xtr = emb(train.Entry, train.Sequence)


def save(name, P, abstain=None):
    labels = [classes[i] for i in P.argmax(1)]
    if abstain is not None:
        labels = ["none" if a else lab for a, lab in zip(abstain, labels)]
    pred = pd.DataFrame({"entry": ec1.Entry, "true": ec1.ec1, "pred": labels,
                         **{f"p_{c}": P[:, i] for i, c in enumerate(classes)}})
    pred.to_csv(out / f"{name}_predictions.csv", index=False)
    s = metrics.summarize(pred, classes, cfg["gate"])
    s["model"] = name
    json.dump(s, open(out / f"{name}_summary.json", "w"), indent=2)
    metrics.confusion(pred, classes).to_csv(out / f"{name}_confusion.csv")
    print(f"{name}: {s['accuracy']}% {s['accuracy_ci95']}, gate {s['passes_gate']}")
    return pred


def logreg(Xa, ya, Xb):
    m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced"))
    m.fit(Xa, ya)
    return m.predict_proba(Xb)[:, [list(m.classes_).index(c) for c in classes]]


preds = {"esm-probe": save("esm-probe", logreg(Xtr, train.ec1, Xte))}
A = Xtr / np.linalg.norm(Xtr, axis=1, keepdims=True)
B = Xte / np.linalg.norm(Xte, axis=1, keepdims=True)
sim = B @ A.T
nn_idx = sim.argmax(1)
P = np.zeros((len(ec1), len(classes)))
P[np.arange(len(ec1)), [classes.index(c) for c in train.ec1.iloc[nn_idx]]] = 1
preds["esm-knn"] = save("esm-knn", P)
log["esm_knn_top_cosine_quartiles"] = [round(float(q), 3) for q in np.percentile(sim.max(1), [25, 50, 75])]
F_tr, _ = features.matrix(train.Sequence, 1)
F_te, _ = features.matrix(ec1.Sequence, 1)
preds["logreg-ctx1-care"] = save("logreg-ctx1-care", logreg(F_tr, train.ec1, F_te))

# ---------------------------------------------------------------- hard subset and hybrid
top = hits.sort_values("bits", ascending=False).drop_duplicates("query").set_index("query").fident
top_id = ec1.Entry.map(top)
hard = (top_id.isna() | (top_id < E["hard_identity"])).to_numpy()
log["hard_n"], log["hard_no_hit"] = int(hard.sum()), int(top_id.isna().sum())
for r in ["nn-homologs", "jev-homologs", "jev-ctx2", "logreg-ctx2", "jev", "laya-english"]:
    preds[r] = pd.read_csv(out / f"{r}_predictions.csv", dtype=str)
nn = preds["nn-homologs"].set_index("entry").loc[ec1.Entry]
Pp = preds["esm-probe"][[f"p_{c}" for c in classes]].to_numpy(float)
Ph = np.where(hard[:, None], Pp, nn[[f"p_{c}" for c in classes]].to_numpy(float))
hyb = save("hybrid-nn-esm", Ph)
preds["hybrid-nn-esm"] = hyb


def correct(r):
    p = preds[r].set_index("entry").loc[ec1.Entry]
    return (p["true"].astype(str) == p["pred"].astype(str)).to_numpy()


def paired(a, b, mask):
    ca, cb = correct(a)[mask], correct(b)[mask]
    w, l = int((ca & ~cb).sum()), int((~ca & cb).sum())
    p = stats.binomtest(w, w + l, 0.5).pvalue if w + l else 1.0
    return {"a": a, "b": b, "n": int(mask.sum()), "acc_a": round(100 * ca.mean(), 1),
            "acc_b": round(100 * cb.mean(), 1), "a_only_correct": w, "b_only_correct": l,
            "mcnemar_p": float(f"{p:.3g}")}


full = np.ones(len(ec1), bool)
by_subset = {}
for name, mask in [("all", full), ("hard", hard), ("easy", ~hard)]:
    by_subset[name] = {r: {"n": int(mask.sum()), "accuracy": round(100 * correct(r)[mask].mean(), 1),
                           "ci95": [round(100 * v, 1) for v in metrics.wilson(int(correct(r)[mask].sum()), int(mask.sum()))]}
                       for r in preds}
tests = [paired("esm-probe", "nn-homologs", hard), paired("hybrid-nn-esm", "nn-homologs", full),
         paired("esm-probe", "logreg-ctx1-care", full), paired("esm-probe", "jev-ctx2", full),
         paired("esm-knn", "nn-homologs", hard)]
s_probe = json.load(open(out / "esm-probe_summary.json"))
report = {**log, "accuracy_by_subset": by_subset, "paired_tests": tests,
          "pre_registered_reading": {
              "esm_reads_function_from_sequence": s_probe["passes_gate"],
              "esm_helps_where_homology_is_weak": bool(tests[0]["mcnemar_p"] < 0.05
                                                       and tests[0]["acc_a"] > tests[0]["acc_b"]),
              "hybrid_beats_nn": bool(tests[1]["mcnemar_p"] < 0.05 and tests[1]["acc_a"] > tests[1]["acc_b"]),
          }}
json.dump(report, open(out / "esm_report.json", "w"), indent=2)
print(json.dumps({k: v for k, v in report.items() if k != "accuracy_by_subset"}, indent=1))
print(pd.DataFrame({k: {r: v["accuracy"] for r, v in d.items()} for k, d in by_subset.items()}).to_string())
