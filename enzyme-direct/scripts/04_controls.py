"""Sequence-based controls for ec1, scored on the same 210 proteins as the models.

    python scripts/04_controls.py

logreg-ctx1 / logreg-ctx2: logistic regression on exactly the information Jev gets at context
    level 1 / 2, trained on pool proteins that share no 30%-identity family with ec1 and are not
    in ec4 (so no benchmark label is ever seen in training).
nn-care-train: EC class of the best MMseqs2 hit in CARE's training set (Swiss-Prot enzymes) with
    the benchmark proteins removed; the practical ceiling for a sequence tool.
"""
import _common  # noqa: F401
import json
import subprocess

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ed import config, data, features, metrics, prompts

cfg = config.load()
classes = list(prompts.EC1_CLASSES)
out = config.RESULTS / "ec1"
out.mkdir(parents=True, exist_ok=True)
ec1 = pd.read_csv(config.BENCH / "ec1.tsv", sep="\t", dtype=str)
ec4 = pd.read_csv(config.BENCH / "ec4.tsv", sep="\t", dtype=str)
log = {}


def save(name: str, P: np.ndarray, extra: dict, abstain=None) -> None:
    """abstain: boolean mask of proteins with no prediction (scored as wrong, uniform probabilities)."""
    labels = [classes[i] for i in P.argmax(1)]
    if abstain is not None:
        labels = ["none" if a else lab for a, lab in zip(abstain, labels)]
    pred = pd.DataFrame({"entry": ec1.Entry, "true": ec1.ec1, "pred": labels,
                         **{f"p_{c}": P[:, i] for i, c in enumerate(classes)}, **extra})
    pred.to_csv(out / f"{name}_predictions.csv", index=False)
    s = metrics.summarize(pred, classes, cfg["gate"])
    s["model"] = name
    json.dump(s, open(out / f"{name}_summary.json", "w"), indent=2)
    metrics.confusion(pred, classes).to_csv(out / f"{name}_confusion.csv")
    print(f"{name}: {s['accuracy']}% {s['accuracy_ci95']}, gate {s['passes_gate']}")


# ---- training pool for the logistic-regression controls
care_acc, care_seq = data.care_exclusion(cfg)
raw, _ = data.fetch_uniprot(cfg)
pool, _ = data.clean_pool(raw, care_acc, care_seq)
cl30 = data.cluster(cfg, pool, cfg["mmseqs"]["cluster_min_seq_id"], "pool30")
test_families = set(pool.Entry[pool.Entry.isin(ec1.Entry)].map(cl30))
train = pool[pool.ec1.isin(classes) & ~pool.Entry.isin(ec1.Entry) & ~pool.Entry.isin(ec4.Entry)
             & ~pool.Entry.map(cl30).isin(test_families)].reset_index(drop=True)
log["logreg_train_n"] = len(train)
log["logreg_train_by_class"] = train.ec1.value_counts().sort_index().to_dict()

for level in cfg["context"]["levels"]:
    Xtr, names = features.matrix(train.Sequence, level)
    Xte, _ = features.matrix(ec1.Sequence, level)
    lr = make_pipeline(StandardScaler(), LogisticRegression(
        C=cfg["controls"]["logreg"]["C"], max_iter=cfg["controls"]["logreg"]["max_iter"],
        class_weight="balanced"))
    lr.fit(Xtr, train.ec1)
    P = lr.predict_proba(Xte)[:, [list(lr.classes_).index(c) for c in classes]]
    log[f"logreg_ctx{level}_features"] = names
    save(f"logreg-ctx{level}", P, {})

# ---- nearest neighbour against CARE's training set (benchmark proteins removed)
ref_path = data.fetch_care(cfg) / cfg["controls"]["nn_reference"]
if not ref_path.exists():
    import urllib.request
    urllib.request.urlretrieve(f"{cfg['care']['raw_base']}/{cfg['care']['commit']}/"
                               f"{cfg['controls']['nn_reference']}", ref_path)  # nosec - pinned
ref = pd.read_csv(ref_path, usecols=["Entry", "Sequence", "EC number"])
bench_all = pd.concat([ec1, ec4])
ref = ref[~ref.Entry.isin(bench_all.Entry) & ~ref.Sequence.isin(bench_all.Sequence)]
ref_ec1 = ref.groupby("Entry")["EC number"].agg(lambda s: sorted({e.split(".")[0] for e in s}))
ref = ref.drop_duplicates("Entry")
log["nn_reference_n"] = len(ref)
work = config.DATA / "nn"
work.mkdir(parents=True, exist_ok=True)
for df, fa in ((ref, "ref.fasta"), (ec1, "q.fasta")):
    with open(work / fa, "w") as f:
        for e, s in zip(df.Entry, df.Sequence):
            f.write(f">{e}\n{s}\n")
fmt = "query,target,fident,qcov,evalue,bits"
subprocess.run([str(data.mmseqs_bin(cfg)), "easy-search", work / "q.fasta", work / "ref.fasta",
                work / "hits.tsv", work / "tmp", "-s", str(cfg["controls"]["nn_sensitivity"]),
                "--format-output", fmt, "--threads", "4"], check=True, stdout=subprocess.DEVNULL)
hits = pd.read_csv(work / "hits.tsv", sep="\t", header=None, names=fmt.split(","))
top = hits.sort_values(["query", "bits"], ascending=[True, False]).drop_duplicates("query").set_index("query")
P = np.full((len(ec1), len(classes)), 1 / len(classes))   # no hit: no prediction, uniform probabilities
fid = np.full(len(ec1), np.nan)
for i, e in enumerate(ec1.Entry):
    if e in top.index:
        hit_classes = [c for c in ref_ec1[top.at[e, "target"]] if c in classes]
        if hit_classes:
            P[i] = 0
            P[i, [classes.index(c) for c in hit_classes]] = 1 / len(hit_classes)
            fid[i] = top.at[e, "fident"]
log["nn_with_hit"] = int(np.isfinite(fid).sum())
log["nn_top_identity_quartiles"] = [round(float(q), 3) for q in np.nanpercentile(fid, [25, 50, 75])]
save("nn-care-train", P, {"top_hit_identity": fid}, abstain=~np.isfinite(fid))
json.dump(log, open(out / "controls_log.json", "w"), indent=2)
print(json.dumps({k: v for k, v in log.items() if "features" not in k}, indent=1))
