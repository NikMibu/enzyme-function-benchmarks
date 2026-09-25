"""Embedding nearest neighbours (UniProt's ProtT5 embeddings) for the exact EC, ec1 + ec4.

    python scripts/08_embedding_nn.py

See the `embeddings` section of configs/default.yaml for methods and the pre-registered reading.
"""
import _common  # noqa: F401
import json

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from ed import config, homologs, metrics

cfg = config.load()
out = config.RESULTS / "embeddings"
out.mkdir(parents=True, exist_ok=True)
ec1 = pd.read_csv(config.BENCH / "ec1.tsv", sep="\t", dtype=str).assign(bench="ec1")
ec4 = pd.read_csv(config.BENCH / "ec4.tsv", sep="\t", dtype=str).assign(bench="ec4")
bench = pd.concat([ec1, ec4], ignore_index=True)
log = {"embedding_file_last_modified": (config.DATA / "raw" / "sprot_prott5_release.txt").read_text().strip()}

ref, ref_ecs = homologs.reference(cfg, bench)
first_ec = ref_ecs.map(lambda e: sorted(e)[0])

# ---------------------------------------------------------------- MMseqs2 nearest neighbour
tops = []
for tag, q in (("ec1", ec1), ("ec4", ec4)):
    h = homologs.search(cfg, q, ref, tag)
    tops.append(h.sort_values("bits", ascending=False).drop_duplicates("query"))
top = pd.concat(tops).set_index("query")
bench["top_identity"] = bench.Entry.map(top.fident)
bench["mmseqs_ec"] = bench.Entry.map(top.target).map(first_ec)

# ---------------------------------------------------------------- embeddings
with h5py.File(config.DATA / cfg["embeddings"]["file"], "r") as f:
    keys = set(f.keys())
    log["embeddings_in_file"] = len(keys)
    ref_in = ref[ref.Entry.isin(keys)].reset_index(drop=True)
    log["reference_n"], log["reference_with_embedding"] = len(ref), len(ref_in)
    missing_q = sorted(set(bench.Entry) - keys)
    log["benchmark_without_embedding"] = missing_q
    R = np.stack([f[e][()] for e in ref_in.Entry]).astype(np.float32)
    Q = np.stack([f[e][()] if e in keys else np.zeros(R.shape[1]) for e in bench.Entry]).astype(np.float32)
R /= np.linalg.norm(R, axis=1, keepdims=True)
Q /= np.maximum(np.linalg.norm(Q, axis=1, keepdims=True), 1e-9)
best_i, best_sim = np.zeros(len(Q), int), np.zeros(len(Q))
for s in range(0, len(R), 20000):                       # chunked to bound memory
    sim = Q @ R[s:s + 20000].T
    j = sim.argmax(1)
    better = sim[np.arange(len(Q)), j] > best_sim
    best_i[better], best_sim[better] = j[better] + s, sim[np.arange(len(Q)), j][better]
bench["emb_ec"] = ref_in.Entry.iloc[best_i].map(first_ec).to_numpy()
bench["emb_cosine"] = best_sim
bench.loc[bench.Entry.isin(missing_q), "emb_ec"] = np.nan
bench["hybrid_ec"] = bench.mmseqs_ec.fillna(bench.emb_ec)

# ---------------------------------------------------------------- scoring
bench["subset"] = np.select([bench.top_identity.isna(), bench.top_identity < 0.3],
                            ["no hit", "hit < 30%"], "hit >= 30%")
METHODS = {"mmseqs-nn": "mmseqs_ec", "emb-nn": "emb_ec", "hybrid": "hybrid_ec"}
for m, col in METHODS.items():
    bench[f"{m}_l4"] = bench[col].eq(bench.ec)
    bench[f"{m}_l1"] = bench[col].str.split(".").str[0].eq(bench.ec1)
bench.drop(columns=["Sequence", "Protein names"]).to_csv(out / "predictions.csv", index=False)


def acc(mask, col):
    k, n = int(bench.loc[mask, col].sum()), int(mask.sum())
    lo, hi = metrics.wilson(k, n) if n else (0, 0)
    return {"n": n, "correct": k, "accuracy": round(100 * k / n, 1) if n else None,
            "ci95": [round(100 * lo, 1), round(100 * hi, 1)]}


table = {}
for bname, bmask in [("ec1+ec4", bench.bench.notna()), ("ec1", bench.bench.eq("ec1")), ("ec4", bench.bench.eq("ec4"))]:
    for sub in ["all", "no hit", "hit < 30%", "hit >= 30%"]:
        mask = bmask & (bench.subset.eq(sub) if sub != "all" else True)
        for lvl in ["l4", "l1"]:
            table[(bname, sub, lvl)] = {m: acc(mask, f"{m}_{lvl}") for m in METHODS}


def paired(a, b, lvl, mask):
    ca, cb = bench.loc[mask, f"{a}_{lvl}"].to_numpy(), bench.loc[mask, f"{b}_{lvl}"].to_numpy()
    w, l = int((ca & ~cb).sum()), int((~ca & cb).sum())
    p = stats.binomtest(w, w + l, 0.5).pvalue if w + l else 1.0
    return {"a": a, "b": b, "level": lvl, "n": int(mask.sum()), "acc_a": round(100 * ca.mean(), 1),
            "acc_b": round(100 * cb.mean(), 1), "a_only_correct": w, "b_only_correct": l,
            "mcnemar_p": float(f"{p:.3g}")}


allm = bench.bench.notna()
tests = [paired("hybrid", "mmseqs-nn", "l4", allm), paired("emb-nn", "mmseqs-nn", "l4", allm),
         paired("hybrid", "mmseqs-nn", "l1", allm), paired("emb-nn", "mmseqs-nn", "l1", allm)]
log["subset_sizes"] = bench.groupby(["bench", "subset"]).size().unstack(fill_value=0).to_dict("index")
log["emb_cosine_by_subset"] = bench.groupby("subset").emb_cosine.median().round(3).to_dict()
report = {**log, "paired_tests": tests,
          "pre_registered_reading": {
              "hybrid_beats_mmseqs_exact_ec": bool(tests[0]["mcnemar_p"] < 0.05 and tests[0]["acc_a"] > tests[0]["acc_b"]),
              "emb_nn_beats_mmseqs_exact_ec": bool(tests[1]["mcnemar_p"] < 0.05 and tests[1]["acc_a"] > tests[1]["acc_b"]),
          },
          "accuracy": {f"{b} | {s} | level {lvl[1]}": v for (b, s, lvl), v in table.items()}}
json.dump(report, open(out / "report.json", "w"), indent=2)
rows = [{"benchmark": b, "subset": s, "level": lvl[1], **{m: v[m]["accuracy"] for m in METHODS},
         "n": v["mmseqs-nn"]["n"]} for (b, s, lvl), v in table.items()]
pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)
print(json.dumps({k: v for k, v in report.items() if k != "accuracy"}, indent=1))
print(pd.DataFrame(rows).to_string(index=False))
