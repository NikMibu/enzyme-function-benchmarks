"""Post-hoc (added after the test run; no model changes): is the learned confidence better than
*any* simple scalar confidence for the nearest-neighbour prediction? Compares risk-coverage
area (rank-based, so scale-free) against identity, bitscore, identity x coverage, and vote share.
"""
import _common  # noqa: F401
import json

import numpy as np
import pandas as pd

from ee import candidates, config, evaluate, pipeline

cfg = config.load()
out = config.RESULTS / "test"
queries, train, masked = pipeline.load_inputs(cfg)
c = candidates.build(masked, queries, train, cfg["candidates"]["top_k_hits"],
                     cfg["candidates"]["max_candidates"])
c = c[c.role == "test"]
top_bits = c.groupby("query").best_bits.max()

nn = pd.read_csv(out / "predictions" / "nearest_neighbour_per_query.csv").set_index("query")
lg = pd.read_csv(out / "predictions" / "lightgbm_per_query.csv").set_index("query").reindex(nn.index)
vote_share = c.set_index(["query", "ec"]).bits_share
nn_share = pd.Series([vote_share.get((q, e), 0.0) for q, e in zip(nn.index, nn.pred)], index=nn.index)

confs = {
    "NN | identity": nn.top_fid.fillna(0),
    "NN | bitscore": nn.index.map(top_bits).fillna(0).to_series(index=nn.index),
    "NN | identity x coverage": nn.conf,
    "NN | vote share of NN EC": nn_share,
}
res = {}
for split in ["30", "30-50", "price", "all"]:
    m = (nn.split.astype(str) == split) if split != "all" else np.ones(len(nn), bool)
    ok_nn, ok_lg = nn.correct[m].to_numpy(float), lg.correct[m].to_numpy(float)
    row = {"LightGBM | P(correct)": round(evaluate.aurc(lg.conf[m].to_numpy(float), ok_lg), 3)}
    for name, s in confs.items():
        row[name] = round(evaluate.aurc(s[m].to_numpy(float), ok_nn), 3)
        row[f"LightGBM vs {name}"] = evaluate.paired_aurc(
            lg.conf[m].to_numpy(float), ok_lg, s[m].to_numpy(float), ok_nn)
    res[split] = row
json.dump(res, open(out / "confidence_analysis.json", "w"), indent=2)
for split, row in res.items():
    print(f"\n[{split}] AURC (lower = better)")
    for k, v in row.items():
        print(f"  {k}: {v}")
