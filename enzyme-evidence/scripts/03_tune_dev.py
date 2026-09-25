"""Choose top-k on the dev set (CARE-protocol replica), run ablations, then freeze the config.

Nothing here reads CARE test labels: test queries are present in the candidate table but
are filtered out before scoring.
"""
import _common  # noqa: F401
import json

import pandas as pd

from ee import candidates, config, evaluate, pipeline

cfg = config.load()
queries, train, masked = pipeline.load_inputs(cfg)
dev_q = queries[queries.role == "dev"]
out_dir = config.RESULTS / "dev"
out_dir.mkdir(parents=True, exist_ok=True)

rows = []
for k in cfg["candidates"]["k_grid"]:
    c, scores, _ = pipeline.score_all(cfg, queries, train, masked, k, "dev")
    for name, s in scores.items():
        p = evaluate.predictions(dev_q, c, s)
        for split, ps in [("pooled", p)] + list(p.groupby("split")):
            rows.append({"k": k, "model": name, "split": split, **evaluate.summarize(ps)})
grid = pd.DataFrame(rows)
grid.to_csv(out_dir / "k_grid.csv", index=False)

pooled = grid[(grid.split == "pooled") & (grid.model == "lightgbm")].set_index("k")
best_k = int(pooled.L4.idxmax())
print(grid[grid.split == "pooled"].pivot(index="k", columns="model", values="L4").to_string())
print("chosen k (dev LightGBM L4):", best_k)

# Ablation on dev: does the EC-frequency prior matter? (CARE tests are EC3-balanced, not natural.)
no_freq = [f for f in candidates.FEATURES if f != "log_train_freq"]
c, scores, model = pipeline.score_all(cfg, queries, train, masked, best_k, "dev", features=no_freq)
abl = evaluate.summarize(evaluate.predictions(dev_q, c, scores["lightgbm"]))
c, scores, model = pipeline.score_all(cfg, queries, train, masked, best_k, "dev")
full = evaluate.summarize(evaluate.predictions(dev_q, c, scores["lightgbm"]))
imp = dict(sorted(zip(candidates.FEATURES, model.booster_.feature_importance("gain").round(1)),
                  key=lambda x: -x[1]))
json.dump({"chosen_k": best_k, "lightgbm": full, "lightgbm_no_freq_prior": abl,
           "feature_gain": {k: float(v) for k, v in imp.items()}},
          open(out_dir / "summary.json", "w"), indent=2)
print("LightGBM dev L4 with / without EC-frequency prior:", full["L4"], abl["L4"])

cfg["candidates"]["top_k_hits"] = best_k
cfg["frozen"] = True
config.save(cfg)
print("config frozen")
