"""Single evaluation on the CARE test splits with the frozen config."""
import _common  # noqa: F401
import json
import subprocess

import pandas as pd

from ee import config, evaluate, pipeline

cfg = config.load()
if not cfg.get("frozen"):
    raise SystemExit("Config is not frozen. Run scripts/03_tune_dev.py first.")

queries, train, masked = pipeline.load_inputs(cfg)
test_q = queries[queries.role == "test"]
k = cfg["candidates"]["top_k_hits"]
c, scores, _ = pipeline.score_all(cfg, queries, train, masked, k, "test")

out = config.RESULTS / "test"
(out / "predictions").mkdir(parents=True, exist_ok=True)
preds = {name: evaluate.predictions(test_q, c, s) for name, s in scores.items()}

summary, paired, bins = [], {}, []
for name, p in preds.items():
    evaluate.to_care_csv(p, out / "predictions" / f"{name}.csv")
    p.drop(columns="ranked").to_csv(out / "predictions" / f"{name}_per_query.csv", index=False)
    for split, ps in list(p.groupby("split")) + [("pooled", p)]:
        summary.append({"model": name, "split": split, **evaluate.summarize(ps)})
        if name not in ("nearest_neighbour", "oracle"):
            ref = preds["nearest_neighbour"]
            paired[f"{name}|{split}"] = evaluate.paired(ps, ref[ref["query"].isin(ps["query"])])
    b = evaluate.by_identity(p).reset_index().assign(model=name)
    bins.append(b)

# Strict subset of the <30% split: no training hit >30% identity over >80% query coverage.
h30 = masked[masked["query"].isin(test_q[test_q.split == "30"].Entry)]
leaky = set(h30[(h30.fident > 0.3) & (h30.qcov > 0.8)]["query"])
for name, p in preds.items():
    ps = p[(p.split == "30") & ~p["query"].isin(leaky)]
    summary.append({"model": name, "split": "30-strict", **evaluate.summarize(ps)})

summary = pd.DataFrame(summary)
summary.to_csv(out / "summary.csv", index=False)
pd.concat(bins).to_csv(out / "by_identity.csv", index=False)
commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
json.dump({"k": k, "code_commit_at_run": commit, "care_commit": cfg["care"]["commit"],
           "n_leaky_30": len(leaky), "paired_vs_nearest_neighbour": paired},
          open(out / "run.json", "w"), indent=2)
cols = ["model", "split", "n", "L4", "L3", "L1", "oracle_L4", "contested_L4", "ECE", "AURC_risk", "sel_acc@50%"]
print(summary[cols].to_string(index=False))
print(json.dumps(paired, indent=1))
