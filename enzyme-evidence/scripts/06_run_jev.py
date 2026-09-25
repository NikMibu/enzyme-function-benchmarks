"""Jev zero-shot integrator on the same candidates. Needs TYPESAFE_API_KEY.

    python scripts/06_run_jev.py --role dev      # prompt checks, ~$0.03
    python scripts/06_run_jev.py --role test     # once, after the config is frozen
"""
import _common  # noqa: F401
import argparse
import json

from ee import care, config, evaluate, pipeline
from ee.models import jev

ap = argparse.ArgumentParser()
ap.add_argument("--role", choices=["dev", "test"], required=True)
args = ap.parse_args()

cfg = config.load()
if args.role == "test" and not cfg.get("frozen"):
    raise SystemExit("Config is not frozen. Run scripts/03_tune_dev.py first.")
queries, train, masked = pipeline.load_inputs(cfg)
q = queries[queries.role == args.role]
c, scores, _ = pipeline.score_all(cfg, queries, train, masked, cfg["candidates"]["top_k_hits"], args.role)

out = config.RESULTS / args.role
(out / "predictions").mkdir(parents=True, exist_ok=True)
s = jev.score(cfg, c, care.load_ec_names(cfg), out / "jev_responses.jsonl")
p = evaluate.predictions(q, c, s)
nn = evaluate.predictions(q, c, scores["nearest_neighbour"])
evaluate.to_care_csv(p, out / "predictions" / "jev.csv")
p.drop(columns="ranked").to_csv(out / "predictions" / "jev_per_query.csv", index=False)
res = {split: {**evaluate.summarize(ps), "vs_nn": evaluate.paired(ps, nn[nn["query"].isin(ps["query"])])}
       for split, ps in list(p.groupby("split")) + [("pooled", p)]}
json.dump(res, open(out / "jev_summary.json", "w"), indent=2)
print(json.dumps(res, indent=1))
