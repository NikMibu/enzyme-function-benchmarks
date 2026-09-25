"""Query one model on one benchmark (sequence only).

    python scripts/02_run.py --stage ec1 --model jev
    python scripts/02_run.py --stage ec1 --model laya-english
    python scripts/02_run.py --stage ec4 --model jev     # only if that model passed the ec1 gate
"""
import _common  # noqa: F401
import argparse
import json

import pandas as pd

from ed import config, models, prompts

ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["ec1", "ec4"], required=True)
ap.add_argument("--model", choices=["jev", "laya-english", "laya-multilingual"], required=True)
ap.add_argument("--pilot", type=int, help="plumbing check on N pool proteins outside every benchmark")
args = ap.parse_args()

cfg = config.load()
bench = pd.read_csv(config.BENCH / f"{args.stage}.tsv", sep="\t", dtype=str)
if args.stage == "ec4":
    gate = config.RESULTS / "ec1" / f"{args.model}_summary.json"
    if not gate.exists() or not json.load(open(gate))["passes_gate"]:
        raise SystemExit(f"{args.model} did not pass the ec1 gate; ec4 is not run for it.")
    labels = pd.read_csv(config.BENCH / "ec4_labels.tsv", sep="\t", dtype=str)
    qs = prompts.ec4_question(list(labels.ec), dict(zip(labels.ec, labels.name)))
else:
    qs = prompts.ec1_question()
cache = config.RESULTS / args.stage / f"{args.model}_responses.jsonl"
if args.pilot:
    used = set(pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in config.BENCH.glob("ec*.tsv")]).Entry)
    pool = pd.read_csv(config.DATA / "raw" / "uniprot_pool.tsv", sep="\t", dtype=str)
    bench = pool[~pool.Entry.isin(used)].sample(args.pilot, random_state=0)
    cache = config.DATA / "pilot" / f"{args.stage}_{args.model}.jsonl"
items = [(e, prompts.state(s), qs) for e, s in zip(bench.Entry, bench.Sequence)]
out = models.run(args.model, cfg, items, cache)
print(f"{len(out)} responses; example: {json.dumps(next(iter(out.values())))[:400]}")
