"""Score every model that has responses for a stage; write summaries, confusion matrices, figure.

    python scripts/03_evaluate.py --stage ec1
"""
import _common  # noqa: F401
import argparse
import json

import pandas as pd

from ed import config, metrics, plots, prompts

MODELS = ["jev", "laya-english", "laya-multilingual", "jev-ctx1", "jev-ctx2"]

ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["ec1", "ec4"], required=True)
args = ap.parse_args()
cfg = config.load()
bench = pd.read_csv(config.BENCH / f"{args.stage}.tsv", sep="\t", dtype=str)
truth_col = "ec1" if args.stage == "ec1" else "ec"
if args.stage == "ec1":
    classes = list(prompts.EC1_CLASSES)
    ticks = [f"{c} {prompts.EC1_CLASSES[c][0]}" for c in classes]
else:
    classes = list(pd.read_csv(config.BENCH / "ec4_labels.tsv", sep="\t", dtype=str).ec)
    ticks = classes
out = config.RESULTS / args.stage

done, rows = [], []
for m in MODELS:
    cache = out / f"{m}_responses.jsonl"
    if not cache.exists():
        continue
    resp = {r["entry"]: r["response"] for r in map(json.loads, cache.read_text().splitlines())}
    missing = set(bench.Entry) - set(resp)
    if missing:
        print(f"{m}: {len(missing)} proteins without a response; skipped")
        continue
    recs = []
    for e, t in zip(bench.Entry, bench[truth_col]):
        ans = resp[e]["answers"][args.stage]
        probs = {prompts.label_to_class(args.stage, k): float(v) for k, v in ans["probabilities"].items()}
        recs.append({"entry": e, "true": t, "pred": max(probs, key=probs.get),
                     "input_tokens": resp[e]["usage"]["input_tokens"],
                     **{f"p_{c}": probs.get(c, 0.0) for c in classes}})
    pred = pd.DataFrame(recs)
    pred.to_csv(out / f"{m}_predictions.csv", index=False)
    s = metrics.summarize(pred, classes, cfg["gate"])
    s["model"] = resp[bench.Entry[0]]["model"] if m == "jev" else f"convaiinnovations/laya ({m})"
    s["max_input_tokens"] = int(pred.input_tokens.max())
    json.dump(s, open(out / f"{m}_summary.json", "w"), indent=2)
    cm = metrics.confusion(pred, classes)
    cm.to_csv(out / f"{m}_confusion.csv")
    done.append((m, s, cm))
    rows.append({"model": m, "n": s["n"], "accuracy": s["accuracy"],
                 "ci95_low": s["accuracy_ci95"][0], "ci95_high": s["accuracy_ci95"][1],
                 "chance": s["chance_accuracy"], "p_vs_chance": s["p_vs_chance_one_sided"],
                 "macro_f1": s["macro_f1"], "log_loss": s["log_loss"],
                 "log_loss_uniform": s["log_loss_uniform"], "passes_gate": s["passes_gate"],
                 **{f"acc_{c}": s["per_class"][c]["accuracy"] for c in classes}})
if not done:
    raise SystemExit("no complete model runs")
table = pd.DataFrame(rows)
table.to_csv(out / "summary.csv", index=False)
print(table.drop(columns=[c for c in table if c.startswith("acc_")]).to_string(index=False))

# ---- figure: zero-shot runs only (context runs are plotted by 05_compare.py)
zero_shot = [d for d in done if "-ctx" not in d[0]]
if zero_shot:
    plots.confusions(zero_shot, ticks, f"{args.stage}: sequence only, confusion (cell text = count)",
                     config.RESULTS / "figures" / f"confusion_{args.stage}.png")
