"""Score every model that has responses for a stage; write summaries, confusion matrices, figure.

    python scripts/03_evaluate.py --stage ec1
"""
import _common  # noqa: F401
import argparse
import json

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from ed import config, metrics, prompts  # noqa: E402

MODELS = ["jev", "laya-english", "laya-multilingual"]
BLUES = LinearSegmentedColormap.from_list("seq", ["#f7fafe", "#cde2fb", "#86b6ef", "#3987e5",
                                                  "#256abf", "#104281", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6a66"

ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["ec1", "ec4"], required=True)
args = ap.parse_args()
cfg = config.load()
bench = pd.read_csv(config.BENCH / f"{args.stage}.tsv", sep="\t", dtype=str)
truth_col = "ec1" if args.stage == "ec1" else "ec"
if args.stage == "ec1":
    classes = list(prompts.EC1_CLASSES)
    ticks = [f"{c} {prompts.EC1_CLASSES[c][0][:5]}." for c in classes]
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

# ---- figure: one row-normalised confusion matrix per model, shared colour scale
K = len(classes)
size = 3.2 if K <= 6 else 4.6
fig, axes = plt.subplots(1, len(done), figsize=(size * len(done) + 0.8, size + 0.9), squeeze=False)
for ax, (m, s, cm) in zip(axes[0], done):
    frac = cm.div(cm.sum(1), axis=0)
    ax.imshow(frac.to_numpy(), cmap=BLUES, vmin=0, vmax=1)
    for i in range(K):
        for j in range(K):
            v = int(cm.iat[i, j])
            if v:
                ax.text(j, i, v, ha="center", va="center", fontsize=8 if K <= 6 else 6,
                        color="white" if frac.iat[i, j] > 0.5 else INK)
    ax.set_xticks(range(K), ticks, rotation=90, fontsize=7, color=MUTED)
    ax.set_yticks(range(K), ticks, fontsize=7, color=MUTED)
    ax.set_xlabel("predicted", fontsize=8, color=MUTED)
    ax.set_title(f"{m}\n{s['accuracy']}% [{s['accuracy_ci95'][0]}, {s['accuracy_ci95'][1]}], "
                 f"chance {s['chance_accuracy']}%", fontsize=9, color=INK)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
axes[0][0].set_ylabel("true class", fontsize=8, color=MUTED)
fig.suptitle(f"{args.stage}: row-normalised confusion (cell text = count)", fontsize=10, color=INK)
fig.tight_layout()
(config.RESULTS / "figures").mkdir(parents=True, exist_ok=True)
fig.savefig(config.RESULTS / "figures" / f"confusion_{args.stage}.png", dpi=160, facecolor="white")
