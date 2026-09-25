"""Context ladder on ec1: Jev at each context level against its matched control, plus the NN ceiling.

    python scripts/05_compare.py      # after 03_evaluate.py --stage ec1 and 04_controls.py
"""
import _common  # noqa: F401
import json

import numpy as np
import pandas as pd
from scipy import stats

from ed import config, metrics, plots, prompts

cfg = config.load()
out = config.RESULTS / "ec1"
classes = list(prompts.EC1_CLASSES)
LABELS = {
    "jev": "Jev, sequence only",
    "laya-english": "Laya English, sequence only",
    "laya-multilingual": "Laya multilingual, sequence only",
    "jev-ctx1": "Jev + properties & composition",
    "logreg-ctx1": "LogReg on properties & composition",
    "jev-ctx2": "Jev + properties, composition & motifs",
    "logreg-ctx2": "LogReg on properties, composition & motifs",
    "logreg-ctx1-care": "LogReg on properties & composition (3,000 CARE proteins)",
    "esm-knn": "ESM-2 nearest neighbour (3,000 CARE proteins)",
    "esm-probe": "ESM-2 650M probe (3,000 CARE proteins, no homologs)",
    "jev-homologs": "Jev + homolog evidence (enzyme-evidence setup)",
    "nn-care-train": "Nearest neighbour, Swiss-Prot (CARE train)",
}
MODELS = {"jev", "laya-english", "laya-multilingual", "jev-ctx1", "jev-ctx2", "jev-homologs"}
PAIRS = [("jev-ctx1", "logreg-ctx1"), ("jev-ctx2", "logreg-ctx2"),
         ("jev-ctx1", "jev"), ("jev-ctx2", "jev"), ("jev-ctx2", "jev-ctx1")]

runs = {r: (json.load(open(out / f"{r}_summary.json")), pd.read_csv(out / f"{r}_predictions.csv", dtype={"true": str, "pred": str}))
        for r in LABELS if (out / f"{r}_summary.json").exists()}


def paired(a: str, b: str, n_boot: int = 10000) -> dict:
    """Exact McNemar test and bootstrap 95% CI for accuracy(a) - accuracy(b), same proteins."""
    pa, pb = runs[a][1].set_index("entry"), runs[b][1].set_index("entry")
    ca = (pa.true == pa.pred).to_numpy()
    cb = (pb.loc[pa.index].true == pb.loc[pa.index].pred).to_numpy()
    wins, losses = int((ca & ~cb).sum()), int((~ca & cb).sum())
    p = stats.binomtest(wins, wins + losses, 0.5).pvalue if wins + losses else 1.0
    rng = np.random.default_rng(cfg["seed"])
    idx = rng.integers(0, len(ca), (n_boot, len(ca)))
    d = 100 * (ca[idx].mean(1) - cb[idx].mean(1))
    return {"a": a, "b": b, "delta_accuracy": round(100 * (ca.mean() - cb.mean()), 1),
            "ci95": [round(float(np.percentile(d, 2.5)), 1), round(float(np.percentile(d, 97.5)), 1)],
            "a_only_correct": wins, "b_only_correct": losses, "mcnemar_p": float(f"{p:.3g}")}


rows = []
for r, (s, _) in runs.items():
    rows.append({"run": r, "label": LABELS[r], "accuracy": s["accuracy"],
                 "ci95_low": s["accuracy_ci95"][0], "ci95_high": s["accuracy_ci95"][1],
                 "p_vs_chance": s["p_vs_chance_one_sided"], "macro_f1": s["macro_f1"],
                 "log_loss": s["log_loss"], "passes_gate": s["passes_gate"],
                 **{f"acc_{c}": s["per_class"][c]["accuracy"] for c in classes}})
table = pd.DataFrame(rows)
table.to_csv(out / "ladder.csv", index=False)
tests = [paired(a, b) for a, b in PAIRS if a in runs and b in runs]
reading = {}
for lvl in cfg["context"]["levels"]:
    j, c = f"jev-ctx{lvl}", f"logreg-ctx{lvl}"
    if j in runs and c in runs:
        t = next(t for t in tests if t["a"] == j and t["b"] == c)
        reading[f"level{lvl}"] = {
            "jev_uses_context": runs[j][0]["passes_gate"],
            "jev_adds_to_context": bool(t["mcnemar_p"] < 0.05 and t["delta_accuracy"] > 0),
            "control_passes_gate": runs[c][0]["passes_gate"],
        }
json.dump({"paired_tests": tests, "pre_registered_reading": reading}, open(out / "ladder_tests.json", "w"), indent=2)
print(table[["run", "accuracy", "ci95_low", "ci95_high", "p_vs_chance", "macro_f1", "log_loss", "passes_gate"]].to_string(index=False))
print(json.dumps(tests, indent=1))
print(json.dumps(reading, indent=1))

fig_dir = config.RESULTS / "figures"
plots.ladder([(LABELS[r], "model" if r in MODELS else "control", runs[r][0]["accuracy"],
               *runs[r][0]["accuracy_ci95"]) for r in LABELS if r in runs],
             1 / len(classes), fig_dir / "ladder_ec1.png")
ctx = [r for r in ["jev-ctx1", "jev-ctx2", "logreg-ctx1", "logreg-ctx2", "nn-care-train"] if r in runs]
ticks = [f"{c} {prompts.EC1_CLASSES[c][0]}" for c in classes]
plots.confusions([(r, runs[r][0], pd.read_csv(out / f"{r}_confusion.csv", index_col=0, dtype={"true": str})
                   .rename(index=str, columns=str)) for r in ctx],
                 ticks, "ec1 with context, and controls: confusion (cell text = count)",
                 fig_dir / "confusion_ec1_context.png")
