"""Two figures: accuracy by measured identity, and selective accuracy vs coverage (<30% split)."""
import _common  # noqa: F401
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ee import config, evaluate  # noqa: E402

# Categorical slots 1-3 of the reference palette (validated all-pairs); oracle is neutral gray.
SERIES = {"nearest_neighbour": ("#2a78d6", "Nearest neighbour"),
          "weighted_vote": ("#eb6834", "Weighted vote"),
          "lightgbm": ("#1baf7a", "LightGBM ranker")}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({"font.size": 11, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})

out = config.RESULTS / "test"
fig_dir = config.RESULTS / "figures"
fig_dir.mkdir(parents=True, exist_ok=True)
preds = {m: pd.read_csv(out / "predictions" / f"{m}_per_query.csv") for m in [*SERIES, "oracle"]}

# --- Figure 1: level-4 accuracy by top-hit identity, all test splits pooled -----------------
order = list(evaluate.ID_LABELS)          # "no hit" proteins are always wrong; reported in the label
fig, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(len(order))
ns = None
for m, (col, label) in SERIES.items():
    p = preds[m].assign(id_bin=preds[m].id_bin.fillna("no hit"))
    g = p.groupby("id_bin").L4.agg(["mean", "size"]).reindex(order)
    ns = g["size"]
    ax.plot(x, 100 * g["mean"], color=col, lw=2, marker="o", ms=8, mec=SURFACE, mew=2, label=label)
po = preds["oracle"].assign(id_bin=preds["oracle"].id_bin.fillna("no hit"))
go = po.groupby("id_bin").oracle.mean().reindex(order)
ax.plot(x, 100 * go, color="#8a8983", lw=2, ls="--", label="Oracle (true EC retrieved)")
ax.set_xticks(x, [f"{o}\nn={int(n)}" for o, n in zip(order, ns)])
n_nohit = int(preds["nearest_neighbour"].id_bin.isna().sum())
ax.set_xlabel(f"Top-hit identity to CARE train (test splits pooled; {n_nohit} no-hit proteins omitted)")
ax.set_ylabel("EC level-4 accuracy (%)")
ax.set_ylim(30, 100)
ax.yaxis.grid(True, color=GRID, lw=1)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left", fontsize=9, labelcolor=INK)
ax.set_title("Top-1 accuracy tracks retrieval identity; integrators barely differ", loc="left", color=INK)
fig.tight_layout()
fig.savefig(fig_dir / "accuracy_by_identity.png", dpi=200)

# --- Figure 2: selective accuracy vs coverage on the <30% split -----------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.2))
curves = {"nearest_neighbour": "NN, confidence = identity x coverage",
          "lightgbm": "LightGBM, confidence = P(correct)"}
for m, label in curves.items():
    p = preds[m][preds[m].split.astype(str) == "30"]
    cov, acc = evaluate.risk_coverage(p.conf.to_numpy(float), p.correct.to_numpy(float))
    keep = cov >= 0.05                     # the first ~20 proteins are noise
    ax.plot(100 * cov[keep], 100 * acc[keep], color=SERIES[m][0], lw=2, label=label)
ax.set_xlabel("Coverage: fraction of <30%-split proteins answered (most confident first, %)")
ax.set_ylabel("Level-4 accuracy on answered (%)")
ax.set_ylim(50, 100)
ax.set_xlim(0, 100)
ax.yaxis.grid(True, color=GRID, lw=1)
ax.set_axisbelow(True)
ax.legend(frameon=False, loc="lower left", fontsize=9, labelcolor=INK)
ax.set_title("Selective prediction below 30% identity (CARE test)", loc="left", color=INK)
fig.tight_layout()
fig.savefig(fig_dir / "selective_accuracy_30.png", dpi=200)
print("figures written to", fig_dir)
