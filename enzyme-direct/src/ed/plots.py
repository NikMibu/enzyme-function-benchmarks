import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

BLUES = LinearSegmentedColormap.from_list("seq", ["#f7fafe", "#cde2fb", "#86b6ef", "#3987e5",
                                                  "#256abf", "#104281", "#0d366b"])
INK, MUTED, GRID, ACCENT = "#1f1f1e", "#6b6a66", "#e4e2dc", "#256abf"


def confusions(panels, ticks, title, path):
    """panels: [(name, summary, confusion DataFrame)]. Rows are divided by the class size, so a
    model that abstains shows lighter rows rather than inflated ones."""
    K = len(ticks)
    size = 3.2 if K <= 6 else 4.6
    fig, axes = plt.subplots(1, len(panels), figsize=(size * len(panels) + 0.8, size + 0.9), squeeze=False)
    for ax, (m, s, cm) in zip(axes[0], panels):
        n_true = [s["per_class"][c]["n"] for c in cm.index]
        frac = cm.div(n_true, axis=0)
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
    fig.suptitle(title, fontsize=10, color=INK)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, facecolor="white")
    plt.close(fig)


def ladder(rows, chance, path):
    """rows: [(label, group, accuracy, lo, hi)] top to bottom; one dot + 95% CI per run."""
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(rows) + 1.2))
    y = list(range(len(rows)))[::-1]
    for yi, (lab, grp, acc, lo, hi) in zip(y, rows):
        col = ACCENT if grp == "model" else MUTED
        ax.plot([lo, hi], [yi, yi], color=col, lw=2, solid_capstyle="round")
        ax.plot(acc, yi, "o", ms=8, color=col, mec="white", mew=2)
        ax.text(hi + 1.2, yi, f"{acc:.1f}%", va="center", fontsize=8, color=INK)
    ax.axvline(100 * chance, color=INK, lw=1, ls=(0, (3, 3)))
    ax.text(100 * chance, len(rows) - 0.35, " chance", fontsize=8, color=MUTED, va="bottom")
    ax.set_yticks(y, [r[0] for r in rows], fontsize=8, color=INK)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, len(rows) - 0.1)
    ax.set_xlabel("EC level-1 accuracy, % (95% Wilson CI), 210 held-out proteins", fontsize=8, color=MUTED)
    ax.tick_params(axis="x", labelsize=8, colors=MUTED, length=0)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor="white")
    plt.close(fig)
