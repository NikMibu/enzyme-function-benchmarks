import numpy as np
import pandas as pd
from scipy import stats


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def summarize(pred: pd.DataFrame, classes: list, gate: dict) -> dict:
    """pred: one row per protein with columns true, pred, and p_<class> probabilities."""
    n, k = len(pred), int((pred.true == pred.pred).sum())
    K = len(classes)
    chance = 1 / K
    lo, hi = wilson(k, n)
    p_binom = stats.binomtest(k, n, chance, alternative="greater").pvalue
    P = pred[[f"p_{c}" for c in classes]].to_numpy(float)
    P = P / P.sum(1, keepdims=True)
    y = pred.true.map({c: i for i, c in enumerate(classes)}).to_numpy()
    onehot = np.eye(K)[y]
    per_class, f1s = {}, []
    for c in classes:
        tp = int(((pred.true == c) & (pred.pred == c)).sum())
        n_true, n_pred = int((pred.true == c).sum()), int((pred.pred == c).sum())
        rec = tp / n_true if n_true else float("nan")
        prec = tp / n_pred if n_pred else 0.0
        f1s.append(0.0 if tp == 0 else 2 * prec * rec / (prec + rec))
        per_class[c] = {"n": n_true, "correct": tp, "accuracy": round(100 * rec, 1),
                        "precision": round(100 * prec, 1), "n_predicted": n_pred}
    return {
        "n": n, "n_classes": K, "correct": k,
        "accuracy": round(100 * k / n, 1),
        "accuracy_ci95": [round(100 * lo, 1), round(100 * hi, 1)],
        "chance_accuracy": round(100 * chance, 1),
        "p_vs_chance_one_sided": float(f"{p_binom:.3g}"),
        "macro_f1": round(100 * float(np.mean(f1s)), 1),
        "log_loss": round(float(-np.mean(np.log(np.clip(P[np.arange(n), y], 1e-12, 1)))), 3),
        "log_loss_uniform": round(float(np.log(K)), 3),
        "brier": round(float(np.mean(((P - onehot) ** 2).sum(1))), 3),
        "brier_uniform": round(float(1 - 1 / K), 3),
        "mean_top_probability": round(float(P.max(1).mean()), 3),
        "passes_gate": bool(lo > chance + gate["margin_over_chance"] and p_binom < gate["alpha"]),
        "per_class": per_class,
    }


def confusion(pred: pd.DataFrame, classes: list) -> pd.DataFrame:
    return (pd.crosstab(pd.Categorical(pred.true, classes), pd.Categorical(pred.pred, classes),
                        rownames=["true"], colnames=["predicted"], dropna=False)
            .reindex(index=classes, columns=classes, fill_value=0))
