"""CARE-native accuracy plus paired comparisons and calibration / selective-prediction metrics."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from .care import ec_level

NO_PRED = "0.0.0.0"   # what CARE's notebook substitutes for a missing prediction
ID_BINS = [-0.01, 0.25, 0.30, 0.40, 0.50, 1.0]
ID_LABELS = ["<25%", "25-30%", "30-40%", "40-50%", ">50%"]


def care_accuracy(pred_ecs: list, true_ecs: str, level: int) -> float:
    """CARE's per-protein accuracy: mean over true ECs of [best prediction reaches `level`]."""
    trues = true_ecs.split(";")
    return float(np.mean([max(ec_level(p, t) for p in pred_ecs) >= level for t in trues]))


def predictions(queries: pd.DataFrame, c: pd.DataFrame, score: pd.Series) -> pd.DataFrame:
    """One row per query (queries without any hit get no prediction and confidence 0)."""
    c = c.assign(score=score.values).sort_values(["query", "score"], ascending=[True, False])
    ranked = c.groupby("query").ec.apply(list)
    top = c.groupby("query").first()
    p = queries[["Entry", "split", "true_ec"]].rename(columns={"Entry": "query"}).set_index("query")
    p["ranked"] = ranked.reindex(p.index)
    p["ranked"] = p.ranked.apply(lambda r: r if isinstance(r, list) else [])
    p["pred"] = p.ranked.apply(lambda r: r[0] if r else NO_PRED)
    p["conf"] = top.score.reindex(p.index).fillna(0.0)
    p["n_cands"] = top.n_cands.reindex(p.index).fillna(0).astype(int)
    p["top_fid"] = top.top_fid.reindex(p.index)
    p["top_qcov"] = top.top_qcov.reindex(p.index)
    p["oracle"] = c.groupby("query").label.any().reindex(p.index).fillna(False)
    for L in (1, 2, 3, 4):
        p[f"L{L}"] = [care_accuracy([pr], t, L) for pr, t in zip(p.pred, p.true_ec)]
    p["correct"] = [pr in t.split(";") for pr, t in zip(p.pred, p.true_ec)]
    p["id_bin"] = pd.cut(p.top_fid.fillna(0), ID_BINS, labels=ID_LABELS)
    p.loc[p.top_fid.isna(), "id_bin"] = np.nan
    return p.reset_index()


def ece(conf: np.ndarray, correct: np.ndarray, n_bins: int = 10) -> float:
    bins = np.minimum((conf * n_bins).astype(int), n_bins - 1)
    return float(sum(abs(correct[bins == b].mean() - conf[bins == b].mean()) * (bins == b).mean()
                     for b in range(n_bins) if (bins == b).any()))


def risk_coverage(conf: np.ndarray, correct: np.ndarray):
    """Coverage and selective accuracy when keeping the most confident predictions first."""
    order = np.argsort(-conf, kind="stable")
    acc = np.cumsum(correct[order]) / np.arange(1, len(order) + 1)
    cov = np.arange(1, len(order) + 1) / len(order)
    return cov, acc


def summarize(p: pd.DataFrame) -> dict:
    conf, correct = p.conf.to_numpy(float), p.correct.to_numpy(float)
    cov, acc = risk_coverage(conf, correct)
    contested = p[p.n_cands >= 2]
    out = {"n": len(p)}
    out.update({f"L{L}": round(100 * p[f"L{L}"].mean(), 1) for L in (1, 2, 3, 4)})
    out.update({
        "oracle_L4": round(100 * p.oracle.mean(), 1),
        "no_hit": round(100 * (p.n_cands == 0).mean(), 1),
        "n_contested": len(contested),
        "contested_L4": round(100 * contested.L4.mean(), 1) if len(contested) else None,
        "ECE": round(ece(conf, correct), 3),
        "Brier": round(float(np.mean((conf - correct) ** 2)), 3),
        "AURC_risk": round(float(np.mean(1 - acc)), 3),       # area under risk-coverage
        "sel_acc@50%": round(100 * float(acc[max(0, int(0.5 * len(acc)) - 1)]), 1),
        "sel_acc@80%": round(100 * float(acc[max(0, int(0.8 * len(acc)) - 1)]), 1),
    })
    return out


def paired(p_model: pd.DataFrame, p_ref: pd.DataFrame, n_boot: int = 2000, seed: int = 0) -> dict:
    """Model minus reference on L4 accuracy: bootstrap 95% CI, exact McNemar, wins / losses."""
    a = p_model.set_index("query").L4
    b = p_ref.set_index("query").L4.reindex(a.index)
    d = (a - b).to_numpy()
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(0, len(d), size=(n_boot, len(d)))].mean(1)
    wins, losses = int((d > 0).sum()), int((d < 0).sum())
    p_mc = binomtest(wins, wins + losses).pvalue if wins + losses else 1.0
    return {"delta_L4": round(100 * d.mean(), 1),
            "ci95": [round(100 * np.percentile(boots, 2.5), 1), round(100 * np.percentile(boots, 97.5), 1)],
            "wins": wins, "losses": losses, "mcnemar_p": round(float(p_mc), 4)}


def by_identity(p: pd.DataFrame) -> pd.DataFrame:
    g = p.assign(id_bin=p.id_bin.astype(object).fillna("no hit")).groupby("id_bin")
    return g.agg(n=("L4", "size"), L4=("L4", "mean"), oracle=("oracle", "mean"))


def to_care_csv(p: pd.DataFrame, path, k: int = 10) -> None:
    """CARE results_summary format: Entry, EC number, 0..k-1 (ranked ECs)."""
    cols = {"Entry": p["query"], "EC number": p.true_ec}
    for i in range(k):
        cols[str(i)] = p.ranked.apply(lambda r: r[i] if len(r) > i else None)
    pd.DataFrame(cols).to_csv(path, index=False)


def aurc(conf: np.ndarray, correct: np.ndarray) -> float:
    return float(np.mean(1 - risk_coverage(conf, correct)[1]))


def paired_aurc(conf_a, correct_a, conf_b, correct_b, n_boot: int = 2000, seed: int = 0) -> dict:
    """AURC(a) - AURC(b) over the same proteins, with a paired bootstrap 95% CI (negative = a better)."""
    rng = np.random.default_rng(seed)
    n = len(conf_a)
    d = [aurc(conf_a[i], correct_a[i]) - aurc(conf_b[i], correct_b[i])
         for i in (rng.integers(0, n, n) for _ in range(n_boot))]
    return {"delta_AURC": round(aurc(conf_a, correct_a) - aurc(conf_b, correct_b), 3),
            "ci95": [round(float(np.percentile(d, 2.5)), 3), round(float(np.percentile(d, 97.5)), 3)]}
