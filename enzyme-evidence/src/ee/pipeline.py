"""Shared glue: load derived data, build candidates, score every integrator."""
import pandas as pd

from . import candidates, care, config
from .models import baselines, lgbm


def load_inputs(cfg):
    queries = pd.read_parquet(config.DERIVED / "queries.parquet")
    train = care.load_train(cfg)
    hits = pd.read_parquet(config.DERIVED / "hits.parquet")
    masked = candidates.mask_hits(hits, queries, train, cfg["retrieval"]["evalue"])
    return queries, train, masked


def score_all(cfg, queries, train, masked, k, eval_role, features=candidates.FEATURES):
    """Fit learned integrators on role=='fit' only, score role==eval_role. Returns (cands, scores)."""
    c = candidates.build(masked, queries, train, k, cfg["candidates"]["max_candidates"])
    c_fit, c_eval = c[c.role == "fit"], c[c.role == eval_role].copy()
    model = lgbm.fit(c_fit, cfg["lgbm"], cfg["seed"], features)
    scores = {
        "nearest_neighbour": baselines.nearest_neighbour(c_eval),
        "weighted_vote": baselines.weighted_vote(c_eval),
        "lightgbm": lgbm.score(model, c_eval, features),
        "oracle": baselines.oracle(c_eval),
    }
    return c_eval, scores, model
