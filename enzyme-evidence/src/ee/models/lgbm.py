"""LightGBM pointwise ranker: P(candidate EC is correct | its evidence features)."""
import lightgbm as lgb
import pandas as pd

from ..candidates import FEATURES


def fit(c_fit: pd.DataFrame, params: dict, seed: int, features=FEATURES) -> lgb.LGBMClassifier:
    model = lgb.LGBMClassifier(objective="binary", random_state=seed, verbose=-1, **params)
    model.fit(c_fit[features], c_fit.label.astype(int))
    return model


def score(model: lgb.LGBMClassifier, c: pd.DataFrame, features=FEATURES) -> pd.Series:
    return pd.Series(model.predict_proba(c[features])[:, 1], index=c.index)
