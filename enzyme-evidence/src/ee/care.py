"""Loading CARE Task 1 data at the protein (entry) level."""
from functools import lru_cache
from pathlib import Path

import pandas as pd

from . import config


def _split_dir(cfg) -> Path:
    return config.care_dir(cfg) / "splits" / "task1"


@lru_cache(maxsize=1)
def _train_cached(split_dir: str) -> pd.DataFrame:
    rows = pd.read_csv(Path(split_dir) / "protein_train.csv")
    # protein_train.csv has one row per (entry, EC); collapse to one row per entry.
    df = (rows.groupby("Entry", sort=True)
              .agg(Sequence=("Sequence", "first"),
                   ecs=("EC number", lambda x: sorted(set(x))),
                   c30=("clusterRes30", "first"),
                   c50=("clusterRes50", "first"))
              .reset_index())
    df["Length"] = df.Sequence.str.len()
    return df


def load_train(cfg) -> pd.DataFrame:
    """One row per training entry: Entry, Sequence, ecs (sorted list), c30, c50, Length."""
    return _train_cached(str(_split_dir(cfg))).copy()


def load_test(cfg, split: str) -> pd.DataFrame:
    """Entry, Sequence, true_ec for a CARE Task 1 test split ('30', '30-50', 'price').

    true_ec is ';'-joined when a protein has several ECs (3 proteins in Price); CARE scores
    those as the mean over true ECs, which evaluate.care_accuracy reproduces.
    """
    df = pd.read_csv(_split_dir(cfg) / f"{split}_protein_test.csv")
    return pd.DataFrame({"Entry": df.Entry, "Sequence": df.Sequence, "true_ec": df["EC number"]})


def load_ec_names(cfg) -> dict:
    """EC -> short enzyme name (last segment of CARE's text2EC hierarchy description)."""
    df = pd.read_csv(config.care_dir(cfg) / "processed_data" / "text2EC.csv")
    return {ec: str(t).split("; ")[-1] for ec, t in zip(df["EC number"], df["Text"])}


def ec_level(pred: str, true: str) -> int:
    """Number of leading EC digits that match (CARE's accuracy level, 0-4)."""
    n = 0
    for p, t in zip(str(pred).split("."), str(true).split(".")):
        if p != t:
            break
        n += 1
    return n
