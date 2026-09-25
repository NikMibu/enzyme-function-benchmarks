"""Query sets: CARE test splits plus a CARE-protocol replica built only from training data.

CARE ships no validation split, so we rebuild its test construction inside the training set:

* dev "30":    single-EC training entries alone in their 30% MMseqs2 cluster, whose EC occurs
               in at least one other entry. At most `max_per_ec3` per EC3 (as in CARE).
* dev "30-50": same, but alone at 50% and *not* alone at 30%.

The remaining eligible entries form the "fit" pool used to train learned integrators.

Every dev/fit query is itself in the retrieval database, so its hits are masked at the same
level CARE held test proteins out: self, identical sequences, and its whole 30% cluster (dev
"30") or 50% cluster (dev "30-50"). Fit queries additionally never see dev-eval entries.
"""
import numpy as np
import pandas as pd

from . import care

QUERY_COLS = ["Entry", "Sequence", "true_ec", "role", "split", "mask_level", "c30", "c50"]


def _eligible(train: pd.DataFrame) -> pd.DataFrame:
    t = train.copy()
    t["m30"] = t.c30.map(t.c30.value_counts())
    t["m50"] = t.c50.map(t.c50.value_counts())
    single = t.ecs.str.len() == 1
    unique_seq = ~t.Sequence.duplicated(keep=False)
    t["ec"] = t.ecs.str[0]
    all_ec_n = t.explode("ecs").ecs.value_counts()
    # EC must survive in the database after this entry is removed.
    t["ec_elsewhere"] = t.ec.map(all_ec_n).fillna(0) >= 2
    t["EC3"] = t.ec.str.rsplit(".", n=1).str[0]
    base = single & unique_seq & t.ec_elsewhere
    t["pool"] = np.select([base & (t.m30 == 1), base & (t.m50 == 1) & (t.m30 > 1)],
                          ["30", "30-50"], default="")
    return t


def build_queries(cfg) -> pd.DataFrame:
    rng_seed = cfg["seed"]
    train = care.load_train(cfg)
    t = _eligible(train)

    dev_parts = []
    for split in ["30", "30-50"]:
        pool = t[t.pool == split]
        dev = (pool.groupby("EC3", group_keys=False)
                   .apply(lambda g: g.sample(min(len(g), cfg["dev"]["max_per_ec3"]),
                                             random_state=rng_seed)))
        dev_parts.append(dev.assign(role="dev", split=split, mask_level=split))
    dev = pd.concat(dev_parts)

    rest = t[(t.pool != "") & ~t.Entry.isin(dev.Entry)]
    n_fit = min(len(rest), cfg["dev"]["fit_pool_max"])
    fit = rest.sample(n_fit, random_state=rng_seed).assign(role="fit")
    fit["split"] = fit.pool
    fit["mask_level"] = fit.pool

    trainq = pd.concat([dev, fit]).rename(columns={"ec": "true_ec"})

    tests = []
    for split in cfg["test_splits"]:
        te = care.load_test(cfg, split)
        tests.append(te.assign(role="test", split=split, mask_level="none", c30=None, c50=None))

    q = pd.concat([trainq[QUERY_COLS]] + [x[QUERY_COLS] for x in tests], ignore_index=True)
    assert q.Entry.is_unique, "query entries must be unique across roles"
    return q


def assert_no_test_leakage(queries: pd.DataFrame, train: pd.DataFrame) -> None:
    """Test proteins (by accession and exact sequence) must be absent from the retrieval DB."""
    test = queries[queries.role == "test"]
    assert not test.Entry.isin(train.Entry).any(), "test accession found in retrieval DB"
    assert not test.Sequence.isin(set(train.Sequence)).any(), "test sequence found in retrieval DB"
