"""Leakage guards: test proteins never in the retrieval DB; dev/fit hits masked like CARE."""
import pandas as pd
import pytest

from ee import candidates, care, config, splits


@pytest.fixture(scope="module")
def data(cfg):
    q = splits.build_queries(cfg)
    return q, care.load_train(cfg)


def test_test_split_absent_from_db(data):
    q, train = data
    splits.assert_no_test_leakage(q, train)


def test_roles_disjoint_and_dev_from_train_only(data):
    q, train = data
    assert q.Entry.is_unique
    assert q[q.role.isin(["dev", "fit"])].Entry.isin(train.Entry).all()


def test_dev_queries_are_isolated_like_care(data):
    q, train = data
    m30 = train.c30.map(train.c30.value_counts())
    m50 = train.c50.map(train.c50.value_counts())
    size30, size50 = dict(zip(train.Entry, m30)), dict(zip(train.Entry, m50))
    d30 = q[(q.role == "dev") & (q.split == "30")].Entry
    d50 = q[(q.role == "dev") & (q.split == "30-50")].Entry
    assert all(size30[e] == 1 for e in d30)
    assert all(size50[e] == 1 and size30[e] > 1 for e in d50)


def test_masking_on_real_hits(cfg, data):
    path = config.DERIVED / "hits.parquet"
    if not path.exists():
        pytest.skip("run scripts/02_retrieve.py first")
    q, train = data
    m = candidates.mask_hits(pd.read_parquet(path), q, train, cfg["retrieval"]["evalue"])
    qi, ti = q.set_index("Entry"), train.set_index("Entry")
    j = m.join(qi, on="query").join(ti[["Sequence", "c30", "c50"]], on="target", rsuffix="_t")
    assert not (j["query"] == j["target"]).any()
    assert not (j.Sequence == j.Sequence_t).any()
    assert not ((j.mask_level == "30") & (j.c30 == j.c30_t)).any()
    assert not ((j.mask_level == "30-50") & (j.c50 == j.c50_t)).any()
    dev = set(q[q.role == "dev"].Entry)
    assert not (j[j.role == "fit"]["target"].isin(dev)).any()
    assert j[j.role == "test"]["target"].isin(train.Entry).all()
