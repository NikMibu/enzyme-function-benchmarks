"""Candidate features and rendering on a tiny synthetic example."""
import pandas as pd

from ee import candidates, render

TRAIN = pd.DataFrame({
    "Entry": ["T1", "T2", "T3", "Q_DEV"],
    "Sequence": ["AAA", "CCC", "DDD", "EEE"],
    "ecs": [["1.1.1.1"], ["1.1.1.1", "2.2.2.2"], ["1.1.1.2"], ["1.1.1.2"]],
    "c30": ["a", "b", "c", "d"], "c50": ["a", "b", "c", "d"],
})
QUERIES = pd.DataFrame({
    "Entry": ["QT", "Q_DEV"], "Sequence": ["GGG", "EEE"], "true_ec": ["1.1.1.2", "1.1.1.2"],
    "role": ["test", "dev"], "split": ["30", "30"], "mask_level": ["none", "30"],
    "c30": [None, "d"], "c50": [None, "d"],
})
HITS = pd.DataFrame({
    "query": ["QT", "QT", "QT", "Q_DEV", "Q_DEV"],
    "target": ["T1", "T2", "T3", "Q_DEV", "T3"],
    "fident": [0.4, 0.3, 0.25, 1.0, 0.5], "alnlen": [100] * 5, "qlen": [100] * 5,
    "tlen": [110, 90, 100, 100, 100], "qcov": [0.9, 0.8, 0.7, 1.0, 0.9], "tcov": [0.9] * 5,
    "evalue": [1e-20, 1e-10, 1e-5, 0, 1e-30], "bits": [200.0, 100.0, 50.0, 500.0, 300.0],
})


def test_masking_and_features():
    m = candidates.mask_hits(HITS, QUERIES, TRAIN, 1e-3)
    assert not ((m["query"] == "Q_DEV") & (m["target"] == "Q_DEV")).any()
    c = candidates.build(m, QUERIES, TRAIN, k=10, max_candidates=10)
    qt = c[c["query"] == "QT"].set_index("ec")
    assert set(qt.index) == {"1.1.1.1", "2.2.2.2", "1.1.1.2"}
    assert qt.loc["1.1.1.1", "n_support"] == 2 and qt.loc["1.1.1.1", "best_rank"] == 1
    # T2 has two ECs, so it splits its 100 bits: (200 + 50) / 350
    assert abs(qt.loc["1.1.1.1", "bits_share"] - 250 / 350) < 1e-9
    assert qt.loc["1.1.1.2", "label"] and not qt.loc["1.1.1.1", "label"]
    # the dev query's own EC count is removed from the frequency prior (2 entries -> 1)
    dv = c[c["query"] == "Q_DEV"].set_index("ec")
    assert abs(dv.loc["1.1.1.2", "log_train_freq"] - __import__("math").log1p(1)) < 1e-9


def test_render_hides_identity_and_shuffles():
    m = candidates.mask_hits(HITS, QUERIES, TRAIN, 1e-3)
    c = candidates.build(m, QUERIES, TRAIN, k=10, max_candidates=10)
    g = c[c["query"] == "QT"]
    r = render.render("QT", g, {"1.1.1.1": "alcohol dehydrogenase"})
    text = r["state"] + str(r["questions"])
    for leak in ["QT", "T1", "T2", "T3", "GGG"]:
        assert leak not in text.replace("QUERY", "")
    crit = r["questions"]["ec"]["criteria"]
    assert render.NONE_LABEL in crit and len(crit) == len(g) + 1
    assert set(r["label_to_ec"].values()) == set(g.ec)
    assert render.render("QT", g, {})["label_to_ec"] == render.render("QT", g, {})["label_to_ec"]
