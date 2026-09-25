"""The benchmarks are held out from CARE's test sets, balanced, non-redundant, and prompts carry
nothing but the sequence."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ed import config, data, prompts  # noqa: E402

CFG = config.load()
EC1 = pd.read_csv(config.BENCH / "ec1.tsv", sep="\t", dtype=str)
EC4 = pd.read_csv(config.BENCH / "ec4.tsv", sep="\t", dtype=str)


@pytest.fixture(scope="module")
def care_test():
    return data.care_exclusion(CFG)


@pytest.mark.parametrize("bench", [EC1, EC4], ids=["ec1", "ec4"])
def test_absent_from_care_test_sets(bench, care_test):
    acc, seqs = care_test
    assert not set(bench.Entry) & acc
    assert not set(bench.Sequence) & seqs


def test_ec1_balanced_and_one_per_family():
    counts = EC1.ec1.value_counts()
    assert sorted(counts.index) == CFG["stages"]["ec1"]["classes"]
    assert counts.nunique() == 1 and counts.iloc[0] == CFG["stages"]["ec1"]["per_class"]
    assert EC1.cluster.is_unique


def test_ec4_balanced_disjoint_and_one_per_family():
    counts = EC4.ec.value_counts()
    assert 10 <= len(counts) <= 20 and counts.nunique() == 1
    assert EC4.cluster.is_unique
    assert not set(EC4.Entry) & set(EC1.Entry)
    assert not set(EC4.Sequence) & set(EC1.Sequence)


def test_single_complete_ec_and_length_window():
    u = CFG["uniprot"]
    for b in (EC1, EC4):
        assert b.ec.str.fullmatch(r"\d+\.\d+\.\d+\.\d+").all()
        assert b.Length.astype(int).between(u["min_len"], u["max_len"]).all()
        assert (b.Sequence.str.len() == b.Length.astype(int)).all()


def test_state_is_sequence_only():
    for _, r in pd.concat([EC1, EC4]).iterrows():
        st = prompts.state(r.Sequence)
        assert st == f"Protein amino-acid sequence ({len(r.Sequence)} residues, one-letter code):\n{r.Sequence}"
        assert r.Entry not in st and r.ec not in st


def test_ec1_options_do_not_depend_on_the_protein():
    q = prompts.ec1_question()["ec1"]
    assert len(q["criteria"]) == 6
    assert set(prompts.LABEL_TO_EC1.values()) == set(prompts.EC1_CLASSES)
