"""Our scorer must reproduce CARE's own published Task 1 numbers from their shipped CSVs."""
import pandas as pd
import pytest

from ee import config
from ee.evaluate import care_accuracy

# performance_evaluation.ipynb, k=1, level 4/3/2/1
EXPECTED = {
    ("BLAST", "30"): (51.4, 60.0, 62.5, 65.7),
    ("BLAST", "30-50"): (81.1, 87.9, 90.7, 92.3),
    ("BLAST", "price"): (35.1, 70.9, 78.4, 78.4),
    ("CLEAN", "30"): (55.1, 68.8, 74.8, 84.5),
}


@pytest.mark.parametrize("method,split", list(EXPECTED))
def test_reproduces_care(cfg, method, split):
    df = pd.read_csv(config.care_dir(cfg) / "task1_baselines" / "results_summary" / method
                     / f"{split}_protein_test_results_df.csv")
    col0 = df["0"].fillna("0.0.0.0").astype(str)
    # CARE expands '; '-separated BLAST/Foldseek outputs (a multi-EC neighbour counts if any EC matches)
    preds = col0.str.split("; ") if method == "BLAST" else col0.apply(lambda x: [x])
    got = tuple(round(100 * sum(care_accuracy(p, t, L) for p, t in zip(preds, df["EC number"])) / len(df), 1)
                for L in (4, 3, 2, 1))
    assert got == EXPECTED[(method, split)]


def test_multi_ec_truth_is_averaged():
    assert care_accuracy(["4.2.1.68"], "4.2.1.68;4.2.1.90", 4) == 0.5
    assert care_accuracy(["4.2.1.68"], "4.2.1.68;4.2.1.90", 3) == 1.0
