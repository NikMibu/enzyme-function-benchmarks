"""Build dev / fit / test query sets and check test leakage against the retrieval DB."""
import _common  # noqa: F401
from ee import care, config, splits

cfg = config.load()
queries = splits.build_queries(cfg)
train = care.load_train(cfg)
splits.assert_no_test_leakage(queries, train)
config.DERIVED.mkdir(parents=True, exist_ok=True)
queries.to_parquet(config.DERIVED / "queries.parquet", index=False)
print(queries.groupby(["role", "split"]).size().to_string())
