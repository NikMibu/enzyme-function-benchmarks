"""MMseqs2 search of all queries against the CARE training set (~20 min on 4 CPU cores)."""
import _common  # noqa: F401
import pandas as pd

from ee import care, config, retrieval

cfg = config.load()
queries = pd.read_parquet(config.DERIVED / "queries.parquet")
train = care.load_train(cfg)
hits = retrieval.search(cfg, queries, train, config.DERIVED / "mmseqs")
hits.to_parquet(config.DERIVED / "hits.parquet", index=False)
print(f"{len(hits)} hits for {hits['query'].nunique()} / {len(queries)} queries")
