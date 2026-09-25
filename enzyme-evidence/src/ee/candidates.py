"""Retrieved hits -> masked top-k neighbours -> one row per (query, candidate EC) with features.

Nothing about the query except its length and its hits enters the features. Labels are only
attached for evaluation / supervised fitting.
"""
import numpy as np
import pandas as pd

FEATURES = [
    "best_rank", "best_bits", "best_fid", "best_qcov", "best_tcov", "len_ratio",
    "n_support", "frac_support", "bits_share", "bits_ratio", "fid_gap", "ec3_share",
    "log_train_freq", "is_prelim", "n_cands", "n_hits", "top_fid", "top_qcov", "qlen",
]


def mask_hits(hits: pd.DataFrame, queries: pd.DataFrame, train: pd.DataFrame,
              evalue: float) -> pd.DataFrame:
    """Drop hits a query must not see (see splits.py) and hits above the e-value cut-off."""
    q = queries[["Entry", "Sequence", "role", "mask_level", "c30", "c50"]].rename(
        columns={"Entry": "query", "Sequence": "q_seq", "c30": "q_c30", "c50": "q_c50"})
    t = train[["Entry", "Sequence", "c30", "c50"]].rename(
        columns={"Entry": "target", "Sequence": "t_seq", "c30": "t_c30", "c50": "t_c50"})
    h = hits[hits.evalue <= evalue].merge(q, on="query").merge(t, on="target")
    dev_entries = set(queries.loc[queries.role == "dev", "Entry"])
    drop = (
        (h["query"] == h["target"])
        | (h.q_seq == h.t_seq)
        | ((h.mask_level == "30") & (h.q_c30 == h.t_c30))
        | ((h.mask_level == "30-50") & (h.q_c50 == h.t_c50))
        | ((h.role == "fit") & h["target"].isin(dev_entries))
    )
    return h.loc[~drop, hits.columns.tolist()]


def build(masked_hits: pd.DataFrame, queries: pd.DataFrame, train: pd.DataFrame,
          k: int, max_candidates: int) -> pd.DataFrame:
    ecs = train.set_index("Entry").ecs
    ec_freq = train.explode("ecs").ecs.value_counts()

    h = masked_hits.sort_values(["query", "bits"], ascending=[True, False])
    h = h.groupby("query").head(k).copy()
    h["rank"] = h.groupby("query").cumcount() + 1
    h["t_ecs"] = h["target"].map(ecs)
    h["w"] = h.bits / h.t_ecs.str.len()           # a multi-EC neighbour splits its vote
    per_q = h.groupby("query").agg(n_hits=("rank", "size"), tot_bits=("bits", "sum"),
                                   top_bits=("bits", "first"), top_fid=("fident", "first"),
                                   top_qcov=("qcov", "first"))

    x = h.explode("t_ecs").rename(columns={"t_ecs": "ec"})
    x["ec3"] = x.ec.str.rsplit(".", n=1).str[0]
    first = x.sort_values("rank").groupby(["query", "ec"]).first()
    c = x.groupby(["query", "ec"]).agg(n_support=("rank", "size"), w_sum=("w", "sum"))
    c["best_rank"] = first["rank"]
    c["best_bits"] = first.bits
    c["best_fid"] = first.fident
    c["best_qcov"] = first.qcov
    c["best_tcov"] = first.tcov
    c["len_ratio"] = first.tlen / first.qlen
    c["qlen"] = first.qlen
    c = c.reset_index()

    # EC3-level agreement: bits share of neighbours carrying any EC in the candidate's EC3.
    ec3_w = (x.drop_duplicates(["query", "target", "ec3"])
              .groupby(["query", "ec3"]).bits.sum().rename("ec3_bits").reset_index())
    c["ec3"] = c.ec.str.rsplit(".", n=1).str[0]
    c = c.merge(ec3_w, on=["query", "ec3"]).join(per_q, on="query")

    c["frac_support"] = c.n_support / c.n_hits
    c["bits_share"] = c.w_sum / c.tot_bits
    c["ec3_share"] = c.ec3_bits / c.tot_bits
    c["bits_ratio"] = c.best_bits / c.top_bits
    c["fid_gap"] = c.top_fid - c.best_fid
    c["is_prelim"] = c.ec.str.contains("n").astype(int)

    # EC frequency as seen by this query: dev/fit queries sit in the DB, so remove their own
    # contribution to reproduce the test-time situation (query absent from the DB).
    qinfo = queries.set_index("Entry")
    own = c["query"].map(qinfo.role).isin(["dev", "fit"]) & (c.ec == c["query"].map(qinfo.true_ec))
    c["log_train_freq"] = np.log1p(c.ec.map(ec_freq).fillna(0) - own.astype(int))

    c = c.sort_values(["query", "best_rank", "bits_share"], ascending=[True, True, False])
    c = c.groupby("query").head(max_candidates).copy()
    c["n_cands"] = c.groupby("query").ec.transform("size")

    truth = c["query"].map(qinfo.true_ec).str.split(";")
    c["label"] = [ec in t for ec, t in zip(c.ec, truth)]
    c["role"] = c["query"].map(qinfo.role)
    c["split"] = c["query"].map(qinfo.split)
    return c[["query", "role", "split", "ec", "label"] + FEATURES].reset_index(drop=True)
