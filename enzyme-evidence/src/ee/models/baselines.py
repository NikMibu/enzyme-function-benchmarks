"""Non-learned integrators over the same candidate table."""
import pandas as pd


def nearest_neighbour(c: pd.DataFrame) -> pd.Series:
    """EC of the top hit. Ties within a multi-EC top hit break on vote share, then EC frequency.

    Score = top-hit identity x query coverage (the 'Conf-0' confidence), identical for every
    candidate of a query except that non-top candidates are pushed below it.
    """
    order = (c.assign(_neg_rank=-c.best_rank)
              .sort_values(["query", "_neg_rank", "bits_share", "log_train_freq"],
                           ascending=[True, False, False, False]))
    is_top = ~order["query"].duplicated()
    conf = (c.top_fid * c.top_qcov).clip(0, 1)
    score = pd.Series(0.0, index=c.index)
    score[order.index[is_top]] = conf[order.index[is_top]]
    score[order.index[~is_top]] = conf[order.index[~is_top]] * 0.5 / c.loc[order.index[~is_top], "best_rank"]
    return score


def weighted_vote(c: pd.DataFrame) -> pd.Series:
    """Bitscore-weighted vote over the top-k neighbours; score = the candidate's vote share."""
    return c.bits_share.clip(0, 1)


def oracle(c: pd.DataFrame) -> pd.Series:
    """Upper bound for any re-ranker over this candidate set."""
    return c.label.astype(float)
