"""Build the held-out benchmarks from recent Swiss-Prot entries that are absent from the CARE test sets.

    python scripts/01_build_benchmarks.py

ec1: balanced EC level-1 set, one protein per 30%-identity cluster.
ec4: 10-20 well-represented exact EC numbers, proteins disjoint from ec1.
Both are written to benchmarks/ and committed, so later UniProt releases do not change them.
"""
import _common  # noqa: F401
import json

import pandas as pd

from ed import config, data

cfg = config.load()
care_acc, care_seq = data.care_exclusion(cfg)
raw, release = data.fetch_uniprot(cfg)
pool, log = data.clean_pool(raw, care_acc, care_seq)
log["uniprot_release"] = release
log["care_accessions_excluded"], log["care_sequences_excluded"] = len(care_acc), len(care_seq)
config.BENCH.mkdir(exist_ok=True)
cols = ["Entry", "ec", "ec1", "Length", "Date of creation", "Organism", "Protein names", "cluster", "Sequence"]

# ---- ec1
s1 = cfg["stages"]["ec1"]
cl30 = data.cluster(cfg, pool, cfg["mmseqs"]["cluster_min_seq_id"], "pool30")
log["ec1_pool_by_class"] = pool.ec1.value_counts().sort_index().to_dict()
log["ec1_clusters_by_class"] = (pool.assign(c=pool.Entry.map(cl30)).groupby("ec1").c.nunique().to_dict())
ec1 = data.pick_one_per_cluster(pool, cl30, "ec1", s1["classes"], s1["per_class"], cfg["seed"])
ec1[cols].to_csv(config.BENCH / "ec1.tsv", sep="\t", index=False)

# ---- ec4 (no protein, and no 30% cluster, shared with ec1)
s4 = cfg["stages"]["ec4"]
rest = pool[~pool.Entry.isin(ec1.Entry) & ~pool.Entry.map(cl30).isin(set(ec1.cluster))].reset_index(drop=True)
cl50 = data.cluster(cfg, rest, s4["cluster_min_seq_id"], "rest50")
n_clusters = rest.assign(c=rest.Entry.map(cl50)).groupby("ec").c.nunique().sort_values(ascending=False)
eligible = n_clusters[n_clusters >= s4["per_ec"]]
ecs = sorted(eligible.index[: s4["n_ecs"]])
if len(ecs) < 10:
    raise SystemExit(f"only {len(ecs)} ECs have {s4['per_ec']} distinct clusters")
log["ec4_ecs"] = {ec: int(n_clusters[ec]) for ec in ecs}
ec4 = data.pick_one_per_cluster(rest, cl50, "ec", ecs, s4["per_ec"], cfg["seed"])
ec4[cols].to_csv(config.BENCH / "ec4.tsv", sep="\t", index=False)

names, log["enzyme_release"] = data.ec_names(cfg)
missing = [ec for ec in ecs if ec not in names]
if missing:
    raise SystemExit(f"no ENZYME name for {missing}")
pd.DataFrame({"ec": ecs, "name": [names[ec] for ec in ecs]}).to_csv(
    config.BENCH / "ec4_labels.tsv", sep="\t", index=False)
log["ec4_names"] = {ec: names[ec] for ec in ecs}
log["n_ec1"], log["n_ec4"] = len(ec1), len(ec4)
json.dump(log, open(config.BENCH / "build_log.json", "w"), indent=2)
print(json.dumps(log, indent=1))
