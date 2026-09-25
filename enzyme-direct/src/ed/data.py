"""Held-out benchmark construction: recent Swiss-Prot enzymes absent from the CARE test sets."""
import io
import platform
import random
import re
import subprocess
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd
import requests

from . import config

UNIPROT = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,id,ec,length,date_created,organism_name,protein_name,sequence"
FULL_EC = re.compile(r"^\d+\.\d+\.\d+\.\d+$")  # excludes partial (1.1.1.-) and preliminary (n) ECs
STANDARD_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


# ---------------------------------------------------------------- CARE (exclusion only)
def fetch_care(cfg) -> Path:
    c = cfg["care"]
    out = config.DATA / "raw" / "CARE" / c["commit"]
    for rel in c["files"]:
        dest = out / rel
        if not dest.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(f"{c['raw_base']}/{c['commit']}/{rel}", dest)  # nosec - pinned
    return out


def care_exclusion(cfg) -> tuple[set, set]:
    """Accessions and exact sequences in the CARE test splits."""
    root = fetch_care(cfg)
    acc, seqs = set(), set()
    for rel in cfg["care"]["files"]:
        if "text2EC" in rel:
            continue
        df = pd.read_csv(root / rel, usecols=["Entry", "Sequence"])
        acc |= set(df.Entry)
        seqs |= set(df.Sequence)
    return acc, seqs


def ec_names(cfg) -> tuple[dict, str]:
    """Accepted EC names from the ExPASy ENZYME database, and its release line."""
    path = config.DATA / "raw" / "enzyme.dat"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(cfg["enzyme_db"], path)  # nosec - fixed public URL
    names, ec, release = {}, None, ""
    for line in open(path):
        if line.startswith("CC   Release"):
            release = line[5:].strip()
        elif line.startswith("ID   "):
            ec = line[5:].strip()
        elif line.startswith("DE   ") and ec:
            names[ec] = names.get(ec, "") + line[5:].strip()
        elif line.startswith("//"):
            ec = None
    return {k: v.rstrip(".") for k, v in names.items()}, release


# ---------------------------------------------------------------- UniProt
def fetch_uniprot(cfg) -> tuple[pd.DataFrame, str]:
    u = cfg["uniprot"]
    cache = config.DATA / "raw" / "uniprot_pool.tsv"
    meta = config.DATA / "raw" / "uniprot_release.txt"
    if cache.exists():
        return pd.read_csv(cache, sep="\t"), meta.read_text().strip()
    frames, release = [], ""
    for cls in range(1, 8):
        q = (f"reviewed:true AND fragment:false AND ec:{cls}.* "
             f"AND date_created:[{u['created_after']} TO *] "
             f"AND length:[{u['min_len']} TO {u['max_len']}]")
        r = requests.get(UNIPROT, params={"query": q, "format": "tsv", "fields": FIELDS}, timeout=600)
        r.raise_for_status()
        release = f"{r.headers.get('X-UniProt-Release', '')} {r.headers.get('X-UniProt-Release-Date', '')}"
        frames.append(pd.read_csv(io.StringIO(r.text), sep="\t"))
    pool = pd.concat(frames).drop_duplicates("Entry").reset_index(drop=True)
    cache.parent.mkdir(parents=True, exist_ok=True)
    pool.to_csv(cache, sep="\t", index=False)
    meta.write_text(release.strip())
    return pool, release.strip()


def clean_pool(pool: pd.DataFrame, care_acc: set, care_seq: set) -> tuple[pd.DataFrame, dict]:
    log = {"uniprot_rows": len(pool)}
    ec = pool["EC number"].fillna("").astype(str)
    single = ec.str.split("; ").map(len).eq(1) & ec.map(lambda e: bool(FULL_EC.match(e)))
    p = pool[single].copy()
    log["single_full_ec"] = len(p)
    p = p[p.Sequence.map(lambda s: bool(STANDARD_AA.match(s)))]
    log["standard_residues"] = len(p)
    p = p[~p.Entry.isin(care_acc) & ~p.Sequence.isin(care_seq)]
    log["not_in_care"] = len(p)
    p = p.drop_duplicates("Sequence")
    log["unique_sequence"] = len(p)
    p["ec"] = p["EC number"]
    p["ec1"] = p.ec.str.split(".").str[0]
    return p.reset_index(drop=True), log


# ---------------------------------------------------------------- MMseqs2 clustering
def mmseqs_bin(cfg) -> Path:
    exe = config.TOOLS / "mmseqs" / "bin" / "mmseqs"
    if exe.exists():
        return exe
    assert platform.system() == "Linux", "auto-download only supports Linux; install mmseqs2 yourself"
    config.TOOLS.mkdir(parents=True, exist_ok=True)
    url = (f"https://github.com/soedinglab/MMseqs2/releases/download/"
           f"{cfg['mmseqs']['release']}/mmseqs-linux-avx2.tar.gz")
    tgz = config.TOOLS / "mmseqs.tar.gz"
    urllib.request.urlretrieve(url, tgz)  # nosec - pinned release
    with tarfile.open(tgz) as f:
        f.extractall(config.TOOLS)  # nosec
    tgz.unlink()
    return exe


def cluster(cfg, df: pd.DataFrame, min_seq_id: float, tag: str) -> pd.Series:
    """Entry -> cluster representative, from mmseqs easy-cluster."""
    work = config.DATA / "cluster" / tag
    work.mkdir(parents=True, exist_ok=True)
    with open(work / "in.fasta", "w") as f:
        for e, s in zip(df.Entry, df.Sequence):
            f.write(f">{e}\n{s}\n")
    subprocess.run([str(mmseqs_bin(cfg)), "easy-cluster", work / "in.fasta", work / "out", work / "tmp",
                    "--min-seq-id", str(min_seq_id), "-c", str(cfg["mmseqs"]["cluster_coverage"]),
                    "--cov-mode", "0", "--threads", "4"], check=True, stdout=subprocess.DEVNULL)
    tab = pd.read_csv(work / "out_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    return tab.set_index("member").rep


def pick_one_per_cluster(df: pd.DataFrame, clusters: pd.Series, group: str, groups: list,
                         n: int, seed: int) -> pd.DataFrame:
    """For each group, draw up to n proteins in a seeded random order, never two from one cluster
    (clusters are shared across groups, so a family cannot appear under two labels either)."""
    order = list(df.index)
    random.Random(seed).shuffle(order)
    used, picked = set(), {g: [] for g in groups}
    for i in order:
        g, c = df.at[i, group], clusters[df.at[i, "Entry"]]
        if g in picked and len(picked[g]) < n and c not in used:
            picked[g].append(i)
            used.add(c)
    short = {g: len(v) for g, v in picked.items() if len(v) < n}
    if short:
        raise SystemExit(f"not enough clusters for {short}")
    out = df.loc[[i for g in groups for i in picked[g]]].copy()
    out["cluster"] = out.Entry.map(clusters)
    return out.reset_index(drop=True)
