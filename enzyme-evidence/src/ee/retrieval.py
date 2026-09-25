"""MMseqs2 search of every query against the CARE training set (the only retrieval database)."""
import platform
import shutil
import subprocess
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

from . import config

FORMAT = "query,target,fident,alnlen,qlen,tlen,qcov,tcov,evalue,bits"


def mmseqs_bin(cfg) -> Path:
    exe = config.TOOLS / "mmseqs" / "bin" / "mmseqs"
    if exe.exists():
        return exe
    if shutil.which("mmseqs"):
        return Path(shutil.which("mmseqs"))
    assert platform.system() == "Linux", "auto-download only supports Linux; install mmseqs2 yourself"
    config.TOOLS.mkdir(parents=True, exist_ok=True)
    rel = cfg["retrieval"]["mmseqs_release"]
    url = f"https://github.com/soedinglab/MMseqs2/releases/download/{rel}/mmseqs-linux-avx2.tar.gz"
    tgz = config.TOOLS / "mmseqs.tar.gz"
    urllib.request.urlretrieve(url, tgz)  # nosec - pinned release
    with tarfile.open(tgz) as f:
        f.extractall(config.TOOLS)  # nosec
    tgz.unlink()
    return exe


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    with open(path, "w") as f:
        for entry, seq in zip(df.Entry, df.Sequence):
            f.write(f">{entry}\n{seq}\n")


def search(cfg, queries: pd.DataFrame, train: pd.DataFrame, workdir: Path) -> pd.DataFrame:
    r = cfg["retrieval"]
    exe = str(mmseqs_bin(cfg))
    workdir.mkdir(parents=True, exist_ok=True)
    write_fasta(train, workdir / "train.fasta")
    write_fasta(queries, workdir / "queries.fasta")

    def run(*args):
        subprocess.run([exe, *map(str, args)], check=True, stdout=subprocess.DEVNULL)

    run("createdb", workdir / "train.fasta", workdir / "trainDB")
    run("createdb", workdir / "queries.fasta", workdir / "queryDB")
    run("search", workdir / "queryDB", workdir / "trainDB", workdir / "res", workdir / "tmp",
        "-s", r["sensitivity"], "--max-seqs", r["max_seqs"], "-e", r["evalue"],
        "--threads", r["threads"])
    run("convertalis", workdir / "queryDB", workdir / "trainDB", workdir / "res",
        workdir / "hits.tsv", "--format-output", FORMAT)
    hits = pd.read_csv(workdir / "hits.tsv", sep="\t", header=None, names=FORMAT.split(","))
    return hits
