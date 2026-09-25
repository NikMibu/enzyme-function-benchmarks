from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
TOOLS = ROOT / "tools"
BENCH = ROOT / "benchmarks"
RESULTS = ROOT / "results"


def load(path: Path = ROOT / "configs" / "default.yaml") -> dict:
    return yaml.safe_load(open(path))
