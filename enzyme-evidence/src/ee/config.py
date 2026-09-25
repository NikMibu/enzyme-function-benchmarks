from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
DERIVED = DATA / "derived"
RESULTS = ROOT / "results"
TOOLS = ROOT / "tools"
CONFIG = ROOT / "configs" / "default.yaml"


def load(path: Path = CONFIG) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def save(cfg: dict, path: Path = CONFIG) -> None:
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)


def care_dir(cfg: dict) -> Path:
    return ROOT / cfg["care"]["dir"]
