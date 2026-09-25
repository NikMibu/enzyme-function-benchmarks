import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ee import config  # noqa: E402


@pytest.fixture(scope="session")
def cfg():
    c = config.load()
    if not (config.care_dir(c) / "splits").exists():
        pytest.skip("CARE data not fetched (run scripts/00_fetch_care.sh)")
    return c
