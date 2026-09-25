#!/usr/bin/env bash
# Clone CARE at the pinned commit into data/raw/CARE.
set -euo pipefail
cd "$(dirname "$0")/.."
COMMIT=$(python3 -c "import yaml;print(yaml.safe_load(open('configs/default.yaml'))['care']['commit'])")
if [ ! -d data/raw/CARE/.git ]; then
  mkdir -p data/raw
  git clone --quiet https://github.com/jsunn-y/CARE data/raw/CARE
fi
git -C data/raw/CARE checkout --quiet "$COMMIT"
echo "CARE at $(git -C data/raw/CARE rev-parse HEAD)"
