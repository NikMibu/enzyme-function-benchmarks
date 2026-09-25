"""TypeSafe Jev as a zero-shot evidence integrator (choice over rendered candidates).

Needs TYPESAFE_API_KEY and network access to api.typesafe.ai. Responses are cached as JSONL so
a split is only ever paid for (and sent) once. The response is parsed as
answers.ec.probabilities (label -> p), falling back to a one-hot answers.ec.choice; this
matches the published request format but was not verified against a live key when written.
"""
import json
import os
import time
from pathlib import Path

import pandas as pd
import requests

from ..render import render


def _call(cfg: dict, payload: dict, key: str) -> dict:
    for attempt in range(5):
        r = requests.post(cfg["jev"]["endpoint"], timeout=60,
                          headers={"Authorization": f"Bearer {key}"},
                          json={"model": cfg["jev"]["model"], **payload})
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
        return r.json()
    r.raise_for_status()
    return r.json()


def _probs(resp: dict) -> dict:
    ans = resp["answers"]["ec"]
    if "probabilities" in ans:
        return {k: float(v) for k, v in ans["probabilities"].items()}
    return {ans["choice"]: 1.0}


def score(cfg: dict, c: pd.DataFrame, ec_names: dict, cache: Path) -> pd.Series:
    key = os.environ.get("TYPESAFE_API_KEY")
    if not key:
        raise SystemExit("TYPESAFE_API_KEY is not set; Jev integrator skipped.")
    done = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            rec = json.loads(line)
            done[rec["query"]] = rec
    out = pd.Series(0.0, index=c.index)
    with open(cache, "a") as f:
        for q, g in c.groupby("query"):
            req = render(q, g, ec_names)
            if q not in done:
                resp = _call(cfg, {"state": req["state"], "questions": req["questions"]}, key)
                done[q] = {"query": q, "label_to_ec": req["label_to_ec"], "response": resp}
                f.write(json.dumps(done[q]) + "\n")
            probs = _probs(done[q]["response"])
            ec_p = {ec: probs.get(lab, 0.0) for lab, ec in done[q]["label_to_ec"].items()}
            # Mass on the "none" option is abstention: it is simply not given to any candidate,
            # so it lowers the top candidate's score (= confidence).
            out[g.index] = g.ec.map(ec_p).fillna(0.0).values
    return out
