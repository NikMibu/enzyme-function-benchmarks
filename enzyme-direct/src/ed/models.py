"""Model backends. Every call is cached as one JSONL line per protein, so a benchmark is only
ever paid for (and sent) once, and a crashed run resumes where it stopped."""
import json
import os
import time
from pathlib import Path

import requests


class Jev:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.key = os.environ.get("TYPESAFE_API_KEY")
        if not self.key:
            raise SystemExit("TYPESAFE_API_KEY is not set.")

    def __call__(self, state: str, questions: dict) -> dict:
        for attempt in range(6):
            r = requests.post(self.cfg["endpoint"], timeout=120,
                              headers={"Authorization": f"Bearer {self.key}"},
                              json={"model": self.cfg["model"], "state": state, "questions": questions})
            if r.status_code == 429 or r.status_code >= 500:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            return r.json()
        r.raise_for_status()


class Laya:
    def __init__(self, cfg: dict):
        from laya import Router
        self.router, self.model, self.max_len = Router(), cfg["router_model"], cfg["max_len"]

    def __call__(self, state: str, questions: dict) -> dict:
        resp = self.router.predict(state, questions, model=self.model, max_len=self.max_len)
        if resp["usage"]["input_tokens"] >= self.max_len:
            raise RuntimeError(f"input truncated at {self.max_len} tokens")
        return resp


def backend(name: str, cfg: dict):
    return Jev(cfg["models"][name]) if name == "jev" else Laya(cfg["models"][name])


def run(name: str, cfg: dict, items: list, cache: Path) -> dict:
    """items: [(entry, state, questions)] -> {entry: response}."""
    done = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            rec = json.loads(line)
            done[rec["entry"]] = rec["response"]
    todo = [it for it in items if it[0] not in done]
    if todo:
        call = backend(name, cfg)
        cache.parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "a") as f:
            for i, (entry, st, qs) in enumerate(todo, 1):
                resp = call(st, qs)
                done[entry] = resp
                f.write(json.dumps({"entry": entry, "response": resp}) + "\n")
                f.flush()
                if i % 25 == 0:
                    print(f"{name}: {i}/{len(todo)}", flush=True)
    return done
