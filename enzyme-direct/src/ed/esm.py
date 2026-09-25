"""Mean-pooled ESM-2 embeddings, cached per protein under data/esm/."""
import numpy as np
import torch

from . import config


class Embedder:
    def __init__(self, name: str):
        from transformers import AutoModel, AutoTokenizer
        torch.set_num_threads(4)
        self.tok = AutoTokenizer.from_pretrained(name)
        self.model = AutoModel.from_pretrained(name).eval()
        self.cache = config.DATA / "esm" / name.split("/")[-1]
        self.cache.mkdir(parents=True, exist_ok=True)

    def __call__(self, entries, seqs) -> np.ndarray:
        out = []
        todo = sum(not (self.cache / f"{e}.npy").exists() for e in entries)
        done = 0
        for e, s in zip(entries, seqs):
            f = self.cache / f"{e}.npy"
            if not f.exists():
                with torch.inference_mode():
                    h = self.model(**self.tok(s, return_tensors="pt")).last_hidden_state[0, 1:-1]
                np.save(f, h.mean(0).numpy().astype(np.float32))
                done += 1
                if done % 100 == 0:
                    print(f"esm: {done}/{todo}", flush=True)
            out.append(np.load(f))
        return np.stack(out)
