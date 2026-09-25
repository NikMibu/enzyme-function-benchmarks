# enzyme-function-benchmarks

Experiments on how well an enzyme's function (its EC number) can be predicted: with homology
search, learned integrators, general-purpose decision models (TypeSafe Jev, Laya) and protein
language models (ESM-2, ProtT5). Both experiments are concluded.

**Overall findings: [`FINDINGS.md`](FINDINGS.md)**

## Results in brief

* **Performance comes from homology search.** Jev and Laya never do better than a simple method
  given the same information. On raw sequence they are at chance; with homolog evidence they tie
  the nearest neighbour.
* **Protein language models read function from sequence:** an ESM-2 probe reaches 51.9% on the EC
  class without a single detectable homolog in its training set.
* **The only significant gain over homology search:** ProtT5 embedding neighbours as a fallback for
  proteins without an MMseqs2 hit, 64.4% instead of 62.5% on the exact EC (p = 0.031).
* **Reliable confidence:** LightGBM over the homolog evidence roughly halves the area under the
  risk–coverage curve (0.092 vs 0.173), so it tells better when a prediction can be trusted.

## Contents

| Folder | Question | Data |
|---|---|---|
| [`enzyme-evidence/`](enzyme-evidence/) | Can a learned model weigh homolog evidence better than the nearest neighbour? | CARE benchmark, Task 1 (1,140 test proteins) |
| [`enzyme-direct/`](enzyme-direct/) | Can Jev and Laya read function from sequence, with context, with homologs? How do protein language models compare? | New benchmark from Swiss-Prot since 2018, separate from the CARE test sets (320 proteins) |

The two folders share no code or data. Each has its own README with the setup, all tables,
limitations and reproduction steps (`requirements.txt`, numbered scripts, `pytest -q`). Jev needs
`TYPESAFE_API_KEY`; everything else runs on CPU.

## Methodology

* Analysis plans and pass criteria are committed to git before each run; post-hoc analyses are
  labelled as such.
* Every stage has a control given exactly the same information; comparisons use paired tests.
* Leakage tests keep training and test data apart; benchmarks are frozen as files.
