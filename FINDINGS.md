# Predicting enzyme function: overall findings

Two concluded experiments, 25 September 2026:

* [`enzyme-evidence/`](enzyme-evidence/): weighing homolog evidence on the CARE benchmark (Task 1)
* [`enzyme-direct/`](enzyme-direct/): function from sequence, with computed context, with homologs and with protein embeddings

Setup details, full tables and reproduction steps are in the two folders' READMEs. Every number
here comes from their `results/` files.

## Key finding

**Performance comes from homology search, not from the decision model.** Jev (TypeSafe,
`jev-1.13.0`) and Laya (open weights) never do better than a simple method given the same
information. The only significant gain over plain homology search comes from a protein language
model used as a fallback for proteins without a hit.

## enzyme-evidence: weighing homologs

CARE test sets (1,140 proteins), exact EC number, MMseqs2 search against CARE's training set.

| Method | Accuracy, pooled | Selective prediction (AURC, pooled, lower is better) |
|---|---|---|
| Nearest neighbour (MMseqs2) | 69.5% | 0.173 |
| LightGBM over 19 evidence features | 69.8% | **0.092** |
| Jev, zero-shot | 69.1% | 0.221 |

* No method beats the nearest neighbour on accuracy, although about 8 points of headroom exist.
* The one real gain is LightGBM's confidence. On the hardest split (<30% identity), answering the
  half it is most confident about, it is right 91.7% of the time, against 80.6% when
  thresholding on identity × coverage. Jev's confidence is worse than identity × coverage on that
  split (AURC 0.269 vs 0.168).
* The MMseqs2 search alone beats CARE's published baselines (62.5% vs 55.1% for CLEAN on the
  <30% split).

## enzyme-direct: function without annotations

A new benchmark of Swiss-Prot entries created since 2018, separate from the CARE test sets, one
protein per sequence family: ec1 (6 EC classes × 35, chance 16.7%) and ec4 (11 exact ECs × 10).

| What the method gets | Jev | Best comparison method |
|---|---|---|
| sequence only | 11.9% (Laya 15.7–16.2%) | ESM-2 probe: **51.9%** |
| + properties & composition | 14.3% | Logistic regression: 29.0% |
| + motifs | 27.1% | Logistic regression: 32.9% |
| + homolog evidence | 87.1% | Nearest neighbour: 86.7% |

EC class, ec1 (210 proteins).

* **Raw sequence:** Jev and Laya are at chance and give nearly the same class for every protein.
* **Context:** Jev only uses clues that name a chemistry (motifs), and even then stays below
  logistic regression.
* **Homologs:** Jev departs from the nearest neighbour on only 19 of 210 proteins, and loses 6 to 1
  on the exact EC among those.
* **Protein language models read function from sequence:** an ESM-2 (650M) linear probe reaches
  51.9% without a single detectable homolog in its training set.
* **Embedding fallback:** when MMseqs2 finds no hit, use the ProtT5 nearest neighbour from
  UniProt's precomputed embeddings. On ec1 + ec4 (320 proteins):

  | | MMseqs2 nearest neighbour | With embedding fallback |
  |---|---|---|
  | Exact EC | 62.5% | **64.4%** (6 gained, 0 lost, p = 0.031) |
  | EC class | 80.3% | **87.5%** (p = 2×10⁻⁷) |

## What follows

1. **Jev is a good reader, not a reasoner.** Given evidence laid out in words, it picks the
   obvious answer. It draws no conclusions beyond the stated clues, and its probabilities are
   poorly calibrated. It is built for text decisions, not for proteins.
2. **Match the tool to the data.** Homology search handles relatives; protein language models
   handle sequences without relatives. A text model reading amino-acid strings is the wrong
   tool, which is why fine-tuning Laya was not pursued.
3. **The value lies in two niches:** reliable confidence (LightGBM) and proteins without a homolog
   (embedding fallback). Both are small, targeted improvements on top of homology search.

## Methodology

**Strengths**

* Analysis plans and pass criteria were committed to git before each run (enzyme-direct: commits
  `1460eab`, `35e770b`, `987ac48`, `16e8082`, `ed4b12f`). Post-hoc analyses are labelled as such.
* Each Jev level has a control given exactly the same information; comparisons use paired tests
  on the same proteins.
* Leakage tests: no overlap with the CARE test sets, and no accessions, names or organisms in
  prompts. Benchmarks are frozen as files.

**Limits**

* **Small samples:** 210 and 320 proteins, 37 without a hit. The embedding gain rests on 6
  proteins and should be confirmed on new data.
* **Design changes before the first model run:** only CARE's test sets were excluded (not its
  training set), 35 instead of 40 proteins per class, and entries created since 2018 instead of
  2020. This was needed because otherwise too few EC 6 families remained, and it happened before
  any model was run.
* **Prompts:** Jev ran with one fixed prompt per level, deliberately untuned on test data. Better
  prompts are conceivable but would hardly close the large gaps.
* **Weak controls:** the first logistic-regression control had only 4 ligases in training, so
  "Jev does not beat its control" is a conservative statement.
* **Embeddings:** ProtT5 was trained on UniRef50 and may have seen these sequences, though without
  function labels.

## Open questions

* Confirm the embedding fallback on a larger benchmark with many proteins without a hit.
* Cosine similarity as confidence for fallback answers (median 0.72 without a hit, 0.96 with one).
* Test LightGBM's confidence from enzyme-evidence on the new benchmark.
