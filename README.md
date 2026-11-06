# Surface protein targets in cancer

A personal project, built in response to an open challenge.

**The question:** using Protein Atlas, TCGA and a TMHMM-type model, which cell-surface
proteins are significantly over-expressed in cancer and absent from healthy tissue?

Starting from the 20,431 reviewed human proteins, this narrows to ten candidates
by asking three things in order — **can a drug reach it**, **is it raised in
tumour**, **is it absent from healthy tissue everywhere else** — and then reports
how much the answer depends on the thresholds I chose.

## The answer

| Gene | Cancer type | S | Tumour TPM | Highest healthy tissue | TPM | Survives thresholds |
| --- | --- | --- | --- | --- | --- | --- |
| CLDN6 | Testicular germ cell | 6.40 | 158 | Brain – Cerebellum | 0.9 | 100% |
| TSPAN10 | Uveal melanoma | 6.20 | 221 | Testis | 2.0 | 99% |
| OR2I1 | Cholangiocarcinoma | 4.98 | 95 | Spleen | 2.1 | 62% |
| ALPG | Testicular germ cell | 4.71 | 49 | Fallopian tube | 0.9 | 84% |
| DCSTAMP | Thyroid carcinoma | 4.32 | 94 | Brain – Amygdala | 3.8 | 85% |
| TRPM1 | Uveal melanoma | 4.08 | 94 | Skin | 4.6 | 87% |
| UPK3B | Mesothelioma | 3.98 | 707 | Adipose – Visceral | 43.9 | 87% |
| LY6G6D | Rectum adenocarcinoma | 3.80 | 101 | Seminal vesicle | 6.3 | 79% |
| CD70 | Kidney clear cell | 3.65 | 52 | Lymph node | 3.2 | 50% |
| DBH | Pheochromocytoma | 3.47 | 2,859 | Adrenal gland | 257.4 | 75% |

**S** = log2(median tumour TPM + 1) − log2(highest healthy tissue TPM + 1). The
log of a ratio: how far the typical tumour sits above the single riskiest healthy
tissue. CLDN6 at 6.40 is roughly 85-fold.

**Survives thresholds** = the share of 270 threshold combinations in which the
gene still reaches the top 20 (see *Robustness*).

Four of the ten are in active clinical development, which the pipeline did not
know: CLDN6 (BNT211), CD70 (CTX130), ALPG (TCA-101) and LY6G6D (bispecific
programmes). Positive controls landed on clinically correct indications
unprompted — EGFR→glioblastoma, MSLN and FOLR1→ovarian, CD19→DLBCL.

## The funnel

| Stage | Genes |
| --- | --- |
| Reviewed human proteins (UniProt) | 20,431 |
| Mapped to an Ensembl gene | 19,400 |
| **Surface candidates** — UniProt topology + DeepTMHMM + GPI route | 4,141 |
| Present in the expression matrix | 3,912 |
| Pass tumour-vs-healthy RNA thresholds | 1,513 |
| Reach 5 TPM in tumour | 1,452 |
| **Eligible under the pre-registered rule** | 1,272 |
| **Ranked** | 10 + 5 reserve |

## How it works

**Surface call.** Two independent sources must agree on topology — UniProt's
curated annotation and DeepTMHMM run locally over the whole proteome. Agreement
gives tier S1, one source gives S2. 43% of UniProt's transmembrane proteins have
no topological domains annotated, which is why a predictor is needed rather than
curation alone. GPI-anchored proteins reach the surface by lipid anchor and have
no transmembrane segment, so they enter by a separate route (tier S3) — without
it, LY6G6D at rank 8 would be invisible.

**Tumour versus healthy.** TCGA tumours against GTEx healthy tissue, both
reprocessed through one pipeline by UCSC Toil — otherwise much of any measured
difference is the two pipelines disagreeing. Each of 32 cancer types is mapped to
the healthy tissue it arises from, graded `direct` / `imperfect` / `lineage_only`,
extending the published pairing of Wang et al. (2018). Wilcoxon per gene per
cancer type, BH-FDR within each comparison. Adjacent-normal tissue is run as a
separate second arm and never pooled with GTEx — Aran et al. (2017) show it is an
intermediate state, not healthy tissue.

**Absence from healthy tissue.** The maximum across 51 GTEx and 40 HPA tissues,
taking the higher of two independent atlases. Maximum rather than an average
because one exposed organ is enough: a 90th-percentile summary discards the top
five tissues of fifty-one, and would score mesothelin's worst tissue — lung — at
5.3 TPM instead of 58.9.

**Ranking.** The selection rule was written and dated *before* the first ranking
ran, and a check verifies that ordering. Eligibility is three conditions; the sort
key is S alone. Protein concordance, comparator grade, prevalence, trajectory and
the surfaceome atlases are attached to every candidate and **never** move the
rank — mass spectrometry under-detects membrane proteins, so letting it rank
would penalise targets for being hard to measure.

## Robustness

Five of the thresholds are conventions rather than measurements. The whole
selection was re-run under **all 270 combinations** of them — extracellular
minimum, FDR, fold change, healthy reference, healthy summary — and each gene
scored on how often it still reaches the top 20.

Seven of ten survive at least 75%. The choice that moves the answer most is **how
healthy tissue is summarised** (28% spread); **FDR changes it by 0%**, because at
these sample sizes 73% of all comparisons are already significant and the
fold-change threshold does the filtering.

## What this does not show

- **Tumour purity is not accounted for.** Bulk tissue mixes tumour cells with
  macrophages, T cells, fibroblasts and vessels. I attempted three ways to call
  this computationally and withdrew all three as weaker than they appeared; the
  limitation is stated rather than papered over with a flag.
- **The lineage confound.** S cannot separate a gene tumour cells genuinely
  overexpress from one whose cell of origin is rare in the healthy comparator.
  TSPAN10 and TRPM1 are melanocyte genes compared against skin, where melanocytes
  are a small minority. [An attempt to measure this failed](docs/decisions.md) —
  a marker composite measures how *differentiated* a lineage is, not how much is
  present, and tumours dedifferentiate. Reported as a negative result.
- **RNA is a proxy for protein**, and nine of the ten finalists sit in cancer
  types with no CPTAC proteomics cohort.
- **Eight of ten come from cohorts under 200 tumours.** A specificity score
  structurally favours rare, lineage-restricted tumours.
- **No benchmark against a published control list.** Deliberately not done — the
  source was never agreed, and building one from memory is exactly what the
  method forbids.
- **Nothing here is validated at the bench.**

## Repository

```
docs/METHODOLOGY.md      the plan, written before the work
docs/decisions.md        82 dated decisions — the backbone of the project
docs/checks.md           74 automated checks, one per claim
data/MANIFEST.md         44 files with version, URL and SHA-256
src/                     one script per step; every number is written by code
results/                 every output, including the ones that were withdrawn
results/surface_proteins_slides.pdf    slide walkthrough of the method
```

Everything in `results/` is written by a script except
[`surface_proteins_slides.pdf`](results/surface_proteins_slides.pdf), which is a
hand-made walkthrough of how the pipeline was built.

**`docs/decisions.md` is the file to read.** Every threshold, gene list and scope
choice is recorded with what, the value, why, and the date — written *before* the
result it affects, which is what makes a threshold a method rather than a way of
fitting the answer. It also records three analyses I built and then withdrew, and
the reasoning for each.

## Running it

```bash
conda env create -f environment.yml && conda activate surfprot
python src/step1_download.py --download      # 44 files, 1.7 GB, checksummed
python src/step2_build_gene_map.py
python src/step3_run_deeptmhmm.py            # local, batched, several hours
python src/step3_build_candidates.py
python src/step4_extract_matrix.py && python src/step4_differential.py
python src/step5_specificity.py
python src/step6_protein_differential.py && python src/step6_concordance.py
python src/step7_rank.py
python src/step9_robustness.py
```

Each step has a `*_check.py` beside it. Raw data is excluded from version control
and reproduced from `data/MANIFEST.md` by checksum.

## Still open

- **Finalist review** — the log is built (`docs/step10_override_log.md`) with HPA
  staining, UniProt and trial links for all 15, and the verdict column is mine to
  fill. Three candidates are queued for it with reasons on record: OR2I1 is the
  only one of 1,263 eligible genes UniProt calls existence-uncertain; DCSTAMP is
  present in 66% of tumours and is a myeloid gene, so infiltrate is the likelier
  confound; DBH has no plasma-membrane annotation.
- **Control benchmark** — see above.

## How this was built

A personal project, done on my own time and not on behalf of any institution.
No patient data was accessed — every dataset is public and listed with its
version and checksum in `data/MANIFEST.md`.

Written with Claude Code under my direction. Every threshold, gene list and scope
decision was proposed to me, decided by me, and recorded in `docs/decisions.md`
with a date before the result it affects existed. Every number in this README
comes from a script in `src/` that writes to `results/`; none is typed by hand.

The parts I would point to as mine are the judgement calls: rejecting a
microenvironment classifier whose threshold had been chosen by eye, declining to
author a vital-organ list I was not qualified to write, refusing to score DLBCL
against spleen because it is not the tissue that tumour comes from, and catching
a published claim that survived a string replacement which had silently matched
nothing.
