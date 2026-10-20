# Decision log

Every choice that could change a result is recorded here **before** the result it
affects is produced. A threshold picked after seeing the output is a way of fitting
the answer; written first, it is a defensible method.

Append only. If a decision changes, add a new dated row and mark the old one
superseded — do not edit history.

## Decided

| # | Decision | Value | Why | Date |
| --- | --- | --- | --- | --- |
| D1 | Question being answered | Surface proteins up in tumors, absent from healthy tissue | the brief, restated | 2026-10-20 |
| D2 | Healthy reference | HPA normal tissue plus GTEx | He named HPA; GTEx adds samples and shares a processing pipeline with TCGA via the UCSC Xena Toil hub, which removes differences caused by software rather than biology. **GTEx is an addition to the sources he listed and must be stated to him.** | 2026-10-20 |
| D3 | Scope | Pan-cancer, reporting the best cancer type per gene | His wording says "in cancer" without naming one | 2026-10-20 |
| D4 | Language and environment manager | Python only; conda/mamba, env named `surfprot` | One language keeps the pipeline defensible by one person. Conda chosen over uv because Step 3 may need bioconda-only tools (DeepTMHMM, SignalP, a GPI predictor) that are awkward to install from PyPI. | 2026-10-20 |
| D5 | Python version | 3.12 | System Python is 3.14.7 (miniforge base). 3.14 wheel coverage is unproven for the single-cell extensions named at the end of the runbook; 3.12 is the safe floor for this stack. | 2026-10-20 |
| D6 | Raw data excluded from git | `data/raw/**` and `data/processed/**` gitignored; `results/**` committed | Datasets are gigabytes and are versioned by SHA-256 in `data/MANIFEST.md`. Results must be committed so every number in the report traces to a file under version control. | 2026-10-20 |

## Proposed, not yet decided

These appear in the runbook as the author's own assumptions, explicitly flagged there
as "not standards". **No code may use any of them until the row is moved to the
Decided table above with my confirmation and a date.**

| # | Proposal | Proposed value | Owned by | Status |
| --- | --- | --- | --- | --- |
| P1 | Minimum extracellular length for a surface call | 20 residues | Step 3 | Not decided. Reasoning is that a tiny stub is unlikely to hold an accessible binding site. Test 10 / 20 / 50 in Step 9. |
| P2 | How to run DeepTMHMM; which signal-peptide and GPI predictors | unknown | Step 3 | Not decided. Current access method unverified — check the official pages. |
| P3 | TCGA cancer type → normal tissue mapping | — | Step 4 | Not decided. A scientific judgment; must be reviewed by hand, not generated. |
| P4 | Minimum samples per group for a comparison | — | Step 4 | Not decided. |
| P5 | FDR threshold | < 0.01 | Step 4 | Not decided. Convention, not law. Test 0.05 / 0.01 in Step 9. |
| P6 | Log2 fold-change threshold | >= 2 (fourfold) | Step 4 | Not decided. Test 1 / 2 / 3 in Step 9. |
| P7 | Vital-organ list | heart, brain regions, lung, liver, kidney, bone marrow/blood | Step 5 | Not decided. Must be written before any specificity result is seen. Tolerated-tissue list decided alongside it, with reasons. |
| P8 | Minimum tumor expression | 5 TPM | Step 5 | Not decided. |
| P9 | Normal-tissue summary statistic | maximum, never the average | Step 5 | Not decided, but strongly argued: one exposed vital organ is enough to make a target unsafe. Test max vs 90th percentile in Step 9. |
| P10 | Primary statistical test | Wilcoxon rank-sum, with limma-voom as optional cross-check | Step 4 | Not decided. |
| P11 | What counts as "detected" in the GTEx proteome; CPTAC significance rule; treatment of unmeasured genes | — | Step 6 | Not decided. |
| P12 | Evidence-tier rules (tier 1 / 2 / 3) | per runbook Step 6 table | Step 6 | Not decided. |
| P13 | Selection rule and per-family limit | tier first, then score S; max 2 per gene family | Step 7 | Not decided. **Must be recorded here with a date earlier than the first ranking file.** |
| P14 | Source of positive and negative control lists | a published source, not memory | Step 8 | Not decided. |

## Superseded

_(none yet)_
