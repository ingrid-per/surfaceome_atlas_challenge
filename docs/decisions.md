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
| D7 | Exact dataset file set and releases | 40 files, ~1.64 GB — see `data/MANIFEST.md` for URLs, releases and checksums | Resolved by probing the live endpoints, not from memory. Releases pinned: UniProt 2026_03 (02-Sep-2026, 20,431 reviewed human entries); HPA 25.1 on Ensembl 109; Xena Toil files dated 2021-04-09 on GENCODE v23; LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized; GTEx proteome Jiang et al. Cell 2020; CSPA 2015; SURFY 2018. | 2026-10-20 |
| D8 | HPA cancer IHC (pathology) excluded | Not downloaded | HPA 25.1 no longer publishes it; it exists only on the archived v23 site. Mixing release 23 into an otherwise 25.1 analysis is worse than omitting it. Nothing in Steps 1–9 needs it — tumor evidence comes from TCGA RNA (Step 4) and CPTAC protein (Step 6) — and Step 10 reviews tumor staining on the HPA website by eye. | 2026-10-20 |
| D9 | Source of the GTEx arm in Step 5 | The 7,429 GTEx `Normal Tissue` samples inside the Toil matrix, **not** HPA's republished GTEx | D2 takes the higher of the HPA and GTEx values as the cautious safety estimate. That is only a real cross-check if the two come from independent processing pipelines. HPA's `rna_tissue_gtex` is GTEx rescaled into HPA's own nTPM, so using it would make both arms HPA-processed and the rule meaningless. Toil also puts GTEx on the same pipeline as the TCGA tumors, which is why D2 chose Toil in the first place. | 2026-10-20 |
| D10 | Samples included from the Toil matrix | TCGA `Primary Tumor` and `Solid Tissue Normal`; GTEx `Normal Tissue`. Excluded: all TARGET samples, the 433 GTEx cell lines, TCGA `Metastatic` and `Recurrent Tumor`, and blood-derived cancers | TARGET is pediatric cancer with different biology and was not asked about. A cell line is not a tissue and cannot speak to on-target toxicity in a patient. Metastatic and recurrent samples would mix disease states into a primary-tumor comparison. | 2026-10-20 |
| D11 | Unit of the Toil expression matrix | Values are `log2(TPM + 0.001)`, not TPM. Convert with `TPM = 2^x − 0.001`, clipped at 0, before applying Step 4's `log2(TPM + 1)` | Verified by inspecting the matrix: the floor value `-9.9658` is exactly `log2(0.001)`. Treating these values as TPM, or as already-log2(TPM+1), would silently corrupt every fold change in Step 4. | 2026-10-20 |
| D12 | CPTAC cohorts usable for tumor-vs-normal protein | 8 of 10: CCRCC, COAD, HNSCC, LSCC, LUAD, OV, PDAC, UCEC | BRCA and GBM have no normal proteomics file in LinkedOmicsKB freeze v1.2_reorganized. They still contribute tumor protein evidence, but must be labelled as having no matched normal rather than counted as "no protein in normal tissue". | 2026-10-20 |
| D13 | CPTAC tumor-vs-normal computed by us, not taken precomputed | Compute in Step 6 from the per-cohort abundance matrices | LinkedOmicsKB ships no precomputed differential results. Computing it ourselves applies one consistent method across all 8 cohorts. | 2026-10-20 |
| D14 | Gene universe (the "spine") | All 20,431 reviewed UniProt human accessions | Step 3 cannot call topology without a curated sequence, and only UniProt guarantees one. Genes drop out at each join instead of being filtered away up front, so the Step 6 funnel shows exactly where each loss happened rather than hiding it in a pre-filter. | 2026-10-20 |
| D15 | Resolving one-to-many ID matches | Keep every pairing; log it in `results/step2_multimapped.csv`; flag the row `ambiguous_mapping` | 1,615 accessions point at more than one Ensembl gene and 199 Ensembl genes at more than one accession. Dropping them would discard real surface candidates invisibly; picking a canonical one automatically would be a judgment that has to be defended gene by gene. Keeping and flagging leaves the ambiguity visible all the way into the final table, where Step 10 can weigh it. | 2026-10-20 |
| D16 | Unreviewed (TrEMBL) UniProt entries | Excluded | TrEMBL entries carry automatic rather than curated annotation, so they would enter Step 3 with predictor evidence only, losing half of the two-source agreement the tiers are built on. Including ~80,000 more entries would also push the candidate count far outside the 1,500–4,500 range the runbook uses as a sanity check, weakening that check. | 2026-10-20 |
| D17 | Join key across all datasets | Version-stripped Ensembl gene ID (`ENSG00000146648.22` → `ENSG00000146648`). Gene symbols are carried for reading and **never** used as a join key | Symbols are renamed over time and the datasets are from different years. Concrete case found while planning: NECTIN4, the target of an approved antibody-drug conjugate, appears absent from the Xena matrix when matched by symbol, because GENCODE v23 (2015) still calls it `PVRL4`. Its gene ID `ENSG00000143217` is present and matches correctly. A symbol join would have dropped a known-good positive control and made Step 8 report a false failure. PVRL1/2/3 have the same problem. | 2026-10-20 |
| D18 | UniProt Ensembl cross-reference field | Re-pull using `xref_ensembl_full` as resource U3; keep U2 for its topology annotation | The field pinned in Step 1 (`xref_ensembl`) returns transcript IDs (`ENST…`), not gene IDs, so U2 contains no gene-level link at all. `xref_ensembl_full` returns `ENSG…`. U2 is otherwise correct and still supplies the topology columns Step 3 needs. | 2026-10-20 |
| D19 | Genes present in HPA/UniProt but absent from the Xena matrix | Kept in the gene map, flagged `in_xena = False` | 374 HPA genes have no Xena row because Xena is built on GENCODE v23 (2015) while HPA is on Ensembl 109. The skew cannot be fixed without reprocessing TCGA. Flagging rather than dropping means Step 4 reports these as "no expression data" instead of letting them look like genes that were tested and failed. | 2026-10-20 |

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
