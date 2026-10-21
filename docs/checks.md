# Check register

One row per automated check. Every step adds its checks here before it is called
done. A check that cannot be run by a command is not a check — it is a hope.

Run everything: `pytest -q` from the project root (once tests exist), plus the
step scripts listed below.

| # | Step | Check | Command | Passes when |
| --- | --- | --- | --- | --- |
| 0.1 | 0 | Project layout, rule files and gitignore are complete | `python src/step0_verify_layout.py` | Exit 0 |
| 0.2 | 0 | Machine is capable and all six data hosts are reachable | `python src/step0_check_machine.py` | Exit 0; report written to `docs/step0_machine_report.md` |
| 0.3 | 0 | Environment solves and the full stack imports | `conda run -n surfprot python -c "import pandas, numpy, scipy, statsmodels, matplotlib, Bio, requests, pyarrow"` | Exit 0 |
| 0.4 | 0 | Git initialised, one commit, raw data excluded | `git log --oneline` and `git check-ignore -v data/raw/probe.csv` | One commit listed; probe path reported as ignored |
| 1.0 | 1 | Nothing downloads before the file list is seen | `python src/step1_download.py` (dry run is the default) | Prints 40 files with live sizes, downloads nothing |
| 1.1 | 1 | Every pinned file is on disk, in the manifest, with a SHA-256 and matching size; no undeclared file in `data/raw` | `python src/step1_check.py` | Exit 0 |
| 1.2 | 1 | Dataset shapes and sample counts are plausible; no group is zero | `python src/step1_inventory.py` then `python src/step1_check.py` | Toil = 60,498 genes × 19,131 samples; TCGA tumors ≥ 9,000; GTEx normals ≥ 7,000 |
| 1.3 | 1 | A moved URL cannot silently produce an HTML file named `.xlsx` | built into `step1_download.py` | Magic bytes checked before a download is accepted |
| 1.4 | 1 | `data/raw` and its contents are read-only | `python src/step1_check.py` | Exit 0 |
| 2.1 | 2 | Mapping coverage printed; a large loss is reported, not hidden | `python src/step2_check.py` | ≥90% of the 20,431 accessions carry an Ensembl gene |
| 2.2 | 2 | Every ambiguous mapping appears in the log | `python src/step2_check.py` | 0 unexplained ambiguous accessions |
| 2.3 | 2 | Spot-check of 11 known genes resolves correctly | `python src/step2_check.py` | EGFR, ERBB2, MSLN, TACSTD2, GAPDH, ACTB, B2M, FOLR1, NECTIN4, CLDN6, CD19 all match HPA and Xena |
| 2.4 | 2 | Renamed-gene regression — a symbol join would lose these | `python src/step2_check.py` | NECTIN1–4 all match Xena under their legacy PVRL1–4 symbols |
| 2.5 | 2 | No accession lost silently | `python src/step2_check.py` | mapped + unmapped = 20,431 |
| 2.6 | 2 | `gene_map.csv` is deterministic | `python src/step2_check.py` | Two runs byte-identical |
| 3.1 | 3 | Every reviewed protein accounted for | `python src/step3_check.py` | 20,431 rows |
| 3.2 | 3 | Candidate count in the expected range | `python src/step3_check.py` | 1,500–4,500 (SURFY predicts 2,886; CSPA measured 1,492) |
| 3.3 | 3 | Positive controls reach S1 | `python src/step3_check.py` | EGFR, ERBB2, TACSTD2, NECTIN4, CLDN6, CD19 all S1 with two sources |
| 3.4 | 3 | GPI controls reach S3 with no numeric length | `python src/step3_check.py` | MSLN, FOLR1 = S3, extracellular `not_applicable` (D20) |
| 3.5 | 3 | Negative controls excluded | `python src/step3_check.py` | GAPDH, ACTB, H3C1, H4C1, CALM1 not candidates |
| 3.6 | 3 | Excluded-by-design proteins documented | `python src/step3_check.py` | B2M absent, with the reason stated (see P6b) |
| 3.7 | 3 | S1 always has two independent sources | `python src/step3_check.py` | 0 S1 rows with `tm_sources < 2` or a missing prediction |
| 3.8 | 3 | Every GPI call records its source | `python src/step3_check.py` | 0 GPI calls without `gpi_source`; predictor-only count reported (D22) |
| 3.9 | 3 | Predictor gaps flagged, never silently negative | `python src/step3_check.py` | No protein reaches S1 without a prediction |
| 3.10 | 3 | Compartment filter applied and reported (D31) | `python src/step3_build_candidates.py` | Internal-membrane-only candidates excluded with evidence; unannotated kept |
| 3.11 | 3 | A curated GPI anchor is not overridden by a predicted TM helix (D32) | `python src/step3_check.py` | Every curated-GPI protein with no UniProt TM is S3; CEACAM5, CD48, FCGR3B, TNFRSF10C, ART1 named explicitly |
| 4.1 | 4 | Back-transform correct; no negative expression (D11) | `python src/step4_check.py` | Matrix floor −9.9658 maps to exactly 0; 3,912 × 19,131 |
| 4.2 | 4 | No comparison ran below the sample minimum (D37) | `python src/step4_check.py` | 0 violations; skipped and low-power comparisons listed |
| 4.3 | 4 | The heavy p<0.01 mass is power, not a batch effect | `python src/step4_check.py` | Significance tracks sample size (ρ≥0.5); dependence gone after the effect-size filter; global median log2FC ≈ 0; both directions represented |
| 4.4 | 4 | FDR computed within each (cancer type, arm) group (D39) | `python src/step4_check.py` | Independent recomputation of one group matches |
| 4.5 | 4 | Arm concordance reported (D34) | `python src/step4_check.py` | ρ and sign agreement printed; fails if ρ < 0.2 |
| 4.6 | 4 | Positive controls pass somewhere | `python src/step4_check.py` | EGFR, ERBB2, TACSTD2, NECTIN4, MSLN, CLDN6, CD19, FOLR1, CEACAM5 each pass in ≥1 cancer type |
| 4.7 | 4 | Differential results are deterministic | `python src/step4_check.py` | Two runs byte-identical |
| 4.8 | 4 | Three-state trajectory reported (D41) | `python src/step4_trajectory.py` | Every gene × cancer-type pair with all three states classified; hits split into tumour-specific vs field-effect |
| 5.1 | 5 | Every carried-forward row has a healthy maximum and named tissue; 5 TPM floor honoured (D43) | `python src/step5_check.py` | 0 missing; no row below 5 TPM |
| 5.2 | 5 | The maximum is never below the 90th percentile | `python src/step5_check.py` | 0 violations |
| 5.3 | 5 | The rule-selected ubiquitous control scores near zero (D47, P6b) | `python src/step5_check.py` | Top-5 control genes all score S ≤ +2 |
| 5.4 | 5 | Known targets characterised, never used to tune the score (D45) | `python src/step5_check.py` | All 11 reported with S and S_relaxed; no pass/fail applied |
| 5.5 | 5 | Tissue-exclusion variants ordered S ≤ S_nobrain ≤ S_relaxed (D46) | `python src/step5_check.py` | 0 violations |
| 5.6 | 5 | Step 5 results are deterministic | `python src/step5_check.py` | Two runs byte-identical |
| 6.1 | 6 | All 8 CPTAC cohorts with matched normals ran (D12) | `python src/step6_check.py` | 8 of 8 |
| 6.2 | 6 | No cohort below the sample minimum (D37) | `python src/step6_check.py` | 0 violations |
| 6.3 | 6 | Protein analysis ran unbiased, not only on candidates (D48) | `python src/step6_check.py` | 13,839 proteins tested vs 2,060 candidates |
| 6.4 | 6 | No global directional shift in protein log2FC (D51) | `python src/step6_check.py` | median ≈ 0; significant split near 50/50 |
| 6.5 | 6 | Concordance measured against the matched adjacent-normal arm (D49) | `python src/step6_check.py` | ρ ≥ 0.2; GTEx arm carried alongside |
| 6.6 | 6 | Nothing filtered on protein evidence (D48) | `python src/step6_check.py` | `rna_not_tested` present as a label, not a deletion |
| 7.1 | 7 | Eligibility rule honoured (D52) | `python src/step7_check.py` | No Flagged tier, all pass RNA, all ≥5 TPM |
| 7.2 | 7 | Ranked by S descending (D52) | `python src/step7_check.py` | Monotone |
| 7.3 | 7 | At most 2 genes per HGNC gene group (D52) | `python src/step7_check.py` | Max 1 observed |
| 7.4 | 7 | Selection rule recorded before the ranking existed (runbook rule 3) | `python src/step7_check.py` | D52 present in `docs/decisions.md`, dated 2026-10-21 before the first ranking run |
| 7.5 | 7 | Exported JSON matches the committed results | `python src/step7_check.py` | Same genes, same order |

## Pending — owned by later steps

| # | Step | Check |
| --- | --- | --- |
| 8.x | 8 | Every control accounted for; blood-cancer-only targets listed separately |
| 9.x | 9 | Default run reproduces Step 7 exactly; run count matches grid size; stability values in [0, 1] |

## Step 9 — robustness grid

| # | Check | Script |
| --- | --- | --- |
| 9.1 | Run count equals the grid size (270) | `src/step9_check.py` |
| 9.2 | Every setting appears exactly once | `src/step9_check.py` |
| 9.3 | The default cell reproduces the Step 7 top 10 exactly | `src/step9_check.py` |
| 9.4 | Stability values lie in [0, 1] | `src/step9_check.py` |
| 9.5 | Stability recomputes from `step9_runs.csv` | `src/step9_check.py` |
| 9.6 | Every run returned between 1 and 20 genes | `src/step9_check.py` |
| 9.7 | The 20 aa re-tier is byte-identical to the committed candidate set | `src/step9_robustness.py` |
| 9.8 | The gene-level tier collapse matches the `tier` column Step 5 stored | `src/step9_robustness.py` |

## Step 7b — the published deliverable

| # | Check | Script |
| --- | --- | --- |
| 7.10 | `__DATA__` matched exactly once and did not survive the build | `src/step7_build_artifact.py` |
| 7.11 | Withdrawn phrasing is absent; every tab has a section and is known to the switcher | `src/step7_build_artifact.py` |
| 7.12 | The page's JavaScript runs to completion and renders what it claims, including that all 270 selectable settings resolve to a stored run | `src/step7_check_artifact.py` |

## Extra — metric sensitivity (not a runbook step)

| # | Check | Script |
| --- | --- | --- |
| E.1 | Ranking count matches the specification (21) | `src/extra_metric_check.py` |
| E.2 | No duplicated (metric, reference) pair | `src/extra_metric_check.py` |
| E.3 | M1 S (both) reproduces the Step 7 top 10 exactly | `src/extra_metric_check.py` |
| E.4 | Stability values lie in [0, 1] | `src/extra_metric_check.py` |
| E.5 | Stability recomputes from `extra_metric_runs.csv` | `src/extra_metric_check.py` |
| E.6 | Every ranking returned between 1 and 20 genes | `src/extra_metric_check.py` |
| E.7 | The strictest prevalence gate actually changes each list | `src/extra_metric_check.py` |
