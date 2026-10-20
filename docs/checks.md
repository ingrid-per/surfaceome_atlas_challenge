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

## Pending — owned by later steps

| # | Step | Check |
| --- | --- | --- |
| 3.x | 3 | Candidate count roughly 1,500–4,500; EGFR/ERBB2 in S1, MSLN in S3; GAPDH/ACTB absent; tier counts printed |
| 4.x | 4 | P-value histograms sensible per cancer type; counts and group sizes printed |
| 5.x | 5 | B2M scores low; 10 genes hand-verified against the tissue table |
| 6.x | 6 | Funnel counts monotone; 5 tier-1 genes hand-verified; unmeasured genes counted as RNA-only |
| 7.x | 7 | Output deterministic across two runs; three near-miss genes explained; rule dated before the first ranking file |
| 8.x | 8 | Every control accounted for; blood-cancer-only targets listed separately |
| 9.x | 9 | Default run reproduces Step 7 exactly; run count matches grid size; stability values in [0, 1] |
