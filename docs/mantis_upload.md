# Uploading to Mantis

Two files, both written by `src/extra_mantis_export.py`.

| file | rows | size | one row per |
| --- | --- | --- | --- |
| `results/mantis_genes.csv` | 1,272 | 0.85 MB | eligible gene, through its best cancer type |
| `results/mantis_per_cancer.csv` | 296 | 0.20 MB | top 10 by S within each cancer type |
| `results/mantis_candidates.csv` | 3,442 | 2.31 MB | eligible gene × cancer type |

**The upload endpoint refuses bodies over roughly 1 MB** (nginx default), so
`mantis_candidates.csv` will fail with a 413. Use one of the other two.

Start with `mantis_genes.csv` — one point per gene is the cleaner first map.
`mantis_per_cancer.csv` is the better second one: the best candidates in every
indication, without the long tail of genes that merely cleared the bar.

## Column typing at upload

Mantis blocks creation until at least one column is marked **semantic**.

| column | type | note |
| --- | --- | --- |
| `label` | **title** | "CLDN6 — Testicular Germ Cell Tumor" |
| `description` | **semantic** | the embedded text; median 554 characters |
| `hpa_link`, `uniprot_link` | **links** | |
| `coordinate1_tumour_log2tpm` | numeric, or **coordinate1** | |
| `coordinate2_healthy_log2tpm` | numeric, or **coordinate2** | |
| `S`, `log2fc`, `tumour_tpm`, `healthy_max_tpm`, `frac_tumours_above_10tpm`, `threshold_stability`, `gtex_tau`, `healthy_breadth_above_10tpm`, `tumour_tpm_q25`, `S_relaxed`, `fdr` | **numeric** | colour and filter by these |
| `selection`, `cancer_type`, `tier`, `comparator_grade`, `compartment`, `protein_existence`, `agreement`, `trajectory`, `gene_group`, `healthy_max_tissue`, `healthy_max_source` | **categoric** | |
| `low_power`, `gpi`, `cspa`, `surfy` | **categoric** | booleans |
| `symbol`, `ensembl_gene`, `accession`, `extracellular_reported` | unused, or categoric | keys |

## Two different maps from the same file

**Embedding layout** — leave the coordinates as plain numerics. Mantis embeds
`description` and places genes by functional similarity. Colour by `S` or by
`selection` to see whether the top 10 clusters in one region of the surfaceome
or scatters across it.

**Explicit scatter** — type the two coordinate columns as `coordinate1` and
`coordinate2`. Axes become tumour expression and worst healthy expression,
which is the figure the runbook asks for at Step 11, interactive rather than
static. Points near the diagonal have no therapeutic window; the useful region
is bottom-right.

## What the description is built from

Four fields: UniProt protein name, HGNC gene family, HPA molecular function
and biological process, and UniProt subcellular location with evidence codes
stripped. Capped at 400 characters; median 296.

Three things are deliberately left out. HPA's "Protein class" is near-identical
across rows ("Predicted membrane proteins"), so it adds length without
separating anything. The cancer type and the highest healthy tissue are already
categoric columns — repeating them in the embedded text would make the map
cluster by cancer type rather than by protein biology, which is the one thing
this map can show that the existing figures cannot.

Gene symbols are deliberately **not** the semantic column. Embedding `CLDN6`
would place it beside `CLDN3` because the strings look alike, which maps gene
naming rather than biology.
