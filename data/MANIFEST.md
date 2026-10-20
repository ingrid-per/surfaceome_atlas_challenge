# Data manifest

**No file enters `data/raw/` without a row in this table.** HPA and TCGA-derived
resources change between releases, and the same code gives different results on
different versions. A result is only reproducible if the exact release is recorded.

`data/raw/` is read-only. Never delete or edit anything in it.

Checksums are SHA-256. Generate with:

```bash
shasum -a 256 data/raw/<file>
```

## Files

| File | Resource | Version / release | Date retrieved | URL | SHA-256 | Size | Rows |
| --- | --- | --- | --- | --- | --- | --- | --- |
| _none yet_ | | | | | | | |

## Expected resources (Step 1 checklist)

Six resources, from the runbook's Step 1 table. Exact file names, URLs and sizes are
**unverified** — Claude Code must find them on the official pages and show them to
I before anything is downloaded.

| # | Resource | What it gives | Status |
| --- | --- | --- | --- |
| R1 | UniProt, reviewed human proteome | Sequences, curated topology, signal peptides, lipid anchors | NOT YET DOWNLOADED |
| R2 | Human Protein Atlas downloads | Normal-tissue RNA (~50 tissues), antibody staining (IHC), gene ID table | NOT YET DOWNLOADED |
| R3 | UCSC Xena Toil hub (TCGA + GTEx) | Tumor, matched normal and GTEx normal RNA in one pipeline | NOT YET DOWNLOADED |
| R4 | CPTAC, via LinkedOmicsKB or the Proteomic Data Commons | Tumor and adjacent-normal protein levels, 10 cancer types | NOT YET DOWNLOADED |
| R5 | GTEx proteome (Jiang et al., Cell 2020) | Protein levels, 12,627 genes across 32 normal tissues | NOT YET DOWNLOADED |
| R6 | Cell Surface Protein Atlas (CSPA) | 1,492 human proteins with experimental surface evidence | NOT YET DOWNLOADED |

Also needed in Step 3, not a Step 1 download: the SURFY predicted surfaceome
(2,886 proteins).

## Access notes

- `rest.uniprot.org` returns **HTTP 403** to requests that send no `User-Agent`
  header, and 200 with one. Verified 2026-10-20. A download script that omits the
  header will appear to fail for the wrong reason.
- Reviewed human entries in UniProt at 2026-10-20: **20,431** (`x-total-results`
  from the REST search endpoint). Use as a sanity figure for R1.
