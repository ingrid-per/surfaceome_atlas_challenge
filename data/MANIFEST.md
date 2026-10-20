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

| File | Resource | Version / release | Date retrieved | URL | SHA-256 | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `uniprot_human_reviewed.fasta.gz` | UniProt | 2026_03 (02-September-2026) | 2026-10-20 | <https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29&format=fasta&compressed=true> | `dc8ef8cbf0fc1b857b3be9a559613533ec49824079fb057e31fdfcf578cb229a` | 7.16 MB |
| `uniprot_human_reviewed_topology.tsv.gz` | UniProt | 2026_03 (02-September-2026) | 2026-10-20 | <https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29&format=tsv&fields=accession,id,protein_name,gene_primary,gene_names,xref_ensembl,length,reviewed,protein_existence,ft_transmem,ft_topo_dom,ft_signal,ft_lipid,ft_intramem,cc_subcellular_location,ft_carbohyd,sequence&compressed=true> | `ef44fded2ba2bdf87a65e63ffb3c784744cd6e333e4b925ed34e9f5a5e334f7f` | 10.13 MB |
| `uniprot_human_reviewed_xrefs.tsv.gz` | UniProt | 2026_03 (02-September-2026) | 2026-10-20 | <https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29&format=tsv&fields=accession,gene_primary,gene_names,xref_ensembl_full,protein_existence&compressed=true> | `e52c01b3bc587a6cf394a821402dabd2d966277223e256c606048f010ef4c1e2` | 1.42 MB |
| `proteinatlas.tsv.zip` | HPA | 25.1 (Ensembl 109) | 2026-10-20 | <https://www.proteinatlas.org/download/proteinatlas.tsv.zip> | `ad401f8519ecdee5b67b0ad3bf28175d402d1ea0ad6c542b90fcab2bee39b470` | 7.11 MB |
| `rna_tissue_hpa.tsv.zip` | HPA | 25.1 (Ensembl 109) | 2026-10-20 | <https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip> | `735d0f5de5374946914e74949a760457f93b9ea4ca6200e0be55f22bb0473565` | 6.35 MB |
| `normal_ihc_data.tsv.zip` | HPA | 25.1 (Ensembl 109) | 2026-10-20 | <https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip> | `fa2315e5a53ecffd9427a5bfe219c8cc41b54089ef6f6be7b5d79b5455e0677f` | 5.47 MB |
| `subcellular_location.tsv.zip` | HPA | 25.1 (Ensembl 109) | 2026-10-20 | <https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip> | `a9466882bd3328416bcde77c32a6bcb2a17d70183af90e98d032f6f5f781d504` | 246 KB |
| `TcgaTargetGtex_rsem_gene_tpm.gz` | Xena Toil | Toil hub, files dated 2021-04-09, GENCODE v23 | 2026-10-20 | <https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz> | `a8c36cb16ef82eccb0b6e1e62a2af3ab9d54c1f56735bc59dfac4cbbe2391d4d` | 1261.95 MB |
| `TcgaTargetGTEX_phenotype.txt.gz` | Xena Toil | Toil hub, files dated 2021-04-09, GENCODE v23 | 2026-10-20 | <https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz> | `ba4d4461cff0fe5e91cc3e58f793aa47ff3c5d6fbf77d65ccff5a07231aaf2db` | 133 KB |
| `gencode.v23.annotation.gene.probemap` | Xena Toil | Toil hub, files dated 2021-04-09, GENCODE v23 | 2026-10-20 | <https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap> | `6783ea58791ae876efb697889a042cc7be8e32e40fc01191c622ef25d9416931` | 3.09 MB |
| `BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `e584331bc77a440af50692c7874b736b3d135d01d96ce95d66af18cc64009b21` | 20.41 MB |
| `BRCA_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_meta.txt> | `16d00a3022fbafa9d29e74c0ee689ae63279911f04aeef2d21344ed51a0a5cb7` | 7 KB |
| `CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `090077c25f0eeb59b737912f39dec98901db949e1f143bbc481d5cd731550f60` | 15.86 MB |
| `CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `1b00d684cd99d870bee199128485c22d6cb079e0a01de3c044b7a4fdb1735618` | 12.42 MB |
| `CCRCC_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_meta.txt> | `c332231f3cbd204094236dfd302977e70bbdddcf0049e862b60dd10244d04dae` | 11 KB |
| `COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `fa2dde5f9165acdf6320d28e01dbad24b735ede320f5fab4a6ceab6c15aa0580` | 11.72 MB |
| `COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `a76a0cc0ab3986635e8745028eb0bcc64470035ab737cc0937ec7be114a09177` | 11.94 MB |
| `COAD_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_meta.txt> | `5242103f57fb5d919ab7e228dd36c97b669915d6bee88ae1e5a43c3aa270295e` | 11 KB |
| `GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `d7393915eb8ec10c1c48bea03d78cd9ee2a3477f7bb4b93f808ffade03c2f01b` | 17.20 MB |
| `GBM_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_meta.txt> | `6617ef2f9010acec1eba09a765ca3e35e4ba2ea7ba75ba9e8734ed9ad94feecc` | 7 KB |
| `HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `6d52dac50c55a4aa30e0aabaf3047a3a099875ea075b65d8f0bc83ac45a5c4b3` | 17.34 MB |
| `HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `0a067886f9b553d3de24d2dbab21e18037c44a50a52b9d60004d22dba3442adc` | 10.05 MB |
| `HNSCC_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_meta.txt> | `3c9e32ff8b51b4a22f967ca753185e01d79f69d860f9530f5587314bd06c80ac` | 12 KB |
| `LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `fd9f9cc0aed40942fd33d750d1d207bcdf9ea231fcddeff3cc8ac8fdf3e1384d` | 19.68 MB |
| `LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `ea3f8e5a1e42f708af87f6cdfc86d3fa63a96bfdc4f6f3316ba54380ea579aaf` | 18.04 MB |
| `LSCC_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_meta.txt> | `6c334bbd8c9c4acfc94709a7fea5fd0181efa0881cee8ce1793085a5281a9d99` | 12 KB |
| `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `746900f54072198481207ea13f8a7ba14dcabd52345bace0bef7438120c03845` | 19.20 MB |
| `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `3cfc62b159276837dea644f457ebf6b212964320e0b3d8ac87168d82eb73747c` | 17.64 MB |
| `LUAD_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_meta.txt> | `f9facf713af130ad097c56794f1fbc9107334e02ad9bf0df3958b40d68b6af10` | 12 KB |
| `OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `274c22512ab70e2a07197ca6e51d268038c90e903bac823c1b27b7950af4ae6f` | 11.93 MB |
| `OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `6fc249b27c9284b598a94860c1e36f4fd45843f65b14c4ccdfbc0153725bb07a` | 2.87 MB |
| `OV_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_meta.txt> | `ab6b3bb70460c682492a98f7620cbbeda61baab6e4d8bc1e25d19c0580390123` | 6 KB |
| `PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `e1be4bed67bfdac6fff7b12487b8af9e0e2a0e93901039f5f9124e180f24f6bd` | 16.02 MB |
| `PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `beb1bc33ee52d10f0f7a6fa0bfc0de97916471c37d80c5cb655ee0d0a2452ce9` | 6.84 MB |
| `PDAC_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_meta.txt> | `ce53f0acf92cbbbf165c8c27425c27b7b6de301c2e95279eaf124c37a2316d33` | 11 KB |
| `UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt> | `aea5d956b3496530faf86052886de19e857e96294a2dbcb3adb9702457865d36` | 15.79 MB |
| `UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt> | `225df3e7827fbeec499deba3d773828935dbe722559730a902c66f9dabe1d229` | 3.20 MB |
| `UCEC_meta.txt` | CPTAC | LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized | 2026-10-20 | <https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_meta.txt> | `4058a8420129dd37f4d9035e1d87beb70a800ce9cb42754f74e7ad3d35f3020b` | 14 KB |
| `Jiang2020_mmc3.xlsx` | GTEx proteome | Jiang et al., Cell 2020, 183(1):269-283 | 2026-10-20 | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx> | `f278d406990059a625cb0d5c4c29e56dd8634052138a38bd053bc0d61124af80` | 117.28 MB |
| `CSPA_S2_File.xlsx` | CSPA | CSPA, 2015 | 2026-10-20 | <https://wollscheidlab.org/CSPA/data/S2_File.xlsx> | `9b02c85d4c1d49847cf37e23755f07e57ea65b79cc0ebcdaeef1752dab8b5c8a` | 1.36 MB |
| `SURFY_table_S3_surfaceome.xlsx` | SURFY | SURFY, 2018 | 2026-10-20 | <https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx> | `2f1b8262463ce1c59a1f945d22f0e9638cb3bfbf5aabe197f43b562a62fb055a` | 6.55 MB |

**41 files, 1676.50 MB total.**

### Deliberately not downloaded

| Resource | Why not |
| --- | --- |
| HPA pathology.tsv.zip (cancer IHC) | D8: not published in HPA 25.1; only on the archived v23 site |
| HPA rna_tissue_gtex.tsv.zip | D9: GTEx rescaled into HPA nTPM would make both arms of the Step 5 safety check HPA-processed, destroying their independence |
| HPA rna_tissue_consensus.tsv.zip | D9: blends HPA and GTEx into one value, which is the opposite of what D2 asks for |
| CPTAC RNA-seq TPM files | The TCGA/GTEx Toil matrix already covers the RNA layer |
| GTEx portal median-TPM file | D9: GTEx medians are computed from the Toil matrix instead, keeping GTEx on the same pipeline as the TCGA tumors |

## Expected resources (Step 1 checklist)

Six resources, from the runbook's Step 1 table. All resolved by probing the live
endpoints on 2026-10-20, reviewed before downloading, and now present with checksums
in the Files table above.

| # | Resource | What it gives | Status |
| --- | --- | --- | --- |
| R1 | UniProt, reviewed human proteome | Sequences, curated topology, signal peptides, lipid anchors | Downloaded 2026-10-20 |
| R2 | Human Protein Atlas downloads | Normal-tissue RNA (~50 tissues), antibody staining (IHC), gene ID table | Downloaded 2026-10-20 |
| R3 | UCSC Xena Toil hub (TCGA + GTEx) | Tumor, matched normal and GTEx normal RNA in one pipeline | Downloaded 2026-10-20 |
| R4 | CPTAC, via LinkedOmicsKB or the Proteomic Data Commons | Tumor and adjacent-normal protein levels, 10 cancer types | Downloaded 2026-10-20 |
| R5 | GTEx proteome (Jiang et al., Cell 2020) | Protein levels, 12,627 genes across 32 normal tissues | Downloaded 2026-10-20 |
| R6 | Cell Surface Protein Atlas (CSPA) | 1,492 human proteins with experimental surface evidence | Downloaded 2026-10-20 |

Also needed in Step 3, not a Step 1 download: the SURFY predicted surfaceome
(2,886 proteins).

## Access notes

- `rest.uniprot.org` returns **HTTP 403** to requests that send no `User-Agent`
  header, and 200 with one. Verified 2026-10-20. A download script that omits the
  header will appear to fail for the wrong reason.
- Reviewed human entries in UniProt at 2026-10-20: **20,431** (`x-total-results`
  from the REST search endpoint). Use as a sanity figure for R1.
