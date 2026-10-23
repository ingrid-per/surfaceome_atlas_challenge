"""The exact files Step 1 downloads, pinned (decision D7, 2026-10-20).

Every URL here was probed live before being written down — none is from memory.
Sizes are the bytes the server reported at probe time and are used only to warn
about a mismatch, never to decide correctness; the checksum does that.

Two access traps are encoded here deliberately:

1. rest.uniprot.org returns 403 to any request with no User-Agent header.
2. The old CSPA/SURFY host wlab.ethz.ch now serves an HTML page with HTTP 200 at
   the former .xlsx paths. A downloader that trusts the status code would write an
   80 KB web page to disk named .xlsx and only fail several steps later. Hence
   `magic`: every file is checked for its format's leading bytes before it counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field

USER_AGENT = "surface-prot-challenge/0.1 (academic target-discovery pipeline)"

# Leading bytes each container format must start with.
MAGIC = {
    "gzip": b"\x1f\x8b",
    "zip": b"PK\x03\x04",   # also covers .xlsx, which is a zip
    "text": None,           # plain text: checked for a non-HTML first line instead
}


@dataclass(frozen=True)
class Resource:
    rid: str
    group: str
    filename: str
    url: str
    expected_bytes: int | None
    magic: str
    note: str
    release: str = ""


# --- UniProt -----------------------------------------------------------------
# Release 2026_03, 02-September-2026. 20,431 reviewed human entries.
_UP = "https://rest.uniprot.org/uniprotkb/stream"
_UP_QUERY = "query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29"
_UP_FIELDS = (
    "accession,id,protein_name,gene_primary,gene_names,xref_ensembl,length,"
    "reviewed,protein_existence,ft_transmem,ft_topo_dom,ft_signal,ft_lipid,"
    "ft_intramem,cc_subcellular_location,ft_carbohyd,sequence"
)
UNIPROT_RELEASE = "2026_03 (02-September-2026)"

# --- hosts -------------------------------------------------------------------
_HPA = "https://www.proteinatlas.org/download"
_XENA = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download"
_CPTAC = ("https://cptac-pancancer-data.s3.us-west-2.amazonaws.com"
          "/data_freeze_v1.2_reorganized")
_ELS = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783"
_WOLL = "https://wollscheidlab.org"

# CPTAC cohorts. Normal proteomics exists for 8 of 10 — not BRCA, not GBM (D12).
CPTAC_COHORTS = ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
                 "LSCC", "LUAD", "OV", "PDAC", "UCEC"]
CPTAC_NO_NORMAL = {"BRCA", "GBM"}
_CPTAC_PROT = "proteomics_gene_abundance_log2_reference_intensity_normalized"


def _build() -> list[Resource]:
    r: list[Resource] = []

    r += [
        Resource("U1", "UniProt", "uniprot_human_reviewed.fasta.gz",
                 f"{_UP}?{_UP_QUERY}&format=fasta&compressed=true",
                 7_508_673, "gzip",
                 "Reviewed human sequences; input to DeepTMHMM in Step 3",
                 UNIPROT_RELEASE),
        Resource("U2", "UniProt", "uniprot_human_reviewed_topology.tsv.gz",
                 f"{_UP}?{_UP_QUERY}&format=tsv&fields={_UP_FIELDS}&compressed=true",
                 10_624_831, "gzip",
                 "Curated topology: transmembrane, topological domain, signal "
                 "peptide, lipid anchor, subcellular location",
                 UNIPROT_RELEASE),
        # U2's `xref_ensembl` turned out to return transcript IDs (ENST), not
        # gene IDs, so it carries no gene-level link at all. `xref_ensembl_full`
        # does. U2 stays — its topology columns are what Step 3 needs. (D18)
        Resource("U3", "UniProt", "uniprot_human_reviewed_xrefs.tsv.gz",
                 f"{_UP}?{_UP_QUERY}&format=tsv"
                 f"&fields=accession,gene_primary,gene_names,xref_ensembl_full,"
                 f"protein_existence&compressed=true",
                 1_487_846, "gzip",
                 "Accession -> Ensembl GENE id (ENSG) — the join key for Step 2",
                 UNIPROT_RELEASE),
    ]

    # HPA 25.1 (Ensembl 109). Only proteinatlas.tsv.zip is documented on the
    # download page; the /download/tsv/ paths are live but undocumented, so they
    # are pinned here rather than discovered at run time.
    hpa = "25.1 (Ensembl 109)"
    r += [
        Resource("H1", "HPA", "proteinatlas.tsv.zip",
                 f"{_HPA}/proteinatlas.tsv.zip", 7_456_000, "zip",
                 "Gene table: Ensembl, symbol, UniProt — the spine of Step 2", hpa),
        Resource("H2", "HPA", "rna_tissue_hpa.tsv.zip",
                 f"{_HPA}/tsv/rna_tissue_hpa.tsv.zip", 6_658_000, "zip",
                 "Normal-tissue nTPM across HPA's own tissues — HPA arm of Step 5", hpa),
        Resource("H3", "HPA", "normal_ihc_data.tsv.zip",
                 f"{_HPA}/tsv/normal_ihc_data.tsv.zip", 5_736_000, "zip",
                 "Antibody staining in 45 normal tissues with reliability — Step 6", hpa),
        Resource("H4", "HPA", "subcellular_location.tsv.zip",
                 f"{_HPA}/tsv/subcellular_location.tsv.zip", 251_000, "zip",
                 "HPA subcellular location calls — third opinion in Step 3", hpa),
    ]

    xena = "Toil hub, files dated 2021-04-09, GENCODE v23"
    r += [
        Resource("X1", "Xena Toil", "TcgaTargetGtex_rsem_gene_tpm.gz",
                 f"{_XENA}/TcgaTargetGtex_rsem_gene_tpm.gz",
                 1_323_366_000, "gzip",
                 "60,498 genes x 19,131 samples in log2(TPM+0.001) — see D11", xena),
        Resource("X2", "Xena Toil", "TcgaTargetGTEX_phenotype.txt.gz",
                 f"{_XENA}/TcgaTargetGTEX_phenotype.txt.gz", 136_000, "gzip",
                 "Per-sample study, sample type and primary site — drives D10", xena),
        Resource("X3", "Xena Toil", "gencode.v23.annotation.gene.probemap",
                 f"{_XENA}/probeMap%2Fgencode.v23.annotation.gene.probemap",
                 3_250_000, "text",
                 "Ensembl gene ID to symbol for the matrix rows", xena),
    ]

    cptac = "LinkedOmicsKB Release 1 (2023-05-08), freeze v1.2_reorganized"
    for c in CPTAC_COHORTS:
        kinds = ["Tumor"] if c in CPTAC_NO_NORMAL else ["Tumor", "Normal"]
        for k in kinds:
            r.append(Resource(
                f"C_{c}_{k}", "CPTAC", f"{c}_{_CPTAC_PROT}_{k}.txt",
                f"{_CPTAC}/{c}/{c}_{_CPTAC_PROT}_{k}.txt", None, "text",
                f"{c} {k.lower()} protein abundance, log2 reference-intensity "
                f"normalized", cptac))
        r.append(Resource(
            f"C_{c}_meta", "CPTAC", f"{c}_meta.txt",
            f"{_CPTAC}/{c}/{c}_meta.txt", None, "text",
            f"{c} sample metadata", cptac))

    r.append(Resource(
        "H5", "HGNC", "hgnc_complete_set.txt",
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
        "hgnc_complete_set.txt",
        17_000_000, "text",
        "Gene groups (families) for the D52 per-family limit",
        "HGNC current release"))

    r += [
        Resource("G1", "GTEx proteome", "Jiang2020_mmc3.xlsx",
                 f"{_ELS}-mmc3.xlsx", 122_974_000, "zip",
                 "Cell 2020 Table S3: protein abundance over 32 tissues. Sheets "
                 "'D protein relative abundance' and 'E protein tissue median'",
                 "Jiang et al., Cell 2020, 183(1):269-283"),
        Resource("S1", "CSPA", "CSPA_S2_File.xlsx",
                 f"{_WOLL}/CSPA/data/S2_File.xlsx", 1_430_000, "zip",
                 "1,492 validated human cell-surface proteins", "CSPA, 2015"),
        Resource("S2", "SURFY", "SURFY_table_S3_surfaceome.xlsx",
                 f"{_WOLL}/SURFY/table_S3_surfaceome.xlsx", 6_870_000, "zip",
                 "2,886-protein predicted surfaceome", "SURFY, 2018"),
        # Downloaded for the tumour-purity analysis that D61 withdrew. The
        # ANALYSIS is gone; the FILE is still in data/raw, so it still needs a
        # manifest row. It is declared here rather than hand-written into the
        # manifest because write_manifest() rebuilds that file from this list,
        # and anything not declared here gets silently dropped the next time a
        # download runs — which is exactly what happened to it.
        Resource("A1", "Tumour purity (withdrawn, see D61)",
                 "aran2015_tumour_purity.xlsx",
                 "https://static-content.springer.com/esm/art%3A10.1038%2F"
                 "ncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx",
                 553_082, "zip",
                 "Aran et al. Supplementary Data 1. Kept on disk and recorded; "
                 "no step reads it (D61)",
                 "Aran et al., Nat Commun 6:8971 (2015), Supplementary Data 1"),
        # Step 12. The documented host bio-bigdata.hrbmu.edu.cn answers 301 and
        # then 403; this mirror serves the same file and is what was actually
        # fetched, so it is what the manifest records.
        Resource("M1", "CellMarker 2.0", "CellMarker2_Cell_marker_Human.xlsx",
                 "http://117.50.127.228/CellMarker/CellMarker_download_files"
                 "/file/Cell_marker_Human.xlsx", 7_982_475, "zip",
                 "Curated human cell-type marker genes, used in Step 12 to "
                 "build the lineage composite",
                 "CellMarker 2.0 (mirror host)"),
    ]
    return r


RESOURCES: list[Resource] = _build()

# Deliberately NOT downloaded, with the reason. Recorded so the omission is a
# decision rather than something nobody thought of.
EXCLUDED = [
    ("HPA pathology.tsv.zip (cancer IHC)",
     "D8: not published in HPA 25.1; only on the archived v23 site"),
    ("HPA rna_tissue_gtex.tsv.zip",
     "D9: GTEx rescaled into HPA nTPM would make both arms of the Step 5 "
     "safety check HPA-processed, destroying their independence"),
    ("HPA rna_tissue_consensus.tsv.zip",
     "D9: blends HPA and GTEx into one value, which is the opposite of what "
     "D2 asks for"),
    ("CPTAC RNA-seq TPM files",
     "The TCGA/GTEx Toil matrix already covers the RNA layer"),
    ("GTEx portal median-TPM file",
     "D9: GTEx medians are computed from the Toil matrix instead, keeping GTEx "
     "on the same pipeline as the TCGA tumors"),
]
