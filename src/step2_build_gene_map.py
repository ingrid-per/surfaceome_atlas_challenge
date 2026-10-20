"""Step 2: build the master identifier map every later step joins through.

The problem this solves: UniProt is keyed by accession, HPA and Xena by Ensembl
gene, CSPA and SURFY by accession. Join them on the wrong key and genes vanish
with no error at all.

Two rules, both from docs/decisions.md:

  D17  Join only on version-stripped Ensembl gene IDs. Symbols are carried so a
       human can read the table, and are never used as a key. A symbol join
       loses NECTIN4 from the Xena matrix, because GENCODE v23 still calls it
       PVRL4 — and NECTIN4 is the target of an approved drug, so that loss would
       have looked like the pipeline failing to find a known target.

  D15  One-to-many matches are kept and flagged, never dropped and never
       silently resolved.

Writes data/processed/gene_map.csv plus four reports under results/.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

csv.field_size_limit(10**9)

ENSG = re.compile(r"ENSG\d+")
report: list[str] = []


def out(line: str = "") -> None:
    print(line)
    report.append(line)


def strip_version(ensembl_id: str) -> str:
    """ENSG00000146648.22 -> ENSG00000146648.

    Every Ensembl ID entering any join goes through this one function, so the
    rule cannot be applied in one place and forgotten in another.
    """
    return ensembl_id.split(".")[0].strip()


# ------------------------------------------------------------------ inputs --
def load_spine() -> dict[str, dict]:
    """The 20,431 reviewed accessions: symbol, length, protein existence."""
    spine: dict[str, dict] = {}
    with gzip.open(RAW / "uniprot_human_reviewed_topology.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            spine[r["Entry"]] = {
                "symbol": r["Gene Names (primary)"] or "",
                "length": r["Length"],
                "protein_existence": r["Protein existence"],
            }
    out(f"UniProt reviewed accessions (spine)      : {len(spine):,}")
    return spine


def load_accession_to_genes() -> dict[str, set[str]]:
    """Accession -> Ensembl GENE ids, from the xref_ensembl_full field (D18)."""
    m: dict[str, set[str]] = {}
    with gzip.open(RAW / "uniprot_human_reviewed_xrefs.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = {strip_version(g) for g in ENSG.findall(r.get("Ensembl") or "")}
            if genes:
                m[r["Entry"]] = genes
    out(f"  accessions carrying an Ensembl gene    : {len(m):,}")
    return m


def load_hpa() -> dict[str, dict]:
    """Ensembl gene -> HPA symbol and the accessions HPA names for it."""
    path = RAW / "proteinatlas.tsv.zip"
    with zipfile.ZipFile(path) as z, z.open(z.namelist()[0]) as fh:
        rows = list(csv.DictReader(io.TextIOWrapper(fh, encoding="utf8"),
                                   delimiter="\t"))
    # Some HPA headers are quoted; match on the unquoted name.
    up_col = next(c for c in rows[0] if c.strip('"') == "Uniprot")
    hpa = {}
    for r in rows:
        g = strip_version(r["Ensembl"])
        if not g:
            continue
        hpa[g] = {
            "hpa_symbol": r["Gene"],
            "hpa_accessions": {a.strip() for a in (r[up_col] or "").split(",")
                               if a.strip()},
        }
    out(f"HPA gene table rows                      : {len(rows):,} "
        f"({len(hpa):,} distinct Ensembl genes)")
    return hpa


def load_xena() -> dict[str, str]:
    """Ensembl gene -> the symbol GENCODE v23 uses, which is often the old one."""
    xena = {}
    with (RAW / "gencode.v23.annotation.gene.probemap").open() as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            xena[strip_version(f[0])] = f[1]
    out(f"Xena probemap genes (GENCODE v23)        : {len(xena):,}")
    return xena


def load_cspa() -> set[str]:
    """Human accessions with experimental surface evidence. Table B is mouse."""
    df = pd.read_excel(RAW / "CSPA_S2_File.xlsx", sheet_name="Table A")
    human = df[df["organism"].astype(str).str.strip() == "Human"]
    accs = {str(a).strip() for a in human["ID_link"] if str(a).strip()}
    out(f"CSPA human surface proteins              : {len(accs):,}")
    return accs


def load_surfy() -> dict[str, str]:
    """Accession -> SURFY surfaceome label. The real header is the second row."""
    df = pd.read_excel(RAW / "SURFY_table_S3_surfaceome.xlsx",
                       sheet_name="SurfaceomeMasterTable", header=1)
    acc_col = next(c for c in df.columns if "accession" in str(c).lower())
    lab_col = next(c for c in df.columns if str(c).strip() == "Surfaceome Label")
    m = {str(a).strip(): str(l).strip()
         for a, l in zip(df[acc_col], df[lab_col]) if str(a).strip()}
    labels = pd.Series(list(m.values())).value_counts().to_dict()
    out(f"SURFY entries                            : {len(m):,}  labels={labels}")
    return m


# -------------------------------------------------------------------- build --
def main() -> int:
    PROC.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    out("STEP 2 MAPPING REPORT")
    out("Generated by src/step2_build_gene_map.py — do not edit by hand.")
    out("Join key: version-stripped Ensembl gene id. Symbols are never a key (D17).")
    out("")
    out("--- inputs " + "-" * 56)

    spine = load_spine()
    acc2genes = load_accession_to_genes()
    hpa = load_hpa()
    xena = load_xena()
    cspa = load_cspa()
    surfy = load_surfy()

    # Reverse direction, needed to flag genes claimed by several accessions.
    gene2accs: dict[str, set[str]] = defaultdict(set)
    for acc, genes in acc2genes.items():
        for g in genes:
            gene2accs[g].add(acc)

    rows, unmapped, multi, conflicts = [], [], [], []

    for acc in sorted(spine):
        info = spine[acc]
        genes = sorted(acc2genes.get(acc, ()))

        if not genes:
            unmapped.append({"accession": acc, "symbol": info["symbol"],
                             "protein_existence": info["protein_existence"],
                             "reason": "no Ensembl gene cross-reference in UniProt"})
            rows.append({
                "accession": acc, "symbol": info["symbol"], "ensembl_gene": "",
                "hpa_symbol": "", "xena_symbol": "",
                "length": info["length"],
                "protein_existence": info["protein_existence"],
                "n_genes_for_accession": 0, "n_accessions_for_gene": 0,
                "ambiguous_mapping": False, "no_ensembl": True,
                "in_hpa": False, "in_xena": False, "hpa_conflict": False,
                "cspa": acc in cspa,
                "surfy_label": surfy.get(acc, ""),
                "surfy_surface": surfy.get(acc, "") == "surface",
            })
            continue

        if len(genes) > 1:
            multi.append({"kind": "accession -> many genes", "accession": acc,
                          "symbol": info["symbol"], "n": len(genes),
                          "partners": ";".join(genes)})

        for g in genes:
            n_acc = len(gene2accs[g])
            if n_acc > 1 and g == genes[0]:
                multi.append({"kind": "gene -> many accessions", "accession": g,
                              "symbol": hpa.get(g, {}).get("hpa_symbol", ""),
                              "n": n_acc,
                              "partners": ";".join(sorted(gene2accs[g]))})

            h = hpa.get(g)
            conflict = False
            if h and h["hpa_accessions"] and acc not in h["hpa_accessions"]:
                # HPA names a different protein for this gene than UniProt does.
                if not (h["hpa_accessions"] & set(gene2accs[g])):
                    conflict = True
                    conflicts.append({
                        "ensembl_gene": g, "uniprot_accession": acc,
                        "uniprot_symbol": info["symbol"],
                        "hpa_symbol": h["hpa_symbol"],
                        "hpa_accessions": ";".join(sorted(h["hpa_accessions"]))})

            rows.append({
                "accession": acc, "symbol": info["symbol"], "ensembl_gene": g,
                "hpa_symbol": h["hpa_symbol"] if h else "",
                "xena_symbol": xena.get(g, ""),
                "length": info["length"],
                "protein_existence": info["protein_existence"],
                "n_genes_for_accession": len(genes),
                "n_accessions_for_gene": n_acc,
                "ambiguous_mapping": len(genes) > 1 or n_acc > 1,
                "no_ensembl": False,
                "in_hpa": h is not None,
                "in_xena": g in xena,
                "hpa_conflict": conflict,
                "cspa": acc in cspa,
                "surfy_label": surfy.get(acc, ""),
                "surfy_surface": surfy.get(acc, "") == "surface",
            })

    df = pd.DataFrame(rows).sort_values(["accession", "ensembl_gene"],
                                        kind="mergesort")

    # ---------------------------------------------------------- the funnel --
    n_acc = df["accession"].nunique()
    mapped = df[~df["no_ensembl"]]
    out("")
    out("--- coverage " + "-" * 54)
    out(f"accessions in                            : {n_acc:,}")
    out(f"  with an Ensembl gene                   : {mapped['accession'].nunique():,}"
        f"  ({mapped['accession'].nunique() / n_acc:.1%})")
    out(f"  without one (listed in step2_unmapped) : {len(unmapped):,}"
        f"  ({len(unmapped) / n_acc:.1%})")
    out("")
    out(f"distinct Ensembl genes reached           : {mapped['ensembl_gene'].nunique():,}")
    out(f"  found in HPA                           : {mapped[mapped.in_hpa]['ensembl_gene'].nunique():,}")
    out(f"  found in Xena (GENCODE v23)            : {mapped[mapped.in_xena]['ensembl_gene'].nunique():,}")
    out(f"  in both                                : "
        f"{mapped[mapped.in_hpa & mapped.in_xena]['ensembl_gene'].nunique():,}")
    out(f"  in neither                             : "
        f"{mapped[~mapped.in_hpa & ~mapped.in_xena]['ensembl_gene'].nunique():,}")
    out("")
    out(f"rows flagged ambiguous_mapping           : {int(df.ambiguous_mapping.sum()):,}")
    out(f"HPA / UniProt accession conflicts        : {len(conflicts):,}")
    out("")
    out(f"accessions with CSPA surface evidence    : {int(df.groupby('accession').cspa.any().sum()):,}")
    out(f"accessions SURFY calls 'surface'         : {int(df.groupby('accession').surfy_surface.any().sum()):,}")

    # Genes whose current symbol differs from the one GENCODE v23 uses. This is
    # the exact population a symbol-based join would mishandle.
    renamed = mapped[(mapped.xena_symbol != "") & (mapped.hpa_symbol != "")
                     & (mapped.xena_symbol != mapped.hpa_symbol)]
    out("")
    out(f"genes renamed since GENCODE v23          : "
        f"{renamed['ensembl_gene'].nunique():,}")
    out("  (these are exactly the genes a symbol join would lose; D17)")
    for _, r in renamed.head(8).iterrows():
        out(f"    {r.ensembl_gene}  {r.hpa_symbol:<12} was {r.xena_symbol}")

    out("")
    out(f"rows written                             : {len(df):,}")

    # ------------------------------------------------------------- outputs --
    df.to_csv(PROC / "gene_map.csv", index=False)
    pd.DataFrame(unmapped).to_csv(RESULTS / "step2_unmapped.csv", index=False)
    pd.DataFrame(multi).to_csv(RESULTS / "step2_multimapped.csv", index=False)
    pd.DataFrame(conflicts).to_csv(RESULTS / "step2_conflicts.csv", index=False)
    (RESULTS / "step2_mapping_report.txt").write_text("\n".join(report) + "\n")

    print(f"\nWrote data/processed/gene_map.csv ({len(df):,} rows) "
          "and 4 reports under results/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
