"""Step 3c: assign every reviewed protein a surface tier.

Two independent sources decide: UniProt's curated topology (step3a) and
DeepTMHMM's prediction (step3b). Agreement is stronger evidence than either
alone, which is what the tiers encode.

    S1       transmembrane, longest extracellular segment >= MIN_EXTRA,
             called by BOTH sources                                  (high)
    S2       same, but called by only ONE of the two               (medium)
    S3       GPI-anchored, by curated annotation or NetGPI          (medium)
    Flagged  surface call but longest extracellular < MIN_EXTRA        (low)
    -        not a surface candidate

Decisions this implements:

  D20  Default 20 residues, but the measured longest single extracellular
       segment is stored as a number for every protein, so Step 9 re-runs the
       10/20/50 grid from one column instead of re-running the prediction.
       GPI proteins get `not applicable`, never 0 — they sit entirely outside
       the membrane, so the measure is meaningless for them, and a 0 would
       make them fail a length filter they should never face.
  D22  GPI calls record their source (curated / predicted / both), and the
       report states how many rest on prediction alone.
  D30  No SignalP. DeepTMHMM's own signal-peptide calls are used.

NetGPI output is optional: without it, GPI evidence is curated-only and every
affected protein is flagged so the gap is visible rather than assumed away.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_EXTRA_DEFAULT = 20   # D20; a parameter, never hard-coded into the tiers
# Must not be a string pandas reads back as NaN ("NA", "N/A", "null" all are),
# or the distinction D20 asks for is destroyed on the CSV round-trip.
NOT_APPLICABLE = "not_applicable"


# D31. A topology predictor says "outside" relative to whichever membrane the
# protein sits in, not relative to the cell. For a mitochondrial protein
# "outside" is the intermembrane space, which no drug reaches. These patterns
# identify compartment claims in UniProt's subcellular-location text.
PLASMA_PAT = r"cell membrane|plasma membrane|cell surface"
INTERNAL_PAT = (r"mitochondri|endoplasmic reticulum|golgi|lysosom|peroxisom|"
                r"nuclear membrane|endosom|sarcoplasmic|melanosom|vacuol|"
                r"autophagosom|phagosom")


def hpa_plasma_membrane_genes() -> set[str]:
    """Ensembl genes HPA places on the plasma membrane — the independent second
    opinion on compartment (resource H4, downloaded in Step 1 for this)."""
    import io
    import zipfile
    path = RAW / "subcellular_location.tsv.zip"
    if not path.exists():
        return set()
    with zipfile.ZipFile(path) as z, z.open(z.namelist()[0]) as fh:
        hpa = pd.read_csv(io.TextIOWrapper(fh, encoding="utf8"), sep="\t")
    cols = [c for c in hpa.columns if "location" in c.lower()]
    text = hpa[cols].fillna("").agg(";".join, axis=1).str.lower()
    return set(hpa.loc[text.str.contains("plasma membrane"), "Gene"])


def load_netgpi() -> tuple[dict[str, bool], str]:
    """NetGPI output, if I has submitted it yet (D25).

    Accepts the server's tabular output; any file named *netgpi* in data/raw.
    Returns (accession -> is_gpi, provenance note).
    """
    hits = sorted(RAW.glob("*netgpi*"))
    if not hits:
        return {}, "not run"
    path = hits[0]
    calls: dict[str, bool] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t") if "\t" in line else line.split()
        if len(f) < 2:
            continue
        acc = f[0].strip()
        verdict = " ".join(f[1:]).lower()
        if "gpi-anchored" in verdict or verdict.startswith("gpi"):
            calls[acc] = True
        elif "not gpi" in verdict or "non-gpi" in verdict:
            calls[acc] = False
    return calls, path.name


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictor", default="full", choices=["full", "pilot"],
                    help="which DeepTMHMM run to use")
    ap.add_argument("--min-extracellular", type=int, default=MIN_EXTRA_DEFAULT,
                    help="D20 default is 20; Step 9 varies this")
    args = ap.parse_args()

    uni = pd.read_csv(PROC / "step3_uniprot_topology.csv")
    dtm_path = PROC / f"step3_deeptmhmm_{args.predictor}.csv"
    if not dtm_path.exists():
        print(f"{dtm_path} not found — run src/step3_run_deeptmhmm.py first")
        return 1
    dtm = pd.read_csv(dtm_path)
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)

    df = uni.merge(dtm, on="accession", how="left")

    # D21: a protein with no predictor result can never reach S1, which needs
    # two sources. Flag rather than quietly treat absence as a negative call.
    df["predictor_missing"] = df["dtm_type"].isna()
    for col, fill in [("dtm_surface_call", False), ("dtm_n_tm", 0),
                      ("dtm_longest_outside", 0), ("dtm_has_signal", False)]:
        df[col] = df[col].fillna(fill)
    df["dtm_surface_call"] = df["dtm_surface_call"].astype(bool)

    # ---- GPI evidence and its source (D22) ----
    netgpi, netgpi_src = load_netgpi()
    df["gpi_predicted"] = df["accession"].map(netgpi).fillna(False).astype(bool)
    df["gpi_curated"] = df["gpi_curated"].astype(bool)
    df["gpi_any"] = df["gpi_curated"] | df["gpi_predicted"]
    df["gpi_source"] = np.select(
        [df.gpi_curated & df.gpi_predicted, df.gpi_curated, df.gpi_predicted],
        ["both", "curated", "predicted"], default="")

    # ---- transmembrane evidence ----
    uni_tm = df["n_tm_segments"] > 0
    dtm_tm = df["dtm_surface_call"]
    df["tm_sources"] = uni_tm.astype(int) + dtm_tm.astype(int)

    # The measured extracellular length. Take the larger of the two sources:
    # 2,236 UniProt TM proteins carry no topological domains at all, so their
    # curated value is 0 purely for lack of annotation, not for lack of domain.
    df["longest_extracellular"] = df[
        ["longest_extracellular_len", "dtm_longest_outside"]].max(axis=1)
    df["extracellular_source"] = np.where(
        df.dtm_longest_outside > df.longest_extracellular_len,
        "deeptmhmm", np.where(df.longest_extracellular_len > 0, "uniprot", "none"))

    # "Length 0" has two very different meanings. For 2,236 UniProt TM proteins
    # no topological domain is annotated at all, so 0 means "nobody measured it",
    # not "the domain is short". Keeping those apart matters: one is a weak
    # candidate, the other is an unmeasured one.
    df["length_unknown"] = (df.extracellular_source == "none") & (df.tm_sources > 0)
    long_enough = df["longest_extracellular"] >= args.min_extracellular

    # D32: a curated GPI anchor outranks a predictor-ONLY transmembrane call.
    # The C-terminal GPI attachment signal is a hydrophobic stretch, and reading
    # it as a TM helix is a documented failure of TM predictors — it is why the
    # runbook insisted on a separate GPI route. Without this, CEACAM5 (a
    # GPI-anchored ADC target) lands in S2 on a predicted helix UniProt does not
    # annotate. Where UniProt ITSELF annotates a TM alongside the GPI anchor
    # (BST2, RAET1G), that is real dual anchoring and the TM tier stands.
    gpi_beats_tm = df.gpi_curated & (df.n_tm_segments == 0)

    df["tier"] = np.select(
        [
            df.gpi_any & gpi_beats_tm,
            (df.tm_sources == 2) & long_enough,
            (df.tm_sources == 1) & long_enough,
            df.gpi_any & (df.tm_sources == 0),
            (df.tm_sources >= 1) & ~long_enough,
        ],
        ["S3", "S1", "S2", "S3", "Flagged"],
        default="-",
    )
    df["is_candidate"] = df.tier.isin(["S1", "S2", "S3", "Flagged"])

    # D21, enforced rather than left to emerge: S1 means two independent
    # sources agreed. A protein the predictor never saw cannot satisfy that,
    # so it is demoted explicitly and the demotion is counted.
    impossible_s1 = (df.tier == "S1") & df.predictor_missing
    df.loc[impossible_s1, "tier"] = "S2"
    n_demoted = int(impossible_s1.sum())

    # ---- D31: compartment filter ----
    loc = df["subcellular_location"].fillna("").str.lower()
    df["pm_uniprot"] = loc.str.contains(PLASMA_PAT, regex=True)
    df["internal_uniprot"] = loc.str.contains(INTERNAL_PAT, regex=True)
    pm_genes = hpa_plasma_membrane_genes()
    df["pm_hpa"] = df["accession"].map(
        gmap.drop_duplicates("accession").set_index("accession")["ensembl_gene"]
    ).isin(pm_genes)
    df["plasma_membrane_evidence"] = df.pm_uniprot | df.pm_hpa
    df["compartment"] = np.select(
        [df.plasma_membrane_evidence, df.internal_uniprot],
        ["plasma_membrane", "internal_only"], default="unannotated")

    # Excluded only on POSITIVE evidence of being unreachable. A protein with no
    # compartment annotation is kept: dropping it would penalise it for being
    # under-studied, which is where a novel target is most likely to be.
    # D32: a curated cell-surface GPI protein is not excluded by compartment
    # text that happens to mention an internal membrane (this was dropping ART1).
    excluded = ((df.compartment == "internal_only") & df.is_candidate
                & ~df.gpi_curated)
    df.loc[excluded, "tier"] = "internal"
    df.loc[excluded, "is_candidate"] = False
    n_internal = int(excluded.sum())

    # D20: length is meaningless for a GPI-only protein.
    df["extracellular_reported"] = np.where(
        df.tier == "S3", NOT_APPLICABLE, df.longest_extracellular.fillna(0).astype(int).astype(str))

    # ---- external evidence, carried from Step 2 ----
    ext = (gmap.groupby("accession")
           .agg(cspa=("cspa", "any"), surfy=("surfy_surface", "any"),
                ensembl_gene=("ensembl_gene", "first"),
                hpa_symbol=("hpa_symbol", "first"))
           .reset_index())
    df = df.merge(ext, on="accession", how="left")

    out = df[[
        "accession", "symbol", "hpa_symbol", "ensembl_gene", "tier",
        "tm_sources", "n_tm_segments", "dtm_n_tm", "dtm_type",
        "longest_extracellular", "extracellular_reported", "extracellular_source",
        "gpi_any", "gpi_source", "has_signal_peptide", "dtm_has_signal",
        "predictor_missing", "length_unknown", "compartment",
        "plasma_membrane_evidence", "cspa", "surfy", "is_candidate",
    ]].sort_values(["tier", "accession"], kind="mergesort")
    out.to_csv(PROC / "surface_candidates.csv", index=False)

    # ------------------------------------------------------------- report --
    cand = out[out.is_candidate]
    counts = out.tier.value_counts()
    lines = [
        "STEP 3c — surface candidate set",
        "Generated by src/step3_build_candidates.py — do not edit by hand.",
        f"predictor run: {args.predictor}   minimum extracellular: "
        f"{args.min_extracellular} aa (D20)   NetGPI: {netgpi_src}",
        "",
        f"proteins in                              : {len(out):,}",
        f"surface candidates                       : {len(cand):,}",
        "",
        "tier counts:",
    ]
    for tier in ["S1", "S2", "S3", "Flagged", "-"]:
        label = {"-": "not a candidate"}.get(tier, tier)
        lines.append(f"  {label:<38} : {int(counts.get(tier, 0)):,}")

    lines += [
        "",
        f"excluded by the compartment filter (D31)  : {n_internal:,}",
        "  (positive evidence of an internal membrane only; no plasma-membrane",
        "   evidence from UniProt or HPA. Unannotated proteins are kept.)",
        "",
        "compartment of surviving candidates:",
    ]
    for comp, n in cand.compartment.value_counts().items():
        lines.append(f"  {comp:<38} : {n:,}")
    lines += [
        "",
        "source agreement on transmembrane:",
        f"  both UniProt and DeepTMHMM             : {int((out.tm_sources == 2).sum()):,}",
        f"  UniProt only                           : "
        f"{int(((out.n_tm_segments > 0) & (out.dtm_n_tm == 0)).sum()):,}",
        f"  DeepTMHMM only                         : "
        f"{int(((out.n_tm_segments == 0) & (out.dtm_n_tm > 0)).sum()):,}",
        f"  disagreements (exactly one source)     : {int((out.tm_sources == 1).sum()):,}",
        "",
        "GPI evidence (D22):",
        f"  curated only                           : {int((out.gpi_source == 'curated').sum()):,}",
        f"  predicted only                         : {int((out.gpi_source == 'predicted').sum()):,}",
        f"  both                                   : {int((out.gpi_source == 'both').sum()):,}",
        "",
        "where the extracellular length came from:",
    ]
    for src, n in out[out.is_candidate].extracellular_source.value_counts().items():
        lines.append(f"  {src:<38} : {n:,}")

    flagged = out[out.tier == "Flagged"]
    lines += [
        "",
        "Flagged tier, split by reason:",
        f"  extracellular domain genuinely short   : "
        f"{int((~flagged.length_unknown).sum()):,}",
        f"  length never measured by either source : "
        f"{int(flagged.length_unknown.sum()):,}",
    ]
    lines += [
        "",
        f"demoted S1 -> S2, predictor never ran    : {n_demoted:,}",
        f"candidates with no predictor result      : "
        f"{int(cand.predictor_missing.sum()):,}  (capped at S2, flagged)",
        f"candidates with CSPA evidence            : {int(cand.cspa.fillna(False).sum()):,}",
        f"candidates in the SURFY surfaceome       : {int(cand.surfy.fillna(False).sum()):,}",
    ]

    report = "\n".join(lines)
    print(report)
    (RESULTS / "step3_counts.txt").write_text(report + "\n")
    print(f"\nWrote data/processed/surface_candidates.csv ({len(out):,} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
