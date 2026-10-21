# Is S the right metric? A review of what the field uses

Written 2026-10-21, in response to the question "is S the best metric?".
This is a reference memo, not a decision. Nothing in the pipeline changed
because of it. Any change it motivates must be recorded in `docs/decisions.md`
first, per rule 3.

## What S is

```
S = log2(median tumour TPM + 1) − log2(highest healthy tissue TPM + 1)
```

A difference of logs, so S is the **log of a ratio**: tumour expression over
the single worst healthy tissue. S = 2 means roughly four times more in the
tumour than anywhere healthy. It carries no threshold (D45) and is reported
in four tissue-exclusion variants (D46).

In the field's vocabulary this is a **therapeutic-window** measure, and the
window is the central constraint in every ADC and CAR framework reviewed below.
So S is measuring the right thing. The question is whether it measures enough.

## What published frameworks actually do

No framework found ranks on a single expression metric. All combine a window
measure with hard filters on other axes.

| Source | How normal tissue is summarised | What else is required |
| --- | --- | --- |
| **ImmunoTar** (Bioinformatics 2025) | GTEx + pediatric Evo-Devo + healthy proteomics, folded into a weighted sum of ~40 rescaled, "curved" features | Localisation (CIRFESS, COMPARTMENTS, UniProt), DepMap dependency, reagent availability (TTD, ADCdb). Positives/negatives from drug databases; scored by mean average precision |
| **Perna et al.**, Cancer Cell 2017 (AML) | Transcriptomics + surface proteomics of normal tissues; explicit exclusion of normal CD34+CD38−, T cells and vital tissues | Concluded **no single antigen matched CD19's profile**, and moved to combinatorial pairs |
| **ADC target review** ("Evolving Strategies…") | *NormScore*: median per tissue across **45 GTEx normal tissues** | Plots tumour expression as fold change over all GTEx normals. Reports CD70 at ~64-fold tumour/normal; NormScores of marketed targets span 0–25 |
| **Computational ADC selection, breast** (PMC4823055) | Comparison against 10 major organs | `max(normal/cancer) < 2`, `mean(cancer/normal) > 2`, and **overlapping coefficient < 0.6** between the two distributions |
| **In-silico ADC repertoire** (PLOS ONE 2024) | Hard exclusion on **13 critical normal tissues** | Quasi H-score ≥ 150/300 from IHC in ≥1 of 20 tumour types |
| **CAR target criteria** (Hofe 2023, iLABMED) | — | Ranked criteria include **expression level and percent of antigen-positive cells**, and **number of patients with antigen-positive cancers** |
| **Tissue-specificity benchmark** (Kryuchkova-Mostacci & Robinson-Rechavi, Brief Bioinform 2017) | — | Of all tissue-specificity indices, **tau performs best**; Gini comparable; TSI, SPM and EE weaker |

## Where S is strong

- **It is the quantity that governs toxicity.** Every framework above treats
  the tumour-to-normal window as the primary constraint.
- **Maximum, not mean, is the conservative reading** and matches practice:
  both ADC papers exclude on *any* critical tissue rather than on an average.
  Step 9 confirmed this empirically — the 90th-percentile summary is the least
  conservative setting in the entire grid (39% overlap with the baseline).
- **It is interpretable and threshold-free.** A number a reviewer can argue
  with beats a composite score whose weights nobody can defend. ImmunoTar's
  weights, by contrast, are fitted to known drug targets — which would have
  been circular here (the reason D45 rejected calibrating a cut-off on
  approved targets).

## Where S is weak

**1. A ratio ignores absolute level — and antigen density is absolute.**
CAR-T activity has a hard density floor: 28ζ CARs lose control below roughly
1,000–2,000 molecules per cell, and CAR design is tuned to the 1,000–6,000
range (Majzner et al.). A gene at 6 TPM tumour / 0.5 TPM healthy earns the same
S as one at 600 / 50, but only one of them is drugged. Our 5 TPM floor (D43) is
a crude guard on this; S itself is scale-free. **This is the largest conceptual
gap.**

**2. A median hides how many tumours carry the target.** The field asks for
expression "in most if not all tumour cells … and in a large fraction of
patients". We compute this (`frac_tumours_above_10tpm`, `tumour_tpm_q25`) and
show it on every target, but it does not order the list. In our own top 10:

| gene | S | median TPM | lower quartile | % tumours > 10 TPM |
| --- | --- | --- | --- | --- |
| CLDN6 | 6.40 | 158 | 93 | 97% |
| TSPAN10 | 6.20 | 221 | 122 | 100% |
| **DCSTAMP** | **4.32** | **94** | **0.95** | **66%** |
| TRPM1 | 4.08 | 94 | 46 | 99% |

DCSTAMP ranks 5th on S while a quarter of thyroid tumours carry essentially
none of it. S cannot see that; the data already can.

**3. One tissue sets the whole score.** The maximum is a point estimate from a
single tissue's median — it treats "one tissue slightly above" the same as
"ten tissues above", and it is sensitive to one noisy tissue. We compute
`healthy_breadth_above_10tpm` and `gtex_tau` and rank on neither. Of our ten
finalists, eight have breadth 0, UPK3B has 4 and DBH 2.

**4. Two point estimates, not two distributions.** The breast ADC work uses an
**overlapping coefficient** between the tumour and normal expression
distributions. That is strictly more informative than comparing a median to a
maximum, and it is the one methodological upgrade here that is a genuinely
better *form* of the same measurement rather than an extra criterion.

**5. Everything downstream of RNA.** Bulk RNA, no purity correction (D61), and
protein is reported but never ranked on (D48). Already stated as limitations;
unchanged by this review.

**6. We deliberately dropped the critical-tissue list that both ADC papers
use** (D42). That remains defensible — I judged she could not author such
a list — but it is a divergence from published practice, not an alignment with
it, and the report should say so.

## Assessment

S is the right **primary** axis and should stay. It is not sufficient alone,
but nothing in the literature uses a sufficient single metric either — the
field's pattern is *one window measure plus hard filters on density,
prevalence and critical tissues*.

Our pipeline already follows that pattern, arranged differently: eligibility
filters (tier, FDR, fold change, 5 TPM) → rank by S → report everything else
beside each target without letting it move the rank (D52). The honest gap is
that **prevalence is computed, displayed, and ignored by the ordering**, and
prevalence is the criterion the literature is most insistent about.

## Options, in order of what they would cost

1. **Change nothing; cite this memo in the report.** Defensible as it stands.
2. **Metric-sensitivity analysis, in the shape of Step 9.** Re-rank under
   alternative metrics — S, S with a prevalence requirement, overlapping
   coefficient, tau-weighted, ImmunoTar-style weighted sum — and report which
   targets survive all of them. Adds information, overturns nothing, and keeps
   the pre-registered rule intact.
3. **Promote prevalence to an eligibility filter** (e.g. require
   `frac_tumours_above_10tpm ≥ 0.70`), recorded as a new decision.
4. **Replace S** with a composite or distribution-overlap score.

Option 4 carries a specific cost worth naming: D52 was written and dated
*before* the first ranking ran, and check 7.4 verifies that. Changing the
ranking metric now — after seeing which genes won — forfeits that guarantee.
Option 2 is the version that buys the information without spending the
pre-registration.

## Sources

- ImmunoTar — https://academic.oup.com/bioinformatics/article/41/3/btaf060/8008580
- Perna et al., Cancer Cell 2017 — https://www.cell.com/cancer-cell/pdf/S1535-6108(17)30408-7.pdf
- Kryuchkova-Mostacci & Robinson-Rechavi, Brief Bioinform 2017 — https://academic.oup.com/bib/article/18/2/205/2562739
- Majzner et al., antigen density — https://escholarship.org/content/qt80m5z7zj/qt80m5z7zj.pdf
- Programming CAR T tumour recognition (review) — https://aacrjournals.org/cancerdiscovery/article/13/4/829/724992
- Computational selection of ADC targets, breast — https://pmc.ncbi.nlm.nih.gov/articles/PMC4823055/
- In-silico ADC target repertoire, PLOS ONE 2024 — https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0308604
- Target selection for CAR-T, Hofe 2023 — https://onlinelibrary.wiley.com/doi/full/10.1002/ila2.8
- SCAN-ACT, single-cell target discovery — https://pmc.ncbi.nlm.nih.gov/articles/PMC12351953/
