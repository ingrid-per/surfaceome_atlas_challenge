# Surface Target Pipeline — Methodology

Oct 20, 2026

## How to use this runbook

Work one step per Claude Code session, and never let it choose a threshold or a list silently. This runbook covers the core of the brief: surface proteins up in tumors, absent from healthy tissue, using Protein Atlas, TCGA and a TMHMM-type model. Extensions come after the core is done.

**Six working rules, and why each exists**

1. **Plan first.** Start every step in plan mode and read the plan before approving. *Why:* you learn what will happen, and you catch wrong assumptions before code exists.
2. **Every number comes from a script.** Results are written to files by code in `src/`, never typed by hand. *Why:* the numbers can be re-run and checked, and nothing is invented.
3. **Every choice goes in `docs/decisions.md` before you see the result it affects.** Record what, the value, why, and the date. *Why:* a threshold picked after seeing the output is a way of fitting the answer. Written first, it is a defensible method.
4. **Print counts at every stage.** Rows in, rows out, and why rows were dropped. *Why:* most pipeline bugs show up as a count that makes no sense.
5. **Commit to git after each step.** *Why:* you can go back, and each result ties to a version of the code.
6. **You read the code that makes decisions.** Ask Claude Code to explain any line in plain English. *Why:* in the interview you must be able to defend every choice yourself.

**Where this runs:** on your own computer. My workspace here cannot reach the data hosts (they returned errors), and Claude Code runs wherever you start it, so start it in a terminal on your machine.

**File 1: `CLAUDE.md` in the project root.** Claude Code reads it at the start of every session, so these rules apply without repeating them.

```markdown
# Project: surface-protein targets in cancer
Goal: top 10 surface proteins up in tumors and absent from healthy tissue.
Rules:
- Every number in a table or report must come from a script in src/ that writes to results/. Never type numbers by hand.
- Never pick a threshold, gene list or tissue list silently. Propose it, wait for my OK, then record it in docs/decisions.md (what, value, why, date).
- Record every dataset's name, version or date, source URL and checksum in data/MANIFEST.md.
- At each stage print and save counts: rows in, rows out, and the reason rows were dropped.
- Write a small automated check for every step, listed in docs/checks.md.
- data/raw is read-only. Never delete or edit raw data.
- If a result looks too good or too empty, stop and tell me what to check.
- Explain non-obvious code in plain-English comments.
```

**File 2: the prompt template.** Paste it at the start of each step and fill the brackets.

```markdown
Step: [name]
Goal: [one sentence]
Inputs: [files]
Outputs: [files and a summary table]
Decisions to ask me about before coding: [list]
Checks that must pass: [list]
Do not: change earlier steps' outputs, invent data, or skip the counts.
First: enter plan mode and show me the plan.
```

**Order of work.** Do steps 0 to 8 first: that is a complete answer to the question. Steps 9 to 11 add robustness, your own judgment, and the write-up.

## Step 0: Set up the project

Set up a clean, versioned folder with a decision log, so every later result can be traced to code and to a recorded choice.

**What to do**

1. Create a project folder with this layout: `data/raw/`, `data/processed/`, `src/`, `results/`, `docs/`.
2. Create `CLAUDE.md` (from the section above), `docs/decisions.md`, `docs/checks.md` and `data/MANIFEST.md`.
3. Run `git init` and make a first commit.
4. Create a Python environment with an environment manager (conda or uv). Packages to start with: pandas, numpy, scipy, statsmodels, matplotlib, biopython, requests.
5. Start Claude Code in that folder and ask it to check the machine: free disk space, Python version, and whether it has internet access.

**Why these choices**

- **Python, one language.** It covers the statistics now and the single-cell tools later (scanpy, CELLxGENE). R's limma is a good alternative, but mixing two languages makes the pipeline harder for you to follow and defend.
- **Raw data kept separate and read-only.** Any later surprise can be traced back to the untouched originals.
- **A manifest of dataset versions.** HPA and TCGA-derived resources change between releases, and the same code can give different results on different versions.
- **Disk space check first.** Expression matrices for thousands of samples can be large (several gigabytes is plausible, which I have not verified). Ask Claude Code to check file sizes before downloading.

**Checks before moving on**

- `git log` shows one commit, and the folder layout exists.
- Claude Code can report disk space and confirm it can reach the internet. If it cannot, you download files yourself and put them in `data/raw/`.

**First entries for `docs/decisions.md`** (write these yourself, they are the ones that matter most):

| Decision | Value | Why |
| --- | --- | --- |
| Question being answered | Surface proteins up in tumors, absent from healthy tissue | the brief, restated |
| Healthy reference | HPA normal tissue plus GTEx | He named HPA; GTEx adds samples and shares a pipeline with TCGA (see Step 1) |
| Scope | Pan-cancer, reporting the best cancer type per gene | His wording says "in cancer" without naming one |

Tell in the write-up that GTEx is an addition to the sources the brief listed. It is better to say so than to let him find out.

## Step 1: Get and version the datasets

Download six resources, record exactly which version you got, and confirm each one opens and has plausible sizes. I verified in this chat that most of these resources exist (UniProt is the one I did not check); I did not verify exact file names or download links, so Claude Code should find them on the official pages and you confirm before downloading.

| Resource | What it gives you | Why you need it |
| --- | --- | --- |
| UniProt (reviewed human proteome) | Protein sequences, curated topology, signal peptides, lipid anchors | The input for the surface filter in Step 3 |
| Human Protein Atlas (HPA) downloads | Normal-tissue RNA (about 50 tissues), antibody staining (IHC), gene ID table | The healthy side he named; also the protein-level check |
| UCSC Xena Toil hub (TCGA plus GTEx) | Tumor, matched normal and GTEx normal RNA, processed in one pipeline | The tumor side he named; one pipeline avoids processing differences between tumor and normal |
| CPTAC via LinkedOmicsKB or the Proteomic Data Commons | Tumor and adjacent-normal protein levels, 10 cancer types | Protein evidence for tumors (Step 6) |
| GTEx proteome (Jiang et al., Cell 2020) | Protein levels for 12,627 genes across 32 normal tissues | Protein evidence for healthy tissue (Step 6) |
| Cell Surface Protein Atlas (CSPA) | 1,492 human proteins with experimental evidence of surface expression | A check that predicted surface proteins reach the surface |

**Why these choices**

- **HPA and TCGA are his two named sources**, so they are the backbone. Everything else supports them.
- **GTEx is my addition.** TCGA has few matched normal samples for many cancer types, and GTEx has many healthy samples across tissues. The Toil resource reprocesses both together, which reduces differences caused by software rather than biology.
- **A protein layer is added because RNA is only a proxy.** The GTEx proteome paper itself reports tissue-enriched RNA with no matching protein enrichment, which is the exact worry.
- **CSPA is experimental evidence.** TMHMM-type tools predict, and CSPA measured.

**What Claude Code should do**

1. Find each resource's official download page and show you the exact file names, sizes and release versions before downloading anything.
2. Download into `data/raw/`, then record name, version, date, URL and checksum in `data/MANIFEST.md`.
3. Open each file and print its shape (rows and columns), its first few rows, and the number of samples per cancer type and tissue.

**Checks that must pass**

- Every file is in the manifest with a version and a checksum.
- Sample counts per cancer type are printed, and none is zero.
- You have looked at the numbers and they are plausible. TCGA has thousands of tumors, so a count of 12 means something went wrong.

**Prompt for Claude Code**

```markdown
Step 1: download and version datasets.
Goal: get UniProt human reviewed proteome, HPA downloads, the UCSC Xena Toil TCGA+GTEx gene expression, CPTAC protein data (LinkedOmicsKB or PDC), the GTEx proteome tables, and the Cell Surface Protein Atlas.
Inputs: official download pages only.
Outputs: files in data/raw/, data/MANIFEST.md, results/step1_inventory.csv (file, rows, columns, samples per group).
Decisions to ask me about: which exact files and versions, and total download size.
Checks that must pass: every file in the manifest with checksum; sample counts printed.
Do not: download anything before showing me file names and sizes.
First: enter plan mode and show me the plan.
```

## Step 2: Harmonize gene and protein IDs

Build one mapping table that links every dataset to the same gene, because mismatched IDs silently drop or duplicate genes and spoil every later step.

**The problem.** The same gene has several names in different files. Expression files often use Ensembl gene IDs (sometimes with a version suffix), HPA uses gene names, UniProt uses protein accessions, and older papers use retired symbols. Joining on the wrong one loses genes without any error message.

**What to do**

1. Build a master table: Ensembl gene ID, official gene symbol, UniProt accession. HPA's gene table is a convenient starting point, since it lists these together (from memory, so confirm the columns).
2. Strip version suffixes from Ensembl IDs before joining.
3. Keep protein-coding genes with a reviewed UniProt entry. *Why:* the surface filter in Step 3 needs a curated protein sequence, and non-coding RNAs cannot be surface proteins.
4. Resolve one-to-many and many-to-one matches by a written rule, and log each case.
5. Save the unmapped genes to a file so you can see what was lost.

**Checks that must pass**

- Percent of genes mapped in each dataset is printed. Large losses (over a few percent of protein-coding genes) need an explanation.
- No gene symbol maps to two Ensembl IDs without a recorded rule.
- Spot-check about 20 genes by hand across all tables: for example EGFR, ERBB2, MSLN, TACSTD2 and GAPDH should all be found with consistent IDs.

**Prompt for Claude Code**

```markdown
Step 2: harmonize IDs.
Goal: one master mapping Ensembl gene ID <-> gene symbol <-> UniProt accession, used by every later step.
Inputs: files in data/raw/ from Step 1.
Outputs: data/processed/gene_map.csv, results/step2_unmapped.csv, results/step2_mapping_report.txt.
Decisions to ask me about: how to resolve one-to-many matches; whether to keep non-reviewed UniProt entries.
Checks that must pass: mapped percentage per dataset printed; spot-check table for EGFR, ERBB2, MSLN, TACSTD2, GAPDH.
Do not: drop genes without logging them.
First: enter plan mode and show me the plan.
```

## Step 3: Build the surface-protein candidate set

Keep proteins that have a part outside the cell, using two independent sources (curated annotation and a predictor) plus a separate route for lipid-anchored proteins, and label how confident each call is.

**What a topology tool does.** TMHMM and its newer deep-learning version DeepTMHMM read a protein's amino acid sequence and predict which stretches cross the cell membrane and which parts sit inside or outside the cell. The brief said "TMHMM-type", so DeepTMHMM is a fair choice. It predicts from sequence alone, so it cannot prove a protein actually reaches the surface.

**What to do**

1. **Curated annotation.** From UniProt, take each protein's annotated transmembrane segments, "extracellular" topological domains, signal peptide, and lipid-anchor (GPI) information.
2. **Predictor.** Run DeepTMHMM on the sequences of the reviewed human proteome. Check how to run it (web service or package) when you get to this step, since I did not verify its current access.
3. **Lipid-anchored proteins.** Add a signal-peptide predictor (SignalP) and a GPI-anchor predictor (for example NetGPI). *Why:* TMHMM-type tools miss proteins attached by a lipid anchor instead of a membrane-spanning segment. Mesothelin is an example, so a strict TMHMM filter would drop a well-known target.
4. **Outside evidence.** Mark proteins found in the Cell Surface Protein Atlas (experimental) and in the SURFY predicted surfaceome of 2,886 proteins.
5. **Assign a tier to each protein.**

| Tier | Rule | Confidence |
| --- | --- | --- |
| S1 | Transmembrane with an extracellular part of at least 20 residues, called by both UniProt and DeepTMHMM | High |
| S2 | Same, but called by only one of the two | Medium |
| S3 | GPI-anchored by prediction or annotation | Medium |
| Flagged | Surface call but extracellular part under 20 residues | Low |

**Why these choices**

- **Two independent sources, with agreement shown.** UniProt is curated but incomplete for little-studied proteins. A predictor covers every sequence but can be wrong. Agreement is stronger evidence than either alone.
- **The 20-residue cutoff is my assumption, not a standard.** The reasoning is that a tiny stub is unlikely to hold an accessible binding site. Test 10, 20 and 50 in Step 9 and report how much the answer changes.
- **Tiers, not a yes/no.** It keeps uncertain proteins visible instead of hiding them.

**Checks that must pass**

- The candidate set should land in the low thousands. For scale, SURFY's prediction is 2,886 proteins and CSPA's experimental set is 1,492. A number far outside about 1,500 to 4,500 needs investigation.
- Positive controls appear: EGFR and ERBB2 in S1, MSLN in S3.
- Negative controls are absent: GAPDH, ACTB and histones should not be surface candidates.
- Print how many proteins land in each tier, and how many UniProt and DeepTMHMM disagree on.

**Prompt for Claude Code**

```markdown
Step 3: surface candidate set.
Goal: tier every reviewed human protein as S1, S2, S3, Flagged or not-surface, using UniProt topology, DeepTMHMM, SignalP and a GPI predictor, with CSPA and SURFY as external evidence.
Inputs: data/raw/ UniProt, CSPA, SURFY lists; data/processed/gene_map.csv.
Outputs: data/processed/surface_candidates.csv (gene, tier, sources agreeing, extracellular length, GPI flag, CSPA flag), results/step3_counts.txt.
Decisions to ask me about: the minimum extracellular length (propose 20), how to run DeepTMHMM, which GPI predictor.
Checks that must pass: candidate count between about 1,500 and 4,500; EGFR/ERBB2 in S1, MSLN in S3; GAPDH/ACTB absent; tier counts printed.
Do not: hard-code any gene list.
First: enter plan mode and show me the plan.
```

## Step 4: Tumor versus normal at RNA level

For each surface candidate and each cancer type, test whether tumors express it significantly more than the matching normal tissue, and report how large the difference is.

**What to do**

1. **Match each TCGA cancer type to its normal tissue.** Write a small mapping file, for example lung adenocarcinoma to lung and colon adenocarcinoma to colon. *Why:* comparing a tumor with the wrong tissue creates fake signal. This is a scientific judgment, so you review the file by hand and record it in `docs/decisions.md`.
2. **Transform the values.** Use log2(TPM + 1). *Why:* TPM values are heavily skewed, the log makes differences comparable across the range, and the +1 avoids taking the log of zero.
3. **Test each gene** with a Wilcoxon rank-sum test (also called Mann-Whitney U), tumor against normal, within each cancer type. *Why:* it makes no assumption that expression is normally distributed, it resists outliers, and it is easy to explain. limma-voom is more powerful for small groups but harder to explain, so use it only as an optional cross-check.
4. **Correct for multiple testing** with the Benjamini-Hochberg false discovery rate (FDR), separately per cancer type. *Why:* you are testing thousands of genes, and without correction about 5 in every 100 null genes would look significant at p < 0.05.
5. **Report the effect size** as the difference in median log2 expression, tumor minus normal (log2 fold change). *Why:* with thousands of samples, tiny differences become statistically significant, so significance alone does not mean importance.

**Default thresholds (my proposals, to test in Step 9):** FDR below 0.01 and log2 fold change of at least 2, which is a fourfold increase. These are conventions, not laws, so record them as decisions and show how the answer changes when they move.

**Output per gene and cancer type:** number of tumor and normal samples, median expression in each, log2 fold change, p-value, FDR.

**Known weaknesses to state in your report**

- TCGA tumors come from surgery, and GTEx normals often from post-mortem collection. The Toil reprocessing removes software differences, not differences in sample handling.
- Bulk tissue mixes tumor, stromal and immune cells, so a gene can look tumor-high because of the surrounding cells.
- Some cancer types have few matched normals. Print the counts so you can see which comparisons are weak.

**Checks that must pass**

- Plot a histogram of p-values for each cancer type. A healthy test shows a spike near zero and roughly flat elsewhere. A different shape suggests a bug or a batch problem.
- Print how many genes pass the thresholds per cancer type.
- Pick a few well-documented tumor-up or tumor-down genes from the literature (you choose them from a published source, not from memory) and confirm the direction matches.
- Each comparison has at least a stated minimum number of samples in both groups. Propose that minimum and record it.

**Prompt for Claude Code**

```markdown
Step 4: tumor vs normal RNA.
Goal: for every surface candidate and each TCGA cancer type, compute median log2(TPM+1) in tumor and matched normal, log2 fold change, Wilcoxon rank-sum p-value, and BH-FDR per cancer type.
Inputs: Toil TCGA+GTEx expression, data/processed/surface_candidates.csv, data/processed/gene_map.csv.
Outputs: results/step4_tumor_vs_normal.csv, results/step4_pvalue_histograms.png, results/step4_counts.txt.
Decisions to ask me about: the TCGA-to-normal-tissue mapping file, minimum samples per group, thresholds (propose FDR<0.01, log2FC>=2).
Checks that must pass: p-value histograms look sensible; counts printed per cancer type; sample numbers per group printed.
Do not: use normal tissue from a different organ than the mapping file says.
First: enter plan mode and show me the plan.
```

## Step 5: Score how absent a gene is from healthy tissue

Score each gene by its highest expression across vital healthy tissues, never by the average, because one exposed organ is enough to make a target unsafe.

**What to do**

1. **Build the normal-tissue table.** For each gene, collect expression in every tissue from HPA (about 50 tissues) and from GTEx. Compute the maximum per source, and note which tissue gives it.
2. **Define a vital-organ list** in `docs/decisions.md`, before you see results. A starting point to review, not a standard: heart, brain regions, lung, liver, kidney, and bone marrow or blood. Mark other tissues as "tolerated" and record why. Tissues such as testis are often treated as more tolerable, but that is a judgment you make and defend.
3. **Compute the specificity score** for each gene and cancer type:

```latex
S = \log_2(\text{tumor median TPM} + 1) - \log_2(\text{max vital-tissue expression} + 1)
```

4. **Apply a minimum tumor expression** (propose 5 TPM and record it). *Why:* a gene that is absent in healthy tissue but barely expressed in the tumor is not a usable target.
5. **Report, but do not score, expression in tolerated tissues** as a separate column.
6. **Optional second view:** the tau tissue-specificity index (0 means expressed everywhere, 1 means one tissue only). Use it as a cross-check only.

**Why these choices**

- **Worst case, not average.** A protein high in heart but low everywhere else has a low average and is still dangerous. A 2026 colorectal-cancer target framework made the same choice, subtracting HPA and GTEx normal scores individually instead of averaging them.
- **HPA and GTEx kept separate, then the larger value taken.** The two sources use different units and processing (HPA reports normalized TPM, called nTPM). Treating them as one scale would be a hidden assumption, so use the more cautious of the two.
- **Vital-organ weighting is a human judgment.** That is where your bench background matters, so the list is yours to decide, and the report should state it.

**Checks that must pass**

- A protein expected to be everywhere on cell surfaces (for example B2M) should score near zero or negative. If it scores high, something is wrong.
- For ten genes, open the tissue table by hand and confirm the maximum and the tissue it comes from.
- Print how many genes have a maximum vital-tissue value above the detection floor, to see how strict the filter is.

**Prompt for Claude Code**

```markdown
Step 5: normal-tissue specificity.
Goal: for each surface candidate, compute max expression across vital normal tissues in HPA and in GTEx separately, then the specificity score S = log2(tumor median TPM+1) - log2(max vital expression+1), plus the tau index as a cross-check.
Inputs: results/step4_tumor_vs_normal.csv, HPA normal-tissue RNA, GTEx tissue medians.
Outputs: results/step5_specificity.csv (gene, cancer type, S, max vital tissue and value for each source, tolerated-tissue columns, tau).
Decisions to ask me about: the vital-organ list, the tolerated list, the minimum tumor TPM (propose 5).
Checks that must pass: B2M scores low; 10 genes hand-verified; counts printed.
Do not: average across tissues.
First: enter plan mode and show me the plan.
```

## Step 6: Add protein-level evidence

Check each RNA-based hit against real protein measurements, and sort genes into evidence tiers, treating "protein not detected" as a flag and not as proof of absence.

**Why this step exists.** RNA is a proxy for protein. The GTEx proteome paper reports tissue-enriched RNA with no matching protein enrichment. Surface proteins are also hard to detect by mass spectrometry because they are hydrophobic, and one CPTAC table lists roughly 4,700 to 6,200 proteins per study against about 20,000 genes. So protein data can confirm a hit but cannot rule one out.

**What to do**

1. **Tumor protein (CPTAC).** Use LinkedOmicsKB's precomputed tumor-versus-normal protein results, or the Proteomic Data Commons for raw data. Mark a gene "protein up in tumor" if it is significantly higher in tumors than in adjacent normal tissue, in the same direction as the RNA.
2. **Normal protein (GTEx proteome).** Mark whether the protein was detected in each vital tissue. Use detected or not detected, or ranks within the same study. *Why:* abundances from different studies are on different scales and cannot be compared directly.
3. **Antibody staining (HPA IHC).** Record staining level and reliability in vital tissues. I have not verified the exact columns of the HPA file, so ask Claude Code to show them first.
4. **Surface evidence (CSPA).** Mark genes with experimental evidence of surface expression.
5. **Assign evidence tiers.**

| Tier | Rule |
| --- | --- |
| 1: protein-confirmed | Surface tier S1, S2 or S3, passes the RNA thresholds, protein up in tumor, and no protein detected in vital normal tissue |
| 2: RNA-only | Passes the RNA thresholds, but protein not measured or not detected (flagged, not dropped) |
| 3: contradicted | Protein clearly detected in a vital normal tissue (mass spectrometry or high staining) |

**Why tiers.** A strict protein filter would delete real targets just because the instrument missed them. Keeping them in tier 2 is honest about what is known.

**Known weaknesses to state in your report**

- CPTAC covers only 10 cancer types, and some cohorts have few normal samples. One table lists 133 breast tumors against 18 adjacent normals.
- Adjacent normal tissue is not truly healthy, because it sits next to a tumor.
- Bulk tissue still mixes cell types.

**Output:** a funnel table with the number of genes surviving each stage: all genes, protein-coding, surface candidates, passing RNA tests, passing specificity, tier 1, tier 2, tier 3. This table is one of your best exhibits.

**Checks that must pass**

- The funnel counts are monotone (each stage no larger than the one before) and every drop has a logged reason.
- For five tier-1 genes, confirm by hand in the original source that the protein evidence is there.
- Tier 2 should not be tiny. If almost every RNA hit has protein data, check that unmeasured genes are not being treated as measured.

**Prompt for Claude Code**

```markdown
Step 6: protein-level evidence tiers.
Goal: attach protein evidence to every RNA-passing candidate and assign tier 1, 2 or 3.
Inputs: results/step5_specificity.csv, CPTAC protein tumor-vs-normal results, GTEx proteome tables, HPA IHC normal-tissue file, CSPA list.
Outputs: results/step6_evidence_tiers.csv, results/step6_funnel.csv (stage, genes remaining, reason for drop).
Decisions to ask me about: what counts as detected in the GTEx proteome, the CPTAC significance rule, how to treat genes not measured.
Checks that must pass: funnel counts monotone; 5 tier-1 genes hand-verified; unmeasured genes counted as RNA-only, not as confirmed.
Do not: compare raw protein abundances across different studies.
First: enter plan mode and show me the plan.
```

## Step 7: Combine, rank and select the top 10

Write the selection rule down before you look at any ranking, then apply it mechanically, so the top 10 cannot be shaped by what you hoped to see.

**The rule (a proposal to confirm and record in `docs/decisions.md` first)**

1. **Eligible:** surface tier S1, S2 or S3, and passing the RNA thresholds from Step 4 in at least one cancer type.
2. **Sort by:** evidence tier first (tier 1 before tier 2; tier 3 excluded), then by the specificity score S from Step 5, highest first.
3. **Limit per gene family:** at most two genes from the same family (for example the claudins), so the list is not ten near-copies.
4. **For each gene, report** its best cancer type, and in how many cancer types it passes.
5. **Take the first 10.**

**Why write it first.** If you pick the rule after seeing which genes it yields, you are tuning to an answer you already like. A rule fixed in advance can be defended as a method, and you can say exactly when you chose it.

**Why evidence tier before score.** A gene with protein confirmation and a slightly lower score is a safer bet than one that only looks good in RNA, given everything Step 6 showed about RNA being a proxy.

**Output table for each of the 10 genes:**

| Column | Meaning |
| --- | --- |
| Gene and UniProt ID | Identity |
| Best cancer type, and number of cancer types passing | Where it is strongest, and how broad |
| Tumor median, log2 fold change, FDR | The tumor-side evidence |
| Highest vital-tissue expression and which tissue | The safety side |
| Specificity score S | The combined number |
| Surface tier, extracellular length | Why it counts as a surface protein |
| Evidence tier | How well protein data supports it |

**Also keep the next five** (ranks 11 to 15). Step 10 uses them in case you drop a finalist.

**Checks that must pass**

- The same code run twice gives the identical list.
- Pick three genes that just missed the cut and confirm they miss because of the rule, not because of a bug.
- The rule written in `docs/decisions.md` has a date earlier than the first ranking file.

**Prompt for Claude Code**

```markdown
Step 7: select the top 10.
Goal: apply the pre-registered selection rule to produce a top 10 and a reserve of 5.
Inputs: results/step5_specificity.csv, results/step6_evidence_tiers.csv, surface_candidates.csv.
Outputs: results/step7_top10.csv, results/step7_reserve.csv, results/step7_ranking_full.csv.
Decisions to ask me about: the gene-family grouping source and the per-family limit; confirm the sort order. Show me the rule text so I can record it in decisions.md before you run anything.
Checks that must pass: deterministic output; three near-miss genes explained.
Do not: look at the ranking before I confirm the rule is recorded.
First: enter plan mode and show me the plan.
```

## Step 8: Sanity checks with controls

Test the pipeline on cases where you already know the right answer, because a method that cannot rediscover known targets cannot be trusted on new ones.

**What to do**

1. **Build a positive-control list from a published source, not from memory.** Use clinical-stage and approved surface targets from a source such as the Cancer Surfaceome Atlas supplement, ChEMBL or Open Targets. Candidates I would expect to see, to check against that source: ERBB2, EGFR, FOLR1, TACSTD2, NECTIN4, MSLN and members of the claudin family. These come from my general knowledge and are unverified.
2. **Build a negative-control list.** Intracellular housekeeping genes (GAPDH, ACTB) and a widely expressed surface protein (B2M).
3. **Measure:**
   - how many positive controls are in the surface candidate set (Step 3);
   - how many pass the RNA thresholds, and at what rank (Steps 4 and 5);
   - how many negative controls are correctly excluded.
4. **Compare your top 10 and top 50** with published lists: the Cancer Surfaceome Atlas (search it by name) and the ImmunoTar tool paper from 2025. Report the overlap, and explain two or three disagreements.

**Why these checks**

- **Recovery of known targets** is the closest thing to ground truth here. It is the evaluation idea from your study guide: test against an outcome you know, not against your own claims.
- **Negative controls** show the filters are strict enough, since a method that passes everything also "recovers" every target.
- **Comparing with published work** shows you know the field and are not claiming something already done.

**How to read the results**

- No positive controls recovered: a bug, most likely in ID mapping or the tissue mapping.
- Every control recovered and nothing else is new: thresholds may be too loose or the controls too easy.
- Blood cancer targets such as CD19 will not appear in solid-tumor comparisons. Say that rather than counting it as a failure.
- High recovery is not proof of finding new, good targets. State that limit.

**Output:** a small table with each control, whether it was recovered, and its rank.

**Prompt for Claude Code**

```markdown
Step 8: controls and comparison.
Goal: measure recovery of known targets and exclusion of negatives, and compare against published lists.
Inputs: all previous results; control lists I will give you from a published source.
Outputs: results/step8_controls.csv (gene, type, in candidate set, passes RNA, rank, tier), results/step8_overlap_with_published.csv, results/step8_summary.txt.
Decisions to ask me about: which source to take positive controls from; how to load the published lists.
Checks that must pass: every control accounted for; blood-cancer-only targets listed separately.
Do not: make up control lists from memory.
First: enter plan mode and show me the plan.
```

## Step 9: Robustness

Re-run the pipeline under many reasonable settings, and report which targets survive nearly all of them and which depend on one choice.

**What to do**

1. **Define the settings to vary** and record them in `docs/decisions.md`:

| Setting | Values to try |
| --- | --- |
| Minimum extracellular length | 10, 20, 50 residues |
| FDR threshold | 0.05, 0.01 |
| Log2 fold change threshold | 1, 2, 3 |
| Healthy reference | HPA only, GTEx only, both |
| Vital-organ list | with and without testis treated as tolerated; with and without brain |
| Normal-tissue summary | maximum, 90th percentile |

2. **Run every combination.** That is a few hundred runs, which is cheap because the steps are vectorized and only the thresholds change.
3. **For each gene, record the share of runs in which it appears in the top 20.** Call that its stability.
4. **Show the result** as a heatmap or a ranked list of genes by stability.

**Why this step**

- **It exposes the weight of your own choices.** Every threshold in this runbook is a convention. If the top 10 barely moves when they change, the answer is robust. If it reshuffles, you know which choice matters and can say so.
- **It turns a weakness into evidence.** Most analyses present one setting and imply it is the only one. Showing the spread is more honest and harder to attack.
- **It answers the obvious question in advance:** "What if you had used different thresholds?"

**How to read it.** Genes present in most runs are your robust core. Genes that appear in only a few runs are fragile: keep them out of the headline list, or mark them clearly.

**Checks that must pass**

- The default setting is one of the runs, and its top 10 matches Step 7 exactly.
- Run count is printed and matches the grid size.
- Stability values lie between 0 and 1, and the heatmap can be reproduced from the saved table.

**Prompt for Claude Code**

```markdown
Step 9: robustness grid.
Goal: re-run Steps 3-7 across a grid of settings and compute each gene's stability (share of runs where it is in the top 20).
Inputs: the pipeline code from Steps 3-7.
Outputs: results/step9_runs.csv (settings and top 20 per run), results/step9_stability.csv, results/step9_stability_heatmap.png.
Decisions to ask me about: the exact grid values (propose the table in the runbook).
Checks that must pass: the default run reproduces Step 7; run count matches grid size.
Do not: change the pipeline's logic, only its parameters.
First: enter plan mode and show me the plan.
```

## Step 10: Review the finalists yourself

Look at the raw evidence for your top 15 as a bench scientist would, and record every place you disagree with the machine, because this judgment is what an AI-only candidate cannot supply.

**What to do**

1. **Take the top 10 plus the reserve of 5** from Step 7.
2. **For each gene, work through this checklist** (open the original sources yourself, do not ask Claude Code to summarize them):
   - HPA antibody staining images for normal vital tissues and for tumors. Is the staining on the cell membrane or in the cytoplasm?
   - Is the antibody well validated? HPA gives a reliability rating.
   - Is there staining in vital tissues that the RNA data hid?
   - How large is the extracellular part, and does UniProt list glycosylation sites that might shield it?
   - Is the target already in clinical trials or approved drugs? Search the trial registry yourself, and note any reported toxicity.
   - Optional: view the AlphaFold predicted structure of the extracellular region, and check its confidence. A long, low-confidence, disordered region is a poor binding site.
3. **Fill in an override log** with one row per gene.

| Gene | Machine rank | My verdict (keep, drop, flag) | Reason | Evidence link |
| --- | --- | --- | --- | --- |
|  |  |  |  |  |

4. **Replace dropped genes from the reserve**, and say so in the log.

**Why this step**

- **It is where your background shows.** Interpreting staining, antibody quality and accessibility is bench judgment. The computation produces candidates, and you decide which ones deserve belief.
- **It guards against quiet pipeline errors** that passed every automated check.
- **It is your answer to "how do you know it works?"** You can point to specific genes where you checked the primary data and changed the outcome.

**Rules for the log**

- You may override the machine, but never silently. Every change has a reason you could say aloud.
- Say what you would test at the bench before betting on each gene. Do not claim any wet-lab validation you have not done.

**Checks that must pass**

- Every one of the 15 genes has a row, a verdict and a link.
- At least a few genes have a specific note from the images (for example "clean membranous signal in tumor, weak in liver"), not generic text.

This step is manual by design, so there is no Claude Code prompt. You can ask Claude Code to build the log table and gather the links for you, but not to make the judgments.

## Step 11: Write the report

Produce a short report that gives the answer first, shows how far to trust it, and says plainly what is yours and what follows published work.

**Contents, in this order (aim for two pages plus figures)**

1. **The top 10 table** from Step 7 after your Step 10 review, with the override log summarized in one line per change.
2. **The funnel** from Step 6: genes remaining at each stage, with the reason for each drop.
3. **One figure:** every surface candidate as a dot, tumor expression on one axis and worst-case vital-tissue expression on the other, with the top 10 labelled. The ideal region is high tumor and low normal. This shows the logic at a glance.
4. **Method, in plain paragraphs:** the definitions you used for "surface", "significantly overexpressed" and "not in healthy tissue", with the thresholds and why you chose each.
5. **Controls and robustness:** how many known targets you recovered (Step 8), and which top-10 genes are stable across settings (Step 9).
6. **Limitations, stated honestly:**
   - RNA is a proxy for protein, and protein data is shallow for membrane proteins.
   - Bulk tissue mixes cell types, and adjacent normal tissue is not truly healthy.
   - TMHMM-type tools predict topology from sequence and miss lipid-anchored proteins unless you add them, as you did.
   - TCGA and GTEx differ in how samples were collected.
   - Nothing here is validated in the lab.
7. **Prior work in two or three sentences:** this follows approaches in published surface-target analyses such as the Cancer Surfaceome Atlas and the 2025 ImmunoTar tool. Say what you added: the evidence tiers, the robustness grid, and your finalist review.
8. **How you used AI:** one short paragraph. You used Claude Code to write and run the pipeline, every number came from saved scripts, you recorded each decision, and you verified outputs against the original sources and controls.

**Why this shape**

- **Answer first.** The question was direct; the answer should be too.
- **Limitations stated by the author read as judgment.** The same limits found
  by the reader read as gaps.
- **Say how AI was used.** A reader will assume it was. Describing how it was
  controlled turns that assumption into evidence of method.

**Build later, only after all of this is done:** the interactive deliverable,
the state-gated single-cell map for the finalists, and the AlphaFold structure
cards. The core above is a complete answer without them.

**Checks that must pass**

- Every number in the report traces to a file in `results/`.
- A reader outside the field can say, in one sentence, what the report found
  and what its main weakness is.

## Glossary

Plain-language meanings of the terms used above, so you can explain each one without notes.

| Term | Plain meaning |
| --- | --- |
| TCGA | The Cancer Genome Atlas: molecular data from thousands of tumors across about 33 cancer types |
| GTEx | A project measuring gene activity in many healthy tissues from donors |
| HPA | Human Protein Atlas: where genes and proteins are found in tissues, by RNA and antibody staining |
| CPTAC | A consortium that measured proteins in tumors by mass spectrometry |
| TPM | Transcripts per million: a way to scale RNA counts so genes and samples can be compared |
| nTPM | HPA's normalized version of TPM, on its own scale |
| log2 fold change | How many doublings separate two values: 1 means twice as high, 2 means four times |
| Wilcoxon rank-sum test | A test of whether one group tends to have higher values than another, using ranks, not assumed shapes |
| p-value | The chance of seeing a difference this large if there were truly none |
| FDR (Benjamini-Hochberg) | A correction so that among genes called significant, the expected share of false calls is controlled when testing thousands |
| Transmembrane segment | A stretch of a protein that crosses the cell membrane |
| Extracellular domain | The part of a protein that sits outside the cell, where a drug or antibody could bind |
| TMHMM / DeepTMHMM | Programs that predict from a protein's sequence which parts are inside, outside or crossing the membrane |
| GPI anchor | A lipid tether that attaches some proteins to the outer membrane without crossing it |
| Signal peptide | A short tag that sends a protein toward the membrane or out of the cell |
| Surfaceome | The set of all proteins on the cell surface |
| IHC | Immunohistochemistry: staining tissue with antibodies to see where a protein is |
| Adjacent normal tissue | Normal-looking tissue taken next to a tumor, not fully comparable with truly healthy tissue |
| Tau index | A number from 0 to 1 showing how restricted a gene is to a single tissue |
| On-target, off-tumor toxicity | Harm caused when a drug hits its intended target on healthy cells |

## Open items to confirm before relying on them

These are the points in the runbook I could not verify, so treat each as a question for Claude Code or for you to check against the official source.

**Verified in this chat by search (the resource exists; details were not downloaded):** CPTAC covers over 1,000 tumors across 10 cancer types, with many matched normals; the GTEx proteome map covers 12,627 genes in 32 tissues; CSPA reports 1,492 human surface glycoproteins; SURFY predicts a 2,886-protein surfaceome; HPA and the UCSC Xena Toil hub provide the expression data described.

**Not verified (from my general knowledge):**

- Exact download file names and URLs for every resource, and total download sizes.
- Which columns the HPA gene table, normal-tissue IHC file and any pathology file contain.
- That the Toil hub contains TCGA matched normals alongside GTEx in one matrix, and the exact dataset name.
- How to run DeepTMHMM today (web service or package), and which signal-peptide and GPI predictors are current.
- UniProt field names for topology, signal peptide and GPI-anchor annotation.
- Which CPTAC cohorts have matched normals and how many.
- The control lists: build them from a published source.

**My own assumptions, not standards (each must be recorded as a decision and tested in Step 9):**

- 20-residue minimum extracellular length.
- FDR below 0.01, log2 fold change of at least 2, and a 5 TPM minimum tumor expression.
- The default vital-organ and tolerated-tissue lists.
- The Wilcoxon test as the main test, with limma-voom as an optional check.
- The evidence-tier rules and the two-per-family limit.

**Sources I opened or saw in search results**

- [LinkedOmicsKB (CPTAC pan-cancer portal)](https://kb.linkedomics.org/)
- [CPTAC pan-cancer data resource](https://pmc.ncbi.nlm.nih.gov/articles/PMC10506762)
- [GTEx proteome map, Jiang et al., Cell 2020](https://compbio.mit.edu/publications/Jiang_Cell_20.pdf)
- [Cell Surface Protein Atlas](https://pubmed.ncbi.nlm.nih.gov/25894527/)
- [SURFY surfaceome predictor](https://pmc.ncbi.nlm.nih.gov/articles/PMC6243280/)
- [ImmunoTar, Bioinformatics 2025](https://academic.oup.com/bioinformatics/article/41/3/btaf060/8008580)
- [Multi-layered proteogenomic framework for surface targets, Cancers 2026](https://doi.org/10.3390/cancers18162570)
- [otopia, bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.06.24.600275.full.pdf)
