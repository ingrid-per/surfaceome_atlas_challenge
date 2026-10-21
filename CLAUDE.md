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

## Environment

Conda env `surfprot`, Python 3.12 (decided 2026-10-20, see docs/decisions.md).
Activate before running anything in src/:

```bash
conda activate surfprot
```

Spec in `environment.yml`; resolved versions in `environment.lock.yml`.

## Working method

**Assert that text replacements actually matched.** A silent no-op left a
withdrawn claim live in the published deliverable after the analysis behind it
had been removed. When editing a file by string replacement, assert the pattern
was found before writing, and grep the result for the language you intended to
remove. I found it, not the tooling.

**Do not add new analysis without asking first (D60).** If something you notice
suggests a new check, metric, flag or classification, describe the observation
and what could be built, then wait for a recorded decision. Reporting what the
existing data shows is always fine; adding a new column, filter or score to the
pipeline is not something to do unilaterally. This is rule 3 of the runbook
applied to your own ideas, not just to thresholds.


The governing document is `Surface Target Pipeline Claude Code Runbook.md` in this
folder. Work one numbered step per session. Start every step in plan mode. Commit
after each step.

Thresholds listed in the runbook (20-residue extracellular minimum, FDR < 0.01,
log2FC >= 2, 5 TPM floor, vital-organ list, two-per-family limit) are the author's
proposals, NOT yet decided. They sit in the "Proposed, not yet decided" section of
docs/decisions.md and must be confirmed by me before any code uses them.

Note for Step 1: rest.uniprot.org returns 403 to requests with no User-Agent header.
Always send one.
