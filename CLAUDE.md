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

The governing document is `Surface Target Pipeline Claude Code Runbook.md` in this
folder. Work one numbered step per session. Start every step in plan mode. Commit
after each step.

Thresholds listed in the runbook (20-residue extracellular minimum, FDR < 0.01,
log2FC >= 2, 5 TPM floor, vital-organ list, two-per-family limit) are the author's
proposals, NOT yet decided. They sit in the "Proposed, not yet decided" section of
docs/decisions.md and must be confirmed by me before any code uses them.

Note for Step 1: rest.uniprot.org returns 403 to requests with no User-Agent header.
Always send one.
