"""Build the interactive deliverable from the template and the exported JSON.

This replaces an ad-hoc substitution done by hand. A silent no-op in a text
replacement once left a withdrawn claim live in the published page, so every
substitution here asserts that it matched, and the built file is searched for
language that must not survive.

    python src/step7_build_artifact.py [--out PATH]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "docs" / "step7_atlas_template.html"
PAYLOAD = ROOT / "results" / "step7_targets.json"
DEFAULT_OUT = ROOT / "results" / "step7_atlas.html"

# Text that must NOT appear in the built page. Each entry is something that was
# removed deliberately; if it reappears, a template edit has regressed.
FORBIDDEN = [
    "infiltrate flag",          # withdrawn with D61
    "38% of genes",             # withdrawn with D61
    "__DATA__",                 # the placeholder itself
]

# Structures the page needs. Checked here rather than discovered in the browser.
REQUIRED_IDS = ["knobs", "listR", "sbars", "sens", "rbMsg", "rbReset",
                "rbGone", "rbElig", "rbN", "rbN2", "tab-robust"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    tpl = TEMPLATE.read_text()
    data = PAYLOAD.read_text()
    json.loads(data)  # fail here rather than in the browser

    n = tpl.count("__DATA__")
    assert n == 1, f"expected exactly one __DATA__ placeholder, found {n}"
    html = tpl.replace("__DATA__", data)
    assert "__DATA__" not in html, "placeholder survived substitution"

    for phrase in FORBIDDEN:
        assert phrase not in html, f"forbidden text present in build: {phrase!r}"

    missing = [i for i in REQUIRED_IDS if f'id="{i}"' not in html]
    assert not missing, f"template is missing required elements: {missing}"

    # every tab button must have a matching section, and the switcher must know it
    tabs = set(re.findall(r'data-tab="([a-z]+)"', html))
    for t in tabs:
        assert f'id="tab-{t}"' in html, f"tab button '{t}' has no section"
    switch = re.search(r'\[([^\]]+)\]\.forEach\(t=>\s*\n?\s*document\.getElementById\("tab-"', html)
    assert switch, "could not find the tab switcher list"
    known = set(re.findall(r'"([a-z]+)"', switch.group(1)))
    assert tabs == known, f"switcher {sorted(known)} != buttons {sorted(tabs)}"

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)

    d = json.loads(data)
    print(f"built {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}"
          f"  ({len(html)/1024:.0f} KB)")
    print(f"  tabs            : {' '.join(sorted(tabs))}")
    print(f"  robustness runs : {d['robust']['n_runs']}")
    print(f"  genes charted   : {len(d['robust']['symbols'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
