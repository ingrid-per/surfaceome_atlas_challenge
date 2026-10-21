"""Build the interactive deliverable from the template and the exported JSON.

This replaces an ad-hoc substitution done by hand. A silent no-op in a text
replacement once left a withdrawn claim live in the published page, so every
substitution here asserts that it matched, and the built file is searched for
language that must not survive.

    python src/step7_build_artifact.py [--out PATH]
"""

from __future__ import annotations

import argparse
import html
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
    "traces to a script",       # removed from the header at my request
    "fixed before the ranking ran",
    "JavaScript",               # the page must never explain itself in these terms
    "noscript",
]

# Structures the page needs. Checked here rather than discovered in the browser.
REQUIRED_IDS = ["knobs", "listR", "sbars", "sens", "rbMsg", "rbReset",
                "rbGone", "rbElig", "rbN", "rbN2", "tab-robust"]


def static_top10(payload: dict) -> str:
    """Plain HTML for the global top 10, rendered at build time.

    Visible by default and hidden by the page script on boot. A reader whose
    browser or preview window does not run scripts still sees the answer rather
    than an empty shell — which is what happened on iOS. Every value comes from
    the payload, so this cannot drift from the interactive view.
    """
    rows = []
    for i, c in enumerate(payload["top10"], 1):
        flags = []
        if c.get("existence") and c["existence"] != "Evidence at protein level":
            flags.append(f'<b>protein existence {html.escape(c["existence"].lower())}</b>')
        if c.get("stab") is not None and c["stab"] < 0.75:
            flags.append(f'survives {round(c["stab"]*100)}% of threshold settings')
        if c.get("frac10") is not None and c["frac10"] < 0.70:
            flags.append(f'present in only {round(c["frac10"]*100)}% of tumours')
        note = (f'<span class="sub">{" · ".join(flags)}</span>' if flags else "")
        above = ("&mdash;" if c.get("frac_above_healthy") is None
                 else f'{round(c["frac_above_healthy"]*100)}%')
        rows.append(
            f'<tr><td class="n">{i}</td>'
            f'<td><span class="g">{html.escape(c["symbol"])}</span>'
            f'<span class="sub">{html.escape(c["cancer"])}</span>{note}</td>'
            f'<td class="n">{c["S"]:.2f}</td>'
            f'<td class="n">{c["tumour_tpm"]:,.0f}</td>'
            f'<td>{html.escape(c["worst_tissue"])}'
            f'<span class="sub">{c["worst_tpm"]:,.1f} TPM</span></td>'
            f'<td class="n">{above}</td></tr>')
    return (
        '<div class="card"><div class="card-h">'
        '<h2>Top 10 surface targets</h2>'
        '<span class="sub">ranked by specificity score S</span></div>'
        '<div style="padding:13px 16px" class="scroll-x"><table class="stt">'
        '<tr><th>#</th><th>Target</th><th>S</th><th>Tumour TPM</th>'
        '<th>Highest healthy tissue</th><th>Tumours above it</th></tr>'
        + "".join(rows) + '</table></div></div>')


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

    n_fb = html.count("__STATIC_TOP10__")
    assert n_fb == 1, f"expected one __STATIC_TOP10__ slot, found {n_fb}"
    html = html.replace("__STATIC_TOP10__", static_top10(json.loads(data)))
    assert "__STATIC_TOP10__" not in html, "fallback placeholder survived"

    # the fallback is the only thing a script-less reader sees, so verify it
    # carries every finalist rather than trusting the renderer
    fb = html.split('<div id="static-fallback">')[1].split("<section")[0]
    missing_fb = [c["symbol"] for c in json.loads(data)["top10"]
                  if c["symbol"] not in fb]
    assert not missing_fb, f"static fallback is missing {missing_fb}"

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

    # A second, STANDALONE copy with a real document shell.
    #
    # The published artifact is wrapped in <!doctype><head><body> at publish
    # time, so the file above must NOT carry those tags. But that makes it an
    # HTML fragment: opened straight from disk or sent as an attachment it has
    # no <meta viewport>, so phones assume a ~980px screen, never fire the
    # responsive breakpoints, and render the whole page zoomed out. This copy
    # is the one to download, email or open locally.
    standalone = (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        + html.replace("<header class=\"top\">",
                       "</head>\n<body>\n<header class=\"top\">", 1)
        + "\n</body>\n</html>\n"
    )
    assert "</head>" in standalone and "<body>" in standalone, \
        "could not place the document shell — did the header markup change?"
    assert standalone.count("<!doctype html>") == 1
    sa = out.with_name(out.stem + "_standalone" + out.suffix)
    sa.write_text(standalone)

    d = json.loads(data)
    print(f"built {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}"
          f"  ({len(html)/1024:.0f} KB)")
    print(f"  tabs            : {' '.join(sorted(tabs))}")
    print(f"  robustness runs : {d['robust']['n_runs']}")
    print(f"  genes charted   : {len(d['robust']['symbols'])}")
    print(f"  standalone copy : {sa.relative_to(ROOT)}  "
          f"({len(standalone)/1024:.0f} KB, opens on a phone)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
