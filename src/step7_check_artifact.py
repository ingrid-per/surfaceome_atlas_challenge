"""Checks 7.10-7.12: execute the built deliverable's JavaScript and assert it
renders what it claims to.

The page cannot be opened in a browser from here, and a published page that
throws on load looks identical to a working one until someone clicks on it. So
the page script is run under JavaScriptCore against a minimal DOM
(src/artifact_domstub.js) and its output is inspected.

This exists because a silent failure has already happened once on this project:
a string replacement matched nothing and a withdrawn claim stayed live.

Requires jsc, which ships with macOS. Skips with a clear message elsewhere.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAGE = ROOT / "results" / "step7_atlas.html"
STUB = ROOT / "src" / "artifact_domstub.js"
JSC = Path("/System/Library/Frameworks/JavaScriptCore.framework/Versions/A"
           "/Helpers/jsc")

ASSERTIONS = r"""
var hadFail = false;
function ok(cond, msg){
  if (!cond){ print("  FAIL  " + msg); hadFail = true; }
  else print("  pass  " + msg);
}

/* --- the robustness explorer --- */
var knobs = document.getElementById("knobs").innerHTML;
ok(knobs.split('class="knob"').length - 1 === D.robust.axes.length,
   "one control per threshold axis");

var listR = document.getElementById("listR").innerHTML;
var nDefault = D.robust.runs[D.robust.default_key].length;
ok(listR.split('class="row"').length - 1 === nDefault,
   "default run lists all " + nDefault + " genes");
ok(listR.indexOf(D.top10[0].symbol) >= 0,
   "the top-ranked finalist appears in the default run");
ok(document.getElementById("rbGone").cls.indexOf("hidden") >= 0,
   "nothing drops out at the decided settings");
ok(document.getElementById("sbars").innerHTML.split('class="sbar').length - 1 === 25,
   "25 stability bars drawn");
ok(document.getElementById("sens").innerHTML.indexOf("spread") >= 0,
   "sensitivity panel rendered");

/* every reachable setting must resolve to a stored run */
var axes = D.robust.axes, combos = [[]];
axes.forEach(function(a){
  var next = [];
  combos.forEach(function(c){ a.values.forEach(function(v){ next.push(c.concat([v])); }); });
  combos = next;
});
var missing = combos.filter(function(c){ return !D.robust.runs[c.join("|")]; });
ok(missing.length === 0, combos.length + " selectable settings all resolve to a run");

/* moving a threshold must actually change something */
var segs = document.querySelectorAll(".seg");
ok(segs.length === axes.reduce(function(n,a){ return n + a.values.length; }, 0),
   "every axis value has a control");
var alt = segs.filter(function(s){
  return s.dataset.ax === "summary" && s.dataset.v === "p90"; })[0];
ok(!!alt, "the 90th-percentile control exists");
alt.click();
ok(document.getElementById("listR").innerHTML !== listR,
   "moving a threshold changes the list");
ok(document.getElementById("rbGone").cls.indexOf("hidden") < 0,
   "genes that drop out are named");
ok(document.getElementById("rbReset").disabled === false,
   "reset becomes available once a threshold moves");
document.getElementById("rbReset").onclick();
ok(document.getElementById("listR").innerHTML === listR,
   "reset restores the decided list exactly");

/* --- the rest of the page still works --- */
ok(document.getElementById("detailB").innerHTML.indexOf("Survives our thresholds") >= 0,
   "each finalist carries its stability");
ok(document.getElementById("listA").innerHTML.indexOf('class="row"') >= 0,
   "the by-cancer-type tab rendered");
ok(document.getElementById("listC").innerHTML.indexOf('class="row"') >= 0,
   "the pan-cancer tab rendered");
ok(document.getElementById("funnel").innerHTML.indexOf("<tr>") >= 0,
   "the funnel rendered");
var method = document.getElementById("method").innerHTML;
ok(method.indexOf("${") < 0, "no unresolved template expression in the method text");
ok(method.indexOf("undefined") < 0, "no 'undefined' leaked into the method text");

/* tab buttons must not have had their handlers overwritten by each other */
var tabBtns = document.querySelectorAll("nav.tabs button");
ok(tabBtns.length === 5, "five tabs wired");
ok(segs[0].onclick !== tabBtns[0].onclick,
   "threshold controls and tab buttons have separate handlers");

print(hadFail ? "\nARTIFACT CHECKS FAILED" : "\nall artifact checks passed");
if (hadFail) throw new Error("artifact checks failed");
"""


def main() -> int:
    if not JSC.exists():
        print("jsc not available on this platform — skipping artifact checks")
        return 0
    if not PAGE.exists():
        print(f"{PAGE} not found — run src/step7_build_artifact.py first")
        return 1

    html = PAGE.read_text()
    body = html.split("<script>")[0]
    m = re.search(r"<script>(.*)</script>", html, re.S)
    assert m, "no <script> block in the built page"

    harness = (
        "var STATIC_HTML = " + json.dumps(body) + ";\n"
        + STUB.read_text()
        + "\nparseInto(STATIC_HTML, null);\n"
        + m.group(1)
        + "\nprint('  pass  page script ran to completion');\n"
        + ASSERTIONS
    )
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as fh:
        fh.write(harness)
        path = fh.name

    r = subprocess.run([str(JSC), path], capture_output=True, text=True)
    print(r.stdout.rstrip())
    if r.stderr.strip():
        print(r.stderr.rstrip(), file=sys.stderr)
    Path(path).unlink(missing_ok=True)
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
