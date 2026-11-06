"""Step 1: download the pinned datasets and record exactly what arrived.

Defaults to a dry run. Nothing is fetched unless --download is passed, because
the rules say you see the file names and sizes before anything is downloaded.

Usage:
    python src/step1_download.py              # dry run: probe and print the table
    python src/step1_download.py --download   # actually fetch
    python src/step1_download.py --download --only X1,U2
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import stat
import sys
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step1_resources import EXCLUDED, MAGIC, RESOURCES, USER_AGENT, Resource

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
TIMEOUT = 120
CHUNK = 1 << 20

# An HTML page served at a .xlsx path is the failure mode this script exists to
# catch, so these are the bytes we refuse outright.
HTML_STARTS = (b"<!do", b"<!DO", b"<htm", b"<HTM", b"<?xm")


def human(n: int | None) -> str:
    if n is None:
        return "?"
    return f"{n / 1048576:.2f} MB" if n >= 1048576 else f"{n / 1024:.0f} KB"


def _request(url: str, method: str = "GET"):
    return urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT}, method=method
    )


def probe(res: Resource) -> tuple[int | None, str, str]:
    """HEAD the URL. Returns (bytes or None, status, content-type)."""
    try:
        with urllib.request.urlopen(_request(res.url, "HEAD"), timeout=60) as r:
            cl = r.headers.get("Content-Length")
            return (int(cl) if cl else None, str(r.status),
                    (r.headers.get("Content-Type") or "").split(";")[0])
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code}", ""
    except Exception as e:  # noqa: BLE001
        return None, type(e).__name__, ""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def validate(path: Path, res: Resource) -> str | None:
    """Return an error string if the file is not the format it claims to be."""
    with path.open("rb") as fh:
        head = fh.read(512)
    if not head:
        return "file is empty"
    if head[:4] in HTML_STARTS or head.lstrip()[:5].lower() == b"<html":
        return ("server returned an HTML page, not data — the URL has probably "
                "moved (HTTP 200 does not mean you got the file)")
    want = MAGIC.get(res.magic)
    if want and not head.startswith(want):
        return f"expected {res.magic} magic {want!r}, got {head[:4]!r}"
    if res.magic == "text" and b"\t" not in head and b"," not in head:
        return "plain-text file has no tab or comma in its first 512 bytes"
    return None


def unwritable(path: Path) -> None:
    """Strip write bits. data/raw is read-only, per the rules in docs/CLAUDE.md."""
    mode = path.stat().st_mode
    path.chmod(mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def download(res: Resource, dest: Path) -> tuple[int, str]:
    """Fetch to a .part file, validate, then rename. Returns (bytes, sha256)."""
    part = dest.with_suffix(dest.suffix + ".part")
    if part.exists():
        part.chmod(0o644)
        part.unlink()
    total = 0
    with urllib.request.urlopen(_request(res.url), timeout=TIMEOUT) as r, \
            part.open("wb") as out:
        while chunk := r.read(CHUNK):
            out.write(chunk)
            total += len(chunk)
            if total % (64 << 20) < CHUNK:
                print(f"      ... {human(total)}", flush=True)

    err = validate(part, res)
    if err:
        part.unlink()
        raise RuntimeError(err)

    if dest.exists():
        dest.chmod(0o644)
        dest.unlink()
    part.rename(dest)
    unwritable(dest)
    return total, sha256(dest)


FIELDS = ["id", "resource", "release", "filename", "url", "bytes",
          "sha256", "date_retrieved", "note"]


def write_manifest(rows: list[dict]) -> None:
    """Merge these rows into data/MANIFEST.md and the CSV mirror.

    Merging matters: a targeted run such as `--only U3` processes one resource,
    and a plain overwrite here would silently delete the other 40 manifest rows
    while leaving their files on disk.
    """
    (ROOT / "results").mkdir(exist_ok=True)
    csv_path = ROOT / "results" / "step1_manifest.csv"

    merged: dict[str, dict] = {}
    if csv_path.exists():
        with csv_path.open() as fh:
            for r in csv.DictReader(fh):
                r["bytes"] = int(r["bytes"])
                merged[r["id"]] = r
    for r in rows:
        merged[r["id"]] = r

    # Keep the declared order so the manifest reads the same way every time.
    order = {res.rid: i for i, res in enumerate(RESOURCES)}
    rows = sorted(merged.values(), key=lambda r: order.get(r["id"], 10_000))

    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    md = ROOT / "data" / "MANIFEST.md"
    text = md.read_text()
    table = ["| File | Resource | Version / release | Date retrieved | URL | SHA-256 | Size |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        table.append(
            f"| `{r['filename']}` | {r['resource']} | {r['release']} | "
            f"{r['date_retrieved']} | <{r['url']}> | `{r['sha256']}` | "
            f"{human(r['bytes'])} |")
    table.append("")
    table.append(f"**{len(rows)} files, {human(sum(r['bytes'] for r in rows))} total.**")
    table.append("")
    table.append("### Deliberately not downloaded")
    table.append("")
    table.append("| Resource | Why not |")
    table.append("| --- | --- |")
    for name, why in EXCLUDED:
        table.append(f"| {name} | {why} |")

    start = text.index("## Files")
    end = text.index("## Expected resources")
    md.chmod(0o644)
    md.write_text(text[:start] + "## Files\n\n" + "\n".join(table) + "\n\n" + text[end:])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true",
                    help="actually fetch (default is a dry run)")
    ap.add_argument("--only", default="", help="comma-separated resource ids")
    args = ap.parse_args()

    todo = RESOURCES
    if args.only:
        keep = {s.strip() for s in args.only.split(",")}
        todo = [r for r in RESOURCES if r.rid in keep]

    RAW.mkdir(parents=True, exist_ok=True)
    if RAW.exists():
        RAW.chmod(0o755)  # need write access to add files; locked again at the end

    if not args.download:
        print(f"DRY RUN — nothing will be downloaded. {len(todo)} files.\n")
        print(f"{'id':10} {'probed':>11}  {'pinned':>11}  {'status':8} file")
        print("-" * 100)
        total = 0
        for r in todo:
            got, status, _ = probe(r)
            total += got or 0
            flag = ""
            if r.expected_bytes and got and abs(got - r.expected_bytes) / r.expected_bytes > 0.05:
                flag = "  <-- SIZE DIFFERS FROM PINNED"
            if status != "200":
                flag = f"  <-- {status}"
            print(f"{r.rid:10} {human(got):>11}  {human(r.expected_bytes):>11}  "
                  f"{status:8} {r.filename}{flag}")
        print("-" * 100)
        print(f"TOTAL: {human(total)}  across {len(todo)} files")
        print(f"\nRun with --download to fetch into {RAW}")
        return 0

    rows, failures = [], []
    for i, r in enumerate(todo, 1):
        dest = RAW / r.filename
        print(f"[{i}/{len(todo)}] {r.rid}  {r.filename}")
        try:
            if dest.exists():
                n, digest = dest.stat().st_size, sha256(dest)
                print(f"      already present, {human(n)}, reusing")
            else:
                n, digest = download(r, dest)
                print(f"      done, {human(n)}")
            if r.expected_bytes and abs(n - r.expected_bytes) / r.expected_bytes > 0.05:
                print(f"      WARNING: pinned size was {human(r.expected_bytes)}")
            rows.append({"id": r.rid, "resource": r.group, "release": r.release,
                         "filename": r.filename, "url": r.url, "bytes": n,
                         "sha256": digest, "date_retrieved": date.today().isoformat(),
                         "note": r.note})
        except Exception as e:  # noqa: BLE001
            print(f"      FAILED: {type(e).__name__}: {e}")
            failures.append((r.rid, r.filename, str(e)))

    if rows:
        write_manifest(rows)
        print(f"\nManifest updated: {len(rows)} files, "
              f"{human(sum(r['bytes'] for r in rows))}")

    # Lock data/raw. The rule is that raw data is never edited or deleted.
    for p in RAW.iterdir():
        if p.is_file():
            unwritable(p)
    RAW.chmod(0o555)
    print(f"{RAW} set read-only")

    if failures:
        print(f"\n{len(failures)} FAILED:")
        for rid, fn, err in failures:
            print(f"  - {rid} {fn}: {err}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
