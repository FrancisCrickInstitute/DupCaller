#!/usr/bin/env python3
"""Compare every published file of two run_nextflow_pipeline.sh result trees.

Usage: compare_runs.py OLD/results NEW/results

Text and .npz files go through compare_outputs.compare_file (exact tokens,
numbers to 1e-6). The only content allowed to differ is the "Run time:" line
of a *_params.log and a PDF's /CreationDate. A file only NEW has is listed as
a new output, not judged; a file only OLD has, or any other difference, exits 1.
"""
import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_outputs import compare_file  # noqa: E402

# Written by test_mock_pipeline_nextflow.py next to the outputs, not by the pipeline.
IGNORED_SUFFIX = ".mapped.vcf"
PDF_DATE = re.compile(rb"/CreationDate \(D:[^)]*\)")


def differs(a, b):
    if a.read_bytes() == b.read_bytes():
        return []
    name = a.name
    if name.endswith(".log"):
        strip = lambda p: [l for l in p.read_text().splitlines() if not l.startswith("Run time:")]
        return [] if strip(a) == strip(b) else ["log differs beyond its run time"]
    if name.endswith(".pdf"):
        same = PDF_DATE.sub(b"", a.read_bytes()) == PDF_DATE.sub(b"", b.read_bytes())
        return [] if same else ["PDF differs beyond its creation date"]
    if name.endswith(".gz"):
        return [] if gzip.open(a).read() == gzip.open(b).read() else ["decompressed content differs"]
    try:
        return compare_file(a, b)
    except UnicodeDecodeError:
        return ["binary content differs"]


def main(old, new):
    files = lambda root: {
        p.relative_to(root) for p in root.rglob("*") if p.is_file() and not p.name.endswith(IGNORED_SUFFIX)
    }
    a, b = files(old), files(new)
    problems = [f"{p}: only in old" for p in sorted(a - b)]
    for rel in sorted(a & b):
        problems += [f"{rel}: {d}" for d in differs(old / rel, new / rel)[:3]]
    print(f"{len(a & b)} files compared, {len(problems)} differences, {len(b - a)} new outputs")
    for p in problems:
        print(f"  {p}")
    for p in sorted(b - a):
        print(f"  new: {p}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
