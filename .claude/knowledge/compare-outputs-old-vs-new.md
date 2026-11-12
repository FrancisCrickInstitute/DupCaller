# Compare outputs against the previous pipeline · type: pattern
**Rule:** after any change to `nextflow/` (wiring, config, modules, containers, Nextflow version),
run `just compare-to <old-ref> [<its NXF_VER>]` before saying the change leaves results unchanged.
It runs every mock case at both refs and diffs every published file
(`tests/mock_pipeline/compare_runs.py`); it must report 0 differences. Files only the new ref
publishes (a new case, a new output) are listed, not judged: check that list is what the change adds.
**What may differ:** only the `Run time:` line of a `*_params.log` and a PDF's `/CreationDate`.
`compare_runs.py` normalises those two and judges everything else: exact tokens, numbers to 1e-6.
**Why:** `test_mock_pipeline_nextflow.py` compares only the `default`, `metacs`, `lanes` and
`recall` cases with `expected_matched_normal/` (28 files each). The other nine are checked only for
producing files and passing the right option to DupCaller, so a change to their calls still passes.
**Blind spot:** the `udseq` case's outputs are empty (96 families, 0 effective, 0 coverage): the mock
reads carry NanoSeq's 3-base barcode, so no UDSeq pattern change can show up in any output.
**Refs:** a ref before the restructure (`64dac19` and earlier) needs Nextflow `24.10.4` and plain
`docker`; its run script adds `,local` itself. 2026-11-12: `64dac19` against `fd4c56e`, 464 files,
0 differences.
**Applies to:** every PR touching `nextflow/`, `Dockerfile` or an image pin. Asked for by the user
on the restructure (pipelines-technology-ops#157).
