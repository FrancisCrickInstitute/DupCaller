# Compare outputs against the previous pipeline · type: pattern
**Rule:** after any change to `nextflow/` (wiring, config, modules, containers, Nextflow version),
run the old and the new pipeline on the mock data for **all 11 cases** and diff **every published
file**, before saying the change leaves results unchanged.
**How:** check the old commit out in a worktree and run
`tests/mock_pipeline/run_nextflow_pipeline.sh OUTDIR docker all` there, under the Nextflow version
it supports (`NXF_VER`). Run the new pipeline the same way, with `docker,arm` on Apple Silicon. Then
compare `OUTDIR/results/<case>/` file by file, using `compare_file` from
`tests/mock_pipeline/compare_outputs.py` (exact tokens, numbers to 1e-6). Logs carrying paths or
timings may differ; anything else that differs is a finding.
**Why:** `test_mock_pipeline_nextflow.py` compares only the `default` and `metacs` cases with
`expected_matched_normal/` (28 files). The other nine cases are checked only for producing files
and passing the right option to DupCaller, so a change that alters their calls still passes.
**Applies to:** every PR touching `nextflow/`, `Dockerfile` or an image pin. Asked for by the user
on the restructure (pipelines-technology-ops#157), 2026-11-12.
