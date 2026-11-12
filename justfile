set shell := ["bash", "-eu", "-o", "pipefail", "-c"]

# The pipeline needs Nextflow >=26.04; an unpinned older launcher reports a syntax error instead.
export NXF_VER := "26.04.6"

default: test

# .venv with DupCaller installed editable, plus pytest.
dev:
    [ -d .venv ] || uv venv .venv
    uv pip install -p .venv/bin/python -e ".[test]"

# The fast layers: the caller's unit tests, then the pipeline's wiring.
test: test-unit test-wiring

# Needs tabix on PATH.
test-unit:
    .venv/bin/pytest tests/unit -q

# Stub runs: no container, no data.
test-wiring:
    .venv/bin/pytest tests/nextflow -q

# The whole pipeline in containers, 11 cases. On Apple Silicon: just test-mock docker,arm
test-mock profile="":
    PATH="$PWD/.venv/bin:$PATH" NEXTFLOW_PROFILE={{profile}} .venv/bin/pytest tests/mock_pipeline/test_mock_pipeline_nextflow.py -q -rs

# Run all 11 mock cases at REF and at the working tree, then diff every published file.
# REF_NXF is the Nextflow version REF supports (before the restructure: 24.10.4).
compare-to REF REF_NXF="26.04.6" profile="docker":
    #!/usr/bin/env bash
    set -euo pipefail
    out=$(mktemp -d)
    echo "runs in $out"
    git worktree add -q --detach "$out/ref" {{REF}}
    trap 'git worktree remove --force "$out/ref"' EXIT
    PATH="$PWD/.venv/bin:$PATH" NXF_VER={{REF_NXF}} bash "$out/ref/tests/mock_pipeline/run_nextflow_pipeline.sh" "$out/old" {{profile}} all
    PATH="$PWD/.venv/bin:$PATH" bash tests/mock_pipeline/run_nextflow_pipeline.sh "$out/new" {{profile}} all
    .venv/bin/python tests/mock_pipeline/compare_runs.py "$out/old/results" "$out/new/results"
