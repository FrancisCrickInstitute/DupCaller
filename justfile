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
    NEXTFLOW_PROFILE={{profile}} .venv/bin/pytest tests/mock_pipeline/test_mock_pipeline_nextflow.py -q -rs
