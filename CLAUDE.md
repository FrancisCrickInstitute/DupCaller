# DupCaller (Crick fork) — repository guidance

Repo-specific traps and patterns are in `.claude/knowledge/`: read `INDEX.md` before changing the
pipeline.

A fork of [AlexandrovLab/DupCaller](https://github.com/AlexandrovLab/DupCaller). `main` tracks
upstream; `crick` is the base branch here, so branch from it and open pull requests against it.

The Python package in `src/` does the scientific work and belongs to upstream. Change it only for
a fix we intend to send upstream, and keep such a change separate from pipeline work so that
merging `main` into `crick` stays simple. The Crick-specific part is the Nextflow pipeline in
`nextflow/`, which upstream ships as a single script and this fork has split into modules.

## Layout

| Path | Holds |
|---|---|
| `src/DupCaller.py`, `src/DupCaller_sub/` | the caller: `trim`, `index`, `call`, `estimate`, `learn`, `summarize` |
| `nextflow/DupCaller.nf` | parameter checks and the workflow; no process is defined here |
| `nextflow/modules/local/<name>/main.nf` | one process each, with its `container` and a `stub:` block |
| `nextflow/nextflow.config` | every parameter default, the profiles and the manifest |
| `nextflow/conf/base.config` | per-process CPUs, memory and time |
| `nextflow/conf/compute/nemo.config` | NEMO settings on top of nf-core's `crick` profile |
| `nextflow/bin/fetch_genome.sh` | one-off download of hg38 or mm39 into a genome folder |
| `tests/unit/`, `tests/mock_pipeline/` | upstream's tests of the caller and of the whole pipeline |
| `tests/nextflow/` | stub-run wiring and parameter tests |

## Commands

```bash
just dev                    # .venv with the package editable and pytest
just test                   # unit tests, then the wiring tests
just test-wiring            # stub runs only: about a minute, no containers
just test-mock [profile]    # 11 real runs in containers; `docker,arm` on Apple Silicon
just compare-to <ref> [nxf] # all 11 cases at <ref> and here, every published file diffed
```

Run the pipeline as `NXF_VER=26.04.6 nextflow run nextflow/DupCaller.nf -profile crick,nemo ...`.
`nextflow/README.md` is the user documentation.

## Rules for the pipeline

- **Nextflow 26.04 or later only.** The entry script uses a typed `params` block, which 25.10 does
  not parse. A bare `nextflow run` on an older launcher reports a syntax error that names no
  version, so export `NXF_VER=26.04.6`.
- **A new boolean parameter goes in the typed `params` block in `DupCaller.nf`.** From 26.04 an
  untyped `--flag false` binds the string `"false"`, which is truthy. Every other parameter is
  declared in `nextflow.config`; an undeclared one is refused at launch.
- **An absent optional file is `[]`,** or `[[], []]` for a file with its `.tbi`. A module tests the
  input's truth to decide whether to pass the option.
- **A value that depends on other parameters is resolved in `resolveSettings()` and reaches a
  process as an input.** A process that calls a function for it gets a task hash that does not see
  the value, so `-resume` reuses a stale task. `keepTrimmedFastq()` in `bwa_mem2_sort` uses that
  effect deliberately.
- **`BWA_MEM2_SORT` and `MARK_DUPLICATES` clear their inputs** (sparse files of the same size) and
  set `cache 'lenient'`. Keep both together when editing either module.
- **Every process has a `stub:` block and a `script:` block.** After adding a process or a run
  mode, add a row to `tests/nextflow/test_wiring.py`.
- **A stub run checks wiring only, and the mock test compares only two cases with expected
  outputs.** Before calling a pipeline change output-neutral, run `just compare-to <old-ref>`,
  which diffs every published file of all 11 cases (`.claude/knowledge/compare-outputs-old-vs-new.md`).
- **The pipeline does not download references.** `--genome` reads `<genome_base>/<genome>/`, and
  `genomes()` in `DupCaller.nf` must name the files `fetch_genome.sh` leaves there.
- **The DupCaller version is written three times:** `__version__` in `src/DupCaller_sub/__init__.py`,
  and the `dupcaller_container` tag and `manifest.version` in `nextflow/nextflow.config`.
  `release.yml` fails the release when they disagree.

## Not yet settled

- `params.genome_base` has no default yet. hg38 and mm39 are already on NEMO; pointing
  `conf/compute/nemo.config` at them, after checking their file names match `genomes()`, is
  pipelines-technology-ops#162.
- Work items live in `FrancisCrickInstitute/pipelines-technology-ops`, because issues are disabled
  on this fork: #157 for this restructure, #158 for deployment and the Airflow integration.
- The DupCaller image is upstream's `ghcr.io/alexandrovlab/dupcaller`, which is correct while
  `crick` carries no change to `src/`. `release.yml` still publishes to that registry, which this
  fork cannot write to. The first Python change on `crick` needs a Crick image and that workflow
  repointed.
- The pipeline is not under [proteus](https://github.com/FrancisCrickInstitute/proteus), unlike
  chandra, castor and altair. The stub blocks and `tests/nextflow/` are the starting point for it.
