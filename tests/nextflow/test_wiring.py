"""Wiring and parameter checks for nextflow/DupCaller.nf, run with -stub-run.

A stub run executes each process's `stub:` block on empty files: it shows
that every run mode reaches the processes it should and that a bad parameter
combination is refused. It does not show that a process's real command is
right; tests/mock_pipeline/test_mock_pipeline_nextflow.py does that, in
containers.
"""
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import pytest

NEXTFLOW_DIR = Path(__file__).resolve().parents[2] / "nextflow"
PIPELINE = NEXTFLOW_DIR / "DupCaller.nf"
BUILD = NEXTFLOW_DIR / "build_reference.nf"
# NXF_OFFLINE skips the remote nf-core include in nextflow.config.
ENV = {"NXF_ANSI_LOG": "false", "NXF_OFFLINE": "true", "NXF_VER": os.environ.get("NXF_VER", "26.04.6")}
BWA = ["0123", "amb", "ann", "bwt.2bit.64", "pac"]
H5 = ["ref.h5", "tn.h5", "hp.h5", "str.h5", "dbs.h5"]
FULL = {"TRIM_BARCODES": 2, "BWA_MEM2_SORT": 2, "MARK_DUPLICATES": 2, "CALL_VARIANTS": 1, "ESTIMATE_BURDEN": 1}
TUMOR = {**FULL, "TRIM_BARCODES": 1, "BWA_MEM2_SORT": 1, "MARK_DUPLICATES": 1}


def touch(directory, *names):
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        (directory / name).touch()


@pytest.fixture
def run(tmp_path):
    if shutil.which("nextflow") is None:
        if os.environ.get("CI"):
            pytest.fail("no nextflow on PATH under CI")
        pytest.skip("no nextflow on PATH")
    touch(tmp_path, "ref.fa", "ref.fa.fai", *(f"ref.fa.{s}" for s in BWA + H5))
    touch(tmp_path, "t1.fq", "t2.fq", "n1.fq", "n2.fq", "normal.bam", "normal.bam.bai", "repeats.tsv", "barcodes.txt")
    touch(tmp_path, "mask.bed.gz", "mask.bed.gz.tbi", "plain.bed")
    # A genome entry like conf/genomes.config's: core files in one folder, indexes built elsewhere.
    touch(tmp_path / "core", "toy.fa", "germline.vcf.gz", "germline.vcf.gz.tbi", "noise.bed.gz", "noise.bed.gz.tbi")
    core = tmp_path / "core"
    (tmp_path / "toy.config").write_text(
        f'params {{ genomes {{ toy {{ fasta = "{core}/toy.fa"; germline_vcf = "{core}/germline.vcf.gz"; '
        f'noise_mask = ["{core}/noise.bed.gz"]; regions = "chr1" }} }} }}\n')
    (tmp_path / "sample.map").write_text(
        "sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\nS1\tt1.fq\tt2.fq\tn1.fq\tn2.fq\n")
    (tmp_path / "tumor.map").write_text("sample_id\ttumor_fastq_1\ttumor_fastq_2\nS1\tt1.fq\tt2.fq\n")
    # S1's tumour over two lanes plus an empty one; its normal once, repeated on a second row.
    touch(tmp_path, "t1b.fq", "t2b.fq", "empty_1.fq", "empty_2.fq")
    for name in ("t1b.fq", "t2b.fq"):
        (tmp_path / name).write_text("@r\nA\n+\nI\n")
    (tmp_path / "t1.fq").write_text("@r\nA\n+\nI\n")
    (tmp_path / "t2.fq").write_text("@r\nA\n+\nI\n")
    (tmp_path / "n1.fq").write_text("@r\nA\n+\nI\n")
    (tmp_path / "n2.fq").write_text("@r\nA\n+\nI\n")
    (tmp_path / "lanes.map").write_text(
        "sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\n"
        "S1\tt1.fq\tt2.fq\tn1.fq\tn2.fq\n"
        "S1\tt1b.fq\tt2b.fq\tn1.fq\tn2.fq\n"
        "S1\tempty_1.fq\tempty_2.fq\t\t\n")
    touch(tmp_path / "bam", *(f"S1_{t}.mkdped.bam{x}" for t in ("tumor", "normal") for x in ("", ".bai")))
    touch(tmp_path / "bam_tumour_only", "S1_tumor.mkdped.bam", "S1_tumor.mkdped.bam.bai")

    def run(*args, script=PIPELINE):
        result = subprocess.run(
            ["nextflow", "run", str(script), "-stub-run", "-c", "toy.config",
             "--max_cpus", "1", "--max_memory", "2.GB", *args],
            cwd=tmp_path, capture_output=True, text=True, timeout=300, env={**os.environ, **ENV})
        output = result.stdout + result.stderr
        return result.returncode, output, Counter(re.findall(r"\[PROCESS [^\]]+\] (?:\w+:)*(\w+)", output))

    return run


REF = ["--sample_map", "sample.map", "--reference", "ref.fa"]
TUMOR_REF = ["--sample_map", "tumor.map", "--reference", "ref.fa"]
TOY = ["--genome", "toy", "--genome_base", "genomes"]


@pytest.mark.parametrize("args, expected", [
    (REF, FULL),
    ([*REF, "--noise_mask", "mask.bed.gz", "--gene_bed", "mask.bed.gz"], FULL),
    ([*TUMOR_REF, "--tumor_only"], TUMOR),
    ([*TUMOR_REF, "--normal_bam", "normal.bam"], TUMOR),
    # A command-line `false` must arrive as a boolean, not the truthy string "false".
    ([*TUMOR_REF, "--normal_bam", "normal.bam", "--tumor_only", "false"], TUMOR),
    ([*REF, "--protocol", "META-CS", "--barcode_list", "barcodes.txt"], FULL),
    (["--sample_map", "lanes.map", "--reference", "ref.fa"], FULL),
    ([*REF, "--from_bams", "bam"], {"CALL_VARIANTS": 1, "ESTIMATE_BURDEN": 1}),
    ([*TUMOR_REF, "--from_bams", "bam", "--tumor_only"], {"CALL_VARIANTS": 1, "ESTIMATE_BURDEN": 1}),
], ids=["matched_normal", "optional_resources", "tumor_only", "shared_normal_bam", "boolean_false", "metacs",
        "lanes", "from_bams", "from_bams_tumor_only"])
def test_mode_runs_its_processes(run, args, expected):
    code, output, ran = run(*args)
    assert code == 0, output
    assert ran == expected


@pytest.mark.parametrize("args, message", [
    (["--reference", "ref.fa"], "params.sample_map is required"),
    (["--sample_map", "sample.map"], "params.reference is required"),
    ([*REF, "--genome", "hg38", "--genome_base", "genomes"], "params.genome and params.reference are both set"),
    (["--sample_map", "sample.map", "--genome", "toy"], "params.genome needs params.genome_base"),
    (["--sample_map", "sample.map", "--genome", "hg19", "--genome_base", "genomes"], "not in params.genomes"),
    (["--sample_map", "sample.map", "--genome", "hg38", "--genome_base", "genomes"], "params.genomes.hg38.fasta is not set"),
    (["--sample_map", "sample.map", *TOY], "build_reference.nf --genome toy"),
    ([*REF, "--protocol", "nope"], "must be one of renanoseq"),
    ([*REF, "--protocol", "metacs"], "set barcode_list"),
    ([*REF, "--barcode_list", "barcodes.txt"], "has no B"),
    ([*REF, "--tumor_only", "--normal_bam", "normal.bam"], "tumor_only and params.normal_bam are both set"),
    ([*REF, "--noise_mask", "plain.bed"], "must be bgzip-compressed"),
    ([*REF, "--premade", "hg38"], "premade"),
    ([*REF, "--skip_index", "false"], "skip_index"),
    ([*REF, "--from_bams", "nowhere"], "params.from_bams = nowhere: not a folder"),
    ([*REF, "--from_bams", "bam_tumour_only"], "S1_normal.mkdped.bam"),
], ids=["no_sample_map", "no_reference", "genome_and_reference", "genome_without_base", "unknown_genome",
        "genome_without_fasta", "genome_not_built", "unknown_protocol", "metacs_without_list", "list_without_B",
        "tumor_only_and_normal_bam", "unindexed_mask", "removed_premade_param", "removed_skip_index_param",
        "from_bams_not_a_folder", "from_bams_missing_normal"])
def test_bad_parameters_are_refused(run, args, message):
    code, output, ran = run(*args)
    assert code != 0, output
    assert message in output
    assert not ran


def test_build_reference_writes_only_derived_files_then_the_genome_runs(run, tmp_path):
    core_before = sorted(p.name for p in (tmp_path / "core").iterdir())
    code, output, ran = run(*TOY, script=BUILD)
    assert code == 0, output
    assert ran == {"FAIDX": 1, "REPEATS": 1, "INDEX_REFERENCE": 1, "BWA_MEM2_INDEX": 1}
    built = sorted(p.name for p in (tmp_path / "genomes" / "toy").iterdir())
    assert built == sorted(["toy.fa.fai", "toy.fa.repeats.tsv", *(f"toy.fa.{s}" for s in BWA + H5)])
    assert sorted(p.name for p in (tmp_path / "core").iterdir()) == core_before

    code, output, ran = run(*TOY, script=BUILD)
    assert code == 0, output
    assert not ran, "a second build must find everything already built"

    code, output, ran = run("--sample_map", "sample.map", *TOY)
    assert code == 0, output
    assert ran == FULL


def test_lanes_are_joined_and_an_empty_lane_skipped(run, tmp_path):
    code, output, ran = run("--sample_map", "lanes.map", "--reference", "ref.fa")
    assert code == 0, output
    assert "skipping the tumor lane empty_1.fq" in output.replace(str(tmp_path) + "/", "")
    trims = {d.name: d for d in (tmp_path / "work").glob("*/*") if "TRIM_BARCODES" in (d / ".command.run").read_text()}
    staged = sorted(sorted(p.relative_to(d).as_posix() for p in d.glob("r1_*/*")) for d in trims.values())
    assert staged == [["r1_1/n1.fq"], ["r1_1/t1.fq", "r1_2/t1b.fq"]]


def test_publish_bams_false_on_the_command_line_publishes_nothing(run, tmp_path):
    code, output, _ = run(*REF, "--publish_bams", "false", "--outdir", "out")
    assert code == 0, output
    assert not (tmp_path / "out" / "bam").exists()
    code, output, _ = run(*REF, "--outdir", "out2")
    assert code == 0, output
    assert sorted(p.name for p in (tmp_path / "out2" / "bam").iterdir()) == sorted(
        f"S1_{t}.{x}" for t in ("tumor", "normal") for x in ("mkdped.bam", "mkdped.bam.bai", "mkdp_metrics.txt"))
