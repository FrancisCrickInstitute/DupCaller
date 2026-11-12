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

PIPELINE = Path(__file__).resolve().parents[2] / "nextflow" / "DupCaller.nf"
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
    hg38 = ["GRCh38.d1.vd1.fa", "af-only-gnomad.hg38.vcf.gz", "SNP_GRCh38.wgns.bed.gz", "NOISE_GRCh38.wgns.bed.gz"]
    touch(tmp_path / "genomes" / "hg38", "GRCh38.d1.vd1.fa.fai", *hg38, *(f"{hg38[0]}.{s}" for s in BWA + H5),
          *(f"{f}.tbi" for f in hg38[1:]))
    (tmp_path / "sample.map").write_text(
        "sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\nS1\tt1.fq\tt2.fq\tn1.fq\tn2.fq\n")
    (tmp_path / "tumor.map").write_text("sample_id\ttumor_fastq_1\ttumor_fastq_2\nS1\tt1.fq\tt2.fq\n")

    def run(*args):
        result = subprocess.run(
            ["nextflow", "run", str(PIPELINE), "-stub-run", "--max_cpus", "1", "--max_memory", "2.GB", *args],
            cwd=tmp_path, capture_output=True, text=True, timeout=300, env={**os.environ, **ENV})
        output = result.stdout + result.stderr
        return result.returncode, output, Counter(re.findall(r"\[PROCESS [^\]]+\] (\w+)", output))

    return run


REF = ["--sample_map", "sample.map", "--reference", "ref.fa"]
TUMOR_REF = ["--sample_map", "tumor.map", "--reference", "ref.fa"]


@pytest.mark.parametrize("args, expected", [
    (REF, FULL),
    ([*REF, "--noise_mask", "mask.bed.gz", "--gene_bed", "mask.bed.gz"], FULL),
    ([*TUMOR_REF, "--tumor_only"], TUMOR),
    ([*TUMOR_REF, "--normal_bam", "normal.bam"], TUMOR),
    ([*REF, "--skip_index", "false", "--skip_bwa_index", "false", "--repeat_tsv", "repeats.tsv"],
     {**FULL, "INDEX_REFERENCE": 1, "BWA_MEM2_INDEX": 1}),
    (["--sample_map", "sample.map", "--genome", "hg38", "--genome_base", "genomes"], FULL),
    ([*REF, "--protocol", "META-CS", "--barcode_list", "barcodes.txt"], FULL),
], ids=["matched_normal", "optional_resources", "tumor_only", "shared_normal_bam", "build_indexes", "genome", "metacs"])
def test_mode_runs_its_processes(run, args, expected):
    code, output, ran = run(*args)
    assert code == 0, output
    assert ran == expected


@pytest.mark.parametrize("args, message", [
    (["--reference", "ref.fa"], "params.sample_map is required"),
    (["--sample_map", "sample.map"], "params.reference is required"),
    ([*REF, "--genome", "hg38", "--genome_base", "genomes"], "params.genome and params.reference are both set"),
    (["--sample_map", "sample.map", "--genome", "hg38"], "params.genome needs params.genome_base"),
    (["--sample_map", "sample.map", "--genome", "hg19", "--genome_base", "genomes"], "must be one of hg38, mm39"),
    ([*REF, "--protocol", "nope"], "must be one of renanoseq"),
    ([*REF, "--protocol", "metacs"], "set barcode_list"),
    ([*REF, "--barcode_list", "barcodes.txt"], "has no B"),
    ([*REF, "--tumor_only", "--normal_bam", "normal.bam"], "tumor_only and params.normal_bam are both set"),
    # Also proves a command-line `false` arrives as a boolean, not the truthy string "false".
    ([*REF, "--skip_index", "false"], "params.repeat_tsv is required"),
    ([*REF, "--noise_mask", "plain.bed"], "must be bgzip-compressed"),
    ([*REF, "--premade", "hg38"], "premade"),
], ids=["no_sample_map", "no_reference", "genome_and_reference", "genome_without_base", "unknown_genome",
        "unknown_protocol", "metacs_without_list", "list_without_B", "tumor_only_and_normal_bam",
        "index_without_repeats", "unindexed_mask", "removed_premade_param"])
def test_bad_parameters_are_refused(run, args, message):
    code, output, ran = run(*args)
    assert code != 0, output
    assert message in output
    assert not ran
