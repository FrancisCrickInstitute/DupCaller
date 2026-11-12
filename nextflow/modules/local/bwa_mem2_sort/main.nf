// Read through a function on purpose: it stays out of the task hash, so
// changing keep_trimmed_fastq on -resume does not rerun alignment.
def keepTrimmedFastq() { params.keep_trimmed_fastq.toString() == 'true' }

// Once its output is written, this step and MARK_DUPLICATES clear their
// input, which is another task's intermediate file: it is replaced by a
// sparse file of the same size, so it takes no disk space. The producing
// task still counts as cached, and cache 'lenient' hashes inputs by path and
// size only, so -resume still reuses this task. A cleared input starts with
// a zero byte (a FASTQ starts with '@', a BAM with 0x1f) and is refused.
process BWA_MEM2_SORT {
    tag "${sample_id}:${type}"
    cache 'lenient'
    // bwa-mem2 piped into samtools sort needs both in one image: nf-core
    // bwamem2/mem's Seqera Community image (bwa-mem2 2.3, samtools 1.22.1).
    container 'community.wave.seqera.io/library/bwa-mem2_htslib_samtools@sha256:5ebd1290d9680195817ce75915b79ae2e608834c017824b7e2bc7b141509b242'

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path fasta
    path bwa_index

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.bam"),
          path("${sample_id}_${type}.bam.bai")

    script:
    def clear_fastq = keepTrimmedFastq() ? "" : """
    for f in ${read1} ${read2}; do
        real=\$(readlink -f "\$f")
        [ -e "\$(dirname "\$real")/.command.run" ] || continue   # only another task's output
        size=\$(stat -c %s "\$real")
        truncate -s 0 "\$real" && truncate -s "\$size" "\$real" \\
            || echo "WARNING: could not clear \$real" >&2
    done
    """
    """
    set -o pipefail
    for f in ${read1} ${read2}; do
        if [ "\$(head -c 1 "\$f" | od -An -tx1 | tr -d ' ')" = "00" ]; then
            echo "ERROR: \$f was cleared after an earlier run aligned it; delete the TRIM_BARCODES work directory of ${sample_id}:${type} (or run without -resume)" >&2
            exit 1
        fi
    done

    bwa-mem2 mem -C -T 0 \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${sample_id}_${type}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${fasta} ${read1} ${read2} \\
        | samtools sort -@ ${task.cpus} > ${sample_id}_${type}.bam
    samtools index -@ ${task.cpus} ${sample_id}_${type}.bam
    ${clear_fastq}
    """

    stub:
    """
    touch ${sample_id}_${type}.bam ${sample_id}_${type}.bam.bai
    """
}
