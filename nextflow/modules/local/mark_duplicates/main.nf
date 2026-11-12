process MARK_DUPLICATES {
    tag "${sample_id}:${type}"
    cache 'lenient'   // see BWA_MEM2_SORT
    // 4.3.0.0--py36hdfd78af_0
    container 'quay.io/biocontainers/gatk4@sha256:49f91bc0fcb704271c18a519715bf161141c0a00be4cf49c7ef2866de1f59a46'
    // Kept for calling again later without realigning (DupCaller.nf --from_bams).
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.publish_bams.toString() != 'false'

    input:
    tuple val(sample_id), val(type), path(bam), path(bai)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.mkdped.bam"),
          path("${sample_id}_${type}.mkdped.bam.bai"), emit: bam
    path "${sample_id}_${type}.mkdp_metrics.txt", emit: metrics

    script:
    // Sort files spill to the task work dir, not Java's default tmp dir
    // (often a small node-local /tmp).
    def heap_gb = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    """
    if [ "\$(head -c 1 ${bam} | od -An -tx1 | tr -d ' ')" = "00" ]; then
        echo "ERROR: ${bam} was cleared after an earlier run marked duplicates on it; delete the BWA_MEM2_SORT work directory of ${sample_id}:${type} (or run without -resume)" >&2
        exit 1
    fi

    gatk --java-options "-Xmx${heap_gb}g -Djava.io.tmpdir=\$PWD" MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}_${type}.mkdped.bam \\
        -M ${sample_id}_${type}.mkdp_metrics.txt \\
        --TMP_DIR \$PWD \\
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \\
        --DUPLEX_UMI \\
        --TAGGING_POLICY OpticalOnly \\
        --BARCODE_TAG DB \\
        --CREATE_INDEX true

    # The GATK image has no samtools; MarkDuplicates writes X.mkdped.bai.
    mv ${sample_id}_${type}.mkdped.bai ${sample_id}_${type}.mkdped.bam.bai

    for f in ${bam} ${bai}; do
        real=\$(readlink -f "\$f")
        [ -e "\$(dirname "\$real")/.command.run" ] || continue   # only another task's output
        size=\$(stat -c %s "\$real")
        truncate -s 0 "\$real" && truncate -s "\$size" "\$real" \\
            || echo "WARNING: could not clear \$real" >&2
    done
    """

    stub:
    """
    touch ${sample_id}_${type}.mkdped.bam ${sample_id}_${type}.mkdped.bam.bai ${sample_id}_${type}.mkdp_metrics.txt
    """
}
