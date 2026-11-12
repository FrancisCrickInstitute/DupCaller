process TRIM_BARCODES {
    tag "${sample_id}:${type}"
    container "${params.dupcaller_container}"

    input:
    // One file per lane: lane n is staged as r1_n/ and r2_n/, so equal names from
    // different lanes do not collide.
    tuple val(sample_id), val(type), path(read1, stageAs: 'r1_?/*'), path(read2, stageAs: 'r2_?/*')
    path barcode_list   // allowed barcodes for a pattern with B, else []
    val  barcode_pattern

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}_trm_1.fastq"),
          path("${sample_id}_${type}_trm_2.fastq")

    script:
    def bl_arg = barcode_list ? "-bl ${barcode_list}" : ""
    def r1 = [read1].flatten()
    def r2 = [read2].flatten()
    // trim reads its input's first bytes to detect gzip, then reopens it, so
    // several lanes are joined into a file (gzip members concatenate), deleted after.
    def join = r1.size() > 1 ? "cat ${r1.join(' ')} > lanes_1.fq\ncat ${r2.join(' ')} > lanes_2.fq" : ""
    def in1  = r1.size() > 1 ? 'lanes_1.fq' : r1[0]
    def in2  = r2.size() > 1 ? 'lanes_2.fq' : r2[0]
    """
    ${join}
    DupCaller.py trim \\
        -i  ${in1} \\
        -i2 ${in2} \\
        -p  ${barcode_pattern} \\
        ${bl_arg} \\
        -o  ${sample_id}_${type}_trm
    rm -f lanes_1.fq lanes_2.fq
    """

    stub:
    """
    touch ${sample_id}_${type}_trm_1.fastq ${sample_id}_${type}_trm_2.fastq
    """
}
