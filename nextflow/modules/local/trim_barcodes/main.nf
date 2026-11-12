process TRIM_BARCODES {
    tag "${sample_id}:${type}"
    container "${params.dupcaller_container}"

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path barcode_list   // allowed barcodes for a pattern with B, else []
    val  barcode_pattern

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}_trm_1.fastq"),
          path("${sample_id}_${type}_trm_2.fastq")

    script:
    def bl_arg = barcode_list ? "-bl ${barcode_list}" : ""
    """
    DupCaller.py trim \\
        -i  ${read1} \\
        -i2 ${read2} \\
        -p  ${barcode_pattern} \\
        ${bl_arg} \\
        -o  ${sample_id}_${type}_trm
    """

    stub:
    """
    touch ${sample_id}_${type}_trm_1.fastq ${sample_id}_${type}_trm_2.fastq
    """
}
