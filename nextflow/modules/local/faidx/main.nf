process FAIDX {
    container "${params.dupcaller_container}"
    // Used only by build_reference.nf, which keeps every derived file with the genome.
    publishDir "${params.genome_base}/${params.genome}", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.fai"

    script:
    """
    python -c "import pysam; pysam.faidx('${reference}')"
    """

    stub:
    """
    touch ${reference}.fai
    """
}
