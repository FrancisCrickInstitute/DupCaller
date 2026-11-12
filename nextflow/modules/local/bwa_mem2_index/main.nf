process BWA_MEM2_INDEX {
    // Pinned by digest, without a tag: Singularity rejects a reference with both. 2.3--he70b90d_0
    container 'quay.io/biocontainers/bwa-mem2@sha256:374f4b910b4b04d32772fccd4cab1cdcd0758356856a960cfa8b1edebfd38c9f'
    // Used only by build_reference.nf, which keeps every derived file with the genome.
    publishDir "${params.genome_base}/${params.genome}", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.{0123,amb,ann,bwt.2bit.64,pac}"

    script:
    """
    bwa-mem2 index ${reference}
    """

    stub:
    """
    touch ${reference}.{0123,amb,ann,bwt.2bit.64,pac}
    """
}
