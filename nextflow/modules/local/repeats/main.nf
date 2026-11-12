// The STR table DupCaller.py index needs (-rt), with the README's PERF options.
process REPEATS {
    container "${params.dupcaller_container}"
    // Used only by build_reference.nf, which keeps every derived file with the genome.
    publishDir "${params.genome_base}/${params.genome}", mode: 'copy'

    input:
    path reference

    output:
    path "${reference}.repeats.tsv"

    script:
    """
    PERF -m 2 -M 10 -u 2 -i ${reference} -o ${reference}.repeats.tsv
    """

    stub:
    """
    touch ${reference}.repeats.tsv
    """
}
