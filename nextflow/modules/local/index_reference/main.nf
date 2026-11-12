process INDEX_REFERENCE {
    container "${params.dupcaller_container}"

    input:
    path reference
    path repeat_tsv

    output:
    path "${reference}.{ref,tn,hp,str,dbs}.h5"

    script:
    """
    DupCaller.py index -f ${reference} -rt ${repeat_tsv}
    """

    stub:
    """
    touch ${reference}.{ref,tn,hp,str,dbs}.h5
    """
}
