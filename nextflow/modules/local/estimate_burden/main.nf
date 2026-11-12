process ESTIMATE_BURDEN {
    tag "${sample_id}"
    container "${params.dupcaller_container}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(call_dir, stageAs: 'call_in/*')
    path fasta
    path dc_index                               // .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(gene_bed), path(gene_tbi)
    val  regions

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')
    path "${sample_id}_estimate_params.log"

    script:
    def gene_arg   = gene_bed ? "-gb ${gene_bed}" : ""
    def dilute_arg = params.estimate_dilute ? "-d" : ""
    """
    # sigProfilerPlotting caches a template pickle inside site-packages by
    # default, which fails under a container's read-only root filesystem.
    export SIGPROFILERPLOTTING_VOLUME=\$PWD/spp_templates

    # estimate writes into the call directory. The staged one is
    # CALL_VARIANTS' cached output, so work in a mirror of symlinks to it:
    # new files land in this task and the call output stays unchanged for
    # -resume. estimate appends base coverage to _stats.txt, so that file is
    # a real copy, without coverage lines from any earlier estimate.
    cp -rs "\$(readlink -f ${call_dir})" ${sample_id}
    rm ${sample_id}/${sample_id}_stats.txt
    grep -v -e '^SBS Base Coverage' -e '^Indel Base Coverage' -e '^DBS Base Coverage' \\
        "\$(readlink -f ${call_dir})/${sample_id}_stats.txt" > ${sample_id}/${sample_id}_stats.txt

    DupCaller.py estimate \\
        -i ${sample_id} \\
        -f ${fasta} \\
        -r ${regions} \\
        ${gene_arg} \\
        ${dilute_arg}
    """

    stub:
    """
    mkdir ${sample_id}
    touch ${sample_id}_estimate_params.log
    """
}
