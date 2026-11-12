process CALL_VARIANTS {
    tag "${sample_id}"
    container "${params.dupcaller_container}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(tumor_bam),  path(tumor_bai),
          path(normal_bam), path(normal_bai)    // normal is [] with tumor_only
    path fasta
    path dc_index                               // .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(germline_vcf), path(germline_tbi)
    path noise_mask_files                       // mask .bed.gz + .tbi files
    val  noise_mask_names                       // basenames of the masks only, in -m order
    tuple path(target_bed), path(target_tbi)
    tuple path(indel_bed),  path(indel_tbi)
    val  opts                                   // [regions, max_af, trim_template, trim_read]

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')

    script:
    def normal_arg   = normal_bam   ? "-n ${normal_bam}"   : ""
    def germline_arg = germline_vcf ? "-g ${germline_vcf}" : ""
    def noise_arg    = noise_mask_names ? "-m " + noise_mask_names.collect { "'${it}'" }.join(' ') : ""
    def target_arg   = target_bed   ? "-R ${target_bed}"   : ""
    def indel_arg    = indel_bed    ? "-id ${indel_bed}"   : ""
    // Unset, DupCaller picks and logs a random seed, and uses its own -pt default (0.05).
    def seed_arg     = params.seed != null        ? "--seed ${params.seed}"      : ""
    def pt_arg       = params.p_threshold != null ? "-pt ${params.p_threshold}" : ""
    """
    DupCaller.py call \\
        -b  ${tumor_bam} \\
        ${normal_arg} \\
        -f  ${fasta} \\
        -o  ${sample_id} \\
        -p  ${task.cpus} \\
        -r  ${opts.regions} \\
        ${germline_arg} \\
        ${noise_arg} \\
        ${target_arg} \\
        ${indel_arg} \\
        ${seed_arg} \\
        ${pt_arg} \\
        -maf ${opts.max_af} \\
        -gaf ${params.germline_af_cutoff} \\
        -d   ${params.min_n_depth} \\
        -tt  ${opts.trim_template} \\
        -tr  ${opts.trim_read} \\
        -mq  ${params.mapq}
    """

    stub:
    """
    mkdir ${sample_id}
    touch ${sample_id}/${sample_id}_stats.txt
    """
}
