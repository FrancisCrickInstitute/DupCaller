#!/usr/bin/env nextflow

// Run once per genome: derives what DupCaller.nf --genome needs beside a genome's core
// files into <genome_base>/<genome>/, reading the FASTA where conf/genomes.config says
// it lives. A file already in that folder is not built again.
//
//   nextflow run nextflow/build_reference.nf -profile crick --genome hg38 --genome_base DIR

include { genomeSpec; indexPrefix; allExist; dupcallerIndexSuffixes; bwaMem2IndexSuffixes } from './genomes'
include { FAIDX           } from './modules/local/faidx/main'
include { REPEATS         } from './modules/local/repeats/main'
include { INDEX_REFERENCE } from './modules/local/index_reference/main'
include { BWA_MEM2_INDEX  } from './modules/local/bwa_mem2_index/main'

workflow {
    main:
    if (!params.genome) error "params.genome is required: the conf/genomes.config entry to build, e.g. --genome hg38"
    def spec   = genomeSpec()
    def fasta  = file(spec.fasta, checkIfExists: true)
    def prefix = indexPrefix(fasta)
    file("${params.genome_base}/${params.genome}").mkdirs()

    if (!file("${fasta}.fai").exists() && !allExist(prefix, ['.fai']))
        FAIDX(fasta)
    if (!allExist(prefix, dupcallerIndexSuffixes())) {
        // A repeat table already in the folder is reused, so a failed index build resumes there.
        def repeats = spec.repeat_tsv ? file(spec.repeat_tsv, checkIfExists: true)
            : allExist(prefix, ['.repeats.tsv']) ? file("${prefix}.repeats.tsv")
            : REPEATS(fasta)
        INDEX_REFERENCE(fasta, repeats)
    }
    if (!allExist(prefix, bwaMem2IndexSuffixes()))
        BWA_MEM2_INDEX(fasta)

    onComplete:
    log.info(workflow.success
        ? "${params.genome} is ready in ${params.genome_base}/${params.genome}: run DupCaller.nf with --genome ${params.genome} --genome_base ${params.genome_base}"
        : "build of ${params.genome} failed; rerun with -resume")
}
