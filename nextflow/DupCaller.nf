#!/usr/bin/env nextflow

include { genomeSpec; indexPrefix; dupcallerIndexSuffixes; bwaMem2IndexSuffixes } from './genomes'
include { TRIM_BARCODES   } from './modules/local/trim_barcodes/main'
include { BWA_MEM2_SORT   } from './modules/local/bwa_mem2_sort/main'
include { MARK_DUPLICATES } from './modules/local/mark_duplicates/main'
include { CALL_VARIANTS   } from './modules/local/call_variants/main'
include { ESTIMATE_BURDEN } from './modules/local/estimate_burden/main'

// Typed so that a command-line value converts: from Nextflow 26.04 an untyped
// `--tumor_only false` binds the String "false", which is truthy. Every other
// parameter is declared in nextflow.config.
params {
    tumor_only: Boolean         = false
    estimate_dilute: Boolean    = false
    keep_trimmed_fastq: Boolean = false
}

// Library protocols: barcode pattern (trim -p) and end trimming (call
// -tt/-tr), as in the main README's Protocol Settings table, except UDSeq's
// pattern: the Crick lab's libraries carry 7 skipped bases after the 8-base
// barcode (dupcaller-crick's preset), where upstream's preset has none.
def protocols() {
    [
        renanoseq: [barcode_pattern: 'NNNXXXX',      trim_template: 7,  trim_read: 7],
        nanoseqv2: [barcode_pattern: 'NNNXX',        trim_template: 7,  trim_read: 7],
        udseq:     [barcode_pattern: 'NNNNNNNNXXXXXXX', trim_template: 15, trim_read: 8],
        // Tn5 transposon barcode from barcode_list, then the 19 nt ME.
        metacs:    [barcode_pattern: 'B' + 'X' * 19, trim_template: 7,  trim_read: 7],
    ]
}


// Checks the parameters, then resolves those that genome, protocol and
// tumor_only decide. Processes take the result as inputs, never by calling a
// function: -resume's task hash does not see a function's return value.
def resolveSettings() {
    if (!params.sample_map) error "params.sample_map is required"

    def genome = null
    if (params.genome) {
        if (params.reference)
            error "params.genome and params.reference are both set: genome supplies the reference, so leave reference unset"
        genome = genomeSpec()
    } else if (!params.reference) {
        error "params.reference is required (or set params.genome, e.g. hg38)"
    }

    def protocol = null
    if (params.protocol) {
        // Matched case-insensitively, ignoring '-' and '_' (e.g. RE-NanoSeq).
        protocol = protocols()[params.protocol.toString().toLowerCase().replaceAll(/[-_]/, '')]
        if (!protocol)
            error "params.protocol = ${params.protocol}: must be one of ${protocols().keySet().join(', ')}"
    }
    def barcode_pattern = protocol ? protocol.barcode_pattern : params.barcode_pattern
    if (barcode_pattern.contains('B') && !params.barcode_list)
        error "barcode pattern ${barcode_pattern} has B (one barcode from a list): set barcode_list to the file of allowed barcodes, one per line"
    if (params.barcode_list && !barcode_pattern.contains('B'))
        error "params.barcode_list is set but barcode pattern ${barcode_pattern} has no B, so the list would not be used"

    if (params.from_bams && !file(params.from_bams).isDirectory())
        error "params.from_bams = ${params.from_bams}: not a folder; give the bam/ folder an earlier run published"
    if (params.tumor_only && params.normal_bam)
        error "params.tumor_only and params.normal_bam are both set: tumor_only calls without a matched normal, so leave normal_bam unset"
    // Without a matched normal an uncapped max_af would call germline
    // variants, so tumor_only caps it at 0.1 unless it is set below 1.
    def max_af = (params.tumor_only && (params.max_af as double) >= 1) ? 0.1 : params.max_af
    if (max_af != params.max_af)
        log.info "tumor_only with max_af = ${params.max_af}: DupCaller.py call runs with -maf ${max_af}, so variants above that allele fraction (mostly germline) are not called."

    def masks = params.noise_mask ?: genome?.noise_mask ?: []

    [
        reference:       genome ? genome.fasta : params.reference,
        germline_vcf:    params.germline_vcf ?: genome?.germline_vcf,
        noise_mask:      masks instanceof List ? masks : [masks],
        // Without a genome entry, hg38's autosomes, X and Y, as the lab calls them.
        regions:         params.regions ?: genome?.regions ?: (1..22).collect { "chr${it}" }.join(' ') + ' chrX chrY',
        barcode_pattern: barcode_pattern,
        trim_template:   protocol ? protocol.trim_template : params.trim_template,
        trim_read:       protocol ? protocol.trim_read     : params.trim_read,
        max_af:          max_af,
    ]
}

// DupCaller opens every VCF/BED with tabix, so each is staged with its
// .tbi. An unset one is [[], []], which the process reads as absent.
def indexed(p, param_name) {
    if (!p) return [[], []]
    if (!p.toString().endsWith('.gz'))
        error "params.${param_name} = ${p}: must be bgzip-compressed (.gz) with a tabix index (.gz.tbi) next to it"
    [file(p, checkIfExists: true), file("${p}.tbi", checkIfExists: true)]
}

// A sample's reads may span several rows, one per lane: each row adds a tumour
// lane, and a normal lane when its normal_fastq_1/2 are filled. A lane with an
// empty file is skipped with a warning, as are repeats of a lane already listed.
def laneRows(row, types) {
    types.collect { type ->
        def fq1 = row["${type}_fastq_1"]
        def fq2 = row["${type}_fastq_2"]
        if (!fq1 && !fq2) return null
        if (!fq1 || !fq2) error "sample ${row.sample_id}: ${type}_fastq_1 and ${type}_fastq_2 must both be set or both empty"
        def lane = [file(fq1, checkIfExists: true), file(fq2, checkIfExists: true)]
        if (lane.any { it.size() == 0 }) {
            log.warn "sample ${row.sample_id}: skipping the ${type} lane ${fq1}, ${fq2}: a file is empty"
            return null
        }
        [row.sample_id, type] + lane
    }.findAll()
}

workflow alignReads {
    take:
    rows_ch
    types
    fasta
    bwa_index
    barcode_pattern

    main:
    def lanes_ch = rows_ch
        .flatMap { row -> laneRows(row, types) }
        .unique()
        .map { sid, type, fq1, fq2 -> [[sid, type], fq1, fq2] }
        .groupTuple()
        .map { key, fq1s, fq2s -> [key[0], key[1], fq1s, fq2s] }
    def barcode_list = params.barcode_list ? file(params.barcode_list, checkIfExists: true) : []
    def trimmed_ch = TRIM_BARCODES(lanes_ch, barcode_list, barcode_pattern)
    def aligned_ch = BWA_MEM2_SORT(trimmed_ch, fasta, bwa_index)

    emit:
    MARK_DUPLICATES(aligned_ch).bam
}

workflow {
    def s     = resolveSettings()
    def fasta = file(s.reference, checkIfExists: true)

    // With --genome the indexes are in <genome_base>/<genome>/, named after the FASTA
    // (build_reference.nf), and the .fai is the FASTA's own if it has one. With
    // --reference everything sits beside the FASTA. Staged together, so DupCaller and
    // bwa-mem2 find them beside the FASTA's name either way.
    def prefix  = params.genome ? indexPrefix(fasta) : s.reference
    def fai     = file("${s.reference}.fai").exists() ? "${s.reference}.fai" : "${prefix}.fai"
    def wanted  = [fai] + (dupcallerIndexSuffixes() + bwaMem2IndexSuffixes()).collect { "${prefix}${it}" }
    def missing = wanted.findAll { !file(it).exists() }
    if (missing)
        error "missing index files for ${fasta.name}: ${missing.join(', ')}. " +
              (params.genome
                  ? "Build them once with: nextflow run nextflow/build_reference.nf --genome ${params.genome} --genome_base ${params.genome_base}"
                  : "Build them beside the FASTA (samtools faidx, DupCaller.py index, bwa-mem2 index), or use --genome")
    def bwa_index = bwaMem2IndexSuffixes().collect { file("${prefix}${it}") }
    def dc_index  = [file(fai)] + dupcallerIndexSuffixes().collect { file("${prefix}${it}") }

    // Each sample's own normal is aligned unless normal_bam (one normal for every
    // sample) or tumor_only is set.
    def align_normal = !params.normal_bam && !params.tumor_only
    def types = align_normal ? ['tumor', 'normal'] : ['tumor']
    def rows_ch = channel.fromPath(params.sample_map, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)

    // (sample_id, type, BAM, BAI) after duplicate marking: aligned here, or the
    // BAMs an earlier run published to <outdir>/bam/ (from_bams), to call again.
    def markdup_ch = params.from_bams
        ? rows_ch.map { row -> row.sample_id }.unique().flatMap { sid ->
              types.collect { type ->
                  def bam = "${params.from_bams}/${sid}_${type}.mkdped.bam"
                  [sid, type, file(bam, checkIfExists: true), file("${bam}.bai", checkIfExists: true)]
              }
          }
        : alignReads(rows_ch, types, fasta, bwa_index, s.barcode_pattern)

    def tumor_ch = markdup_ch
        .filter { _sid, type, _bam, _bai -> type == 'tumor' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }
    def normal_ch = markdup_ch
        .filter { _sid, type, _bam, _bai -> type == 'normal' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }
    // normal_bam is one pre-aligned, indexed normal shared by every sample.
    def shared_normal = params.normal_bam
        ? [file(params.normal_bam, checkIfExists: true), file("${params.normal_bam}.bai", checkIfExists: true)]
        : [[], []]
    def call_input_ch = align_normal
        ? tumor_ch.join(normal_ch)
        : tumor_ch.map { sid, bam, bai -> [sid, bam, bai] + shared_normal }

    def calls_ch = CALL_VARIANTS(
        call_input_ch,
        fasta,
        dc_index,
        indexed(s.germline_vcf, 'germline_vcf'),
        s.noise_mask.collectMany { indexed(it, 'noise_mask') },
        s.noise_mask.collect { file(it).name },
        indexed(params.target_bed, 'target_bed'),
        indexed(params.indel_bed, 'indel_bed'),
        [regions: s.regions, max_af: s.max_af, trim_template: s.trim_template, trim_read: s.trim_read],
    )

    ESTIMATE_BURDEN(calls_ch, fasta, dc_index, indexed(params.gene_bed, 'gene_bed'), s.regions)
}
