// Shared by DupCaller.nf and build_reference.nf.

def genomeSpec() {
    def spec = params.genomes?.get(params.genome)
    if (!spec)
        error "params.genome = ${params.genome}: not in params.genomes (conf/genomes.config); known: ${params.genomes?.keySet()?.join(', ')}"
    if (!spec.fasta)
        error "params.genomes.${params.genome}.fasta is not set: add the FASTA's path to conf/genomes.config"
    if (!params.genome_base)
        error "params.genome needs params.genome_base, the folder build_reference.nf writes ${params.genome}'s indexes into"
    spec
}

// Where build_reference.nf writes a genome's derived files, named after its FASTA.
def indexPrefix(fasta) { "${params.genome_base}/${params.genome}/${fasta.name}" }

def allExist(prefix, suffixes) { suffixes.every { file("${prefix}${it}").exists() } }

def dupcallerIndexSuffixes() { ['ref', 'tn', 'hp', 'str', 'dbs'].collect { ".${it}.h5" } }
def bwaMem2IndexSuffixes() { ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'].collect { ".${it}" } }
