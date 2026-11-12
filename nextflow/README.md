# DupCaller Nextflow pipeline

Runs DupCaller end to end for one or more samples: barcode trimming,
bwa-mem2 alignment, GATK duplicate marking, `DupCaller.py call` and
`DupCaller.py estimate`. The commands and tool versions are those of the
main README's [Pipeline](../README.md#pipeline) section.

| Step | Process | Tool |
|---|---|---|
| Trim barcodes | `TRIM_BARCODES` | `DupCaller.py trim` |
| Align and sort | `BWA_MEM2_SORT` | bwa-mem2 2.3 (`mem -C -T 0`), samtools 1.22.1 |
| Mark duplicates | `MARK_DUPLICATES` | GATK 4.3.0.0 MarkDuplicates |
| Call mutations | `CALL_VARIANTS` | `DupCaller.py call` |
| Estimate burden | `ESTIMATE_BURDEN` | `DupCaller.py estimate` |

Each process is a module under `modules/local/` and names its own
container. `DupCaller.nf` checks the parameters and connects the modules.
`build_reference.nf` is a separate workflow, run once per genome, that
builds the indexes `DupCaller.nf` reads (see
[Reference genomes](#reference-genomes)).

## Requirements

- Nextflow 26.04.0 or later (`manifest.nextflowVersion`; CI runs 26.04.6).
- Singularity or Docker. Every process runs in a container.
- A genome built once with `build_reference.nf`, or your own reference with its indexes (see [Reference genomes](#reference-genomes)).

## Running on NEMO

```bash
nextflow run nextflow/DupCaller.nf \
    -profile crick \
    --sample_map sample.map \
    --genome hg38 --genome_base /path/to/genomes \
    --protocol renanoseq \
    --outdir results -resume
```

`crick` is nf-core's institutional profile for NEMO
([`conf/crick.config`](https://github.com/nf-core/configs/blob/master/conf/crick.config)),
fetched at launch, so the launching node needs internet access: it sets
the Slurm executor, Singularity with `/nemo` and `/flask` bound, and limits
of 2 TB, 256 CPUs and 168 h, which replace `max_cpus`, `max_memory` and
`max_time`. `nextflow.config` is read automatically from beside
`DupCaller.nf`, so no `-c` file is needed.

Without a matched normal, add `--tumor_only` and leave the normal columns
out of the sample map. `DupCaller.py call` then runs without `-n` and with
`-maf 0.1` unless `max_af` is set below 1, so variants above 10% allele
fraction, which are mostly germline, are not called. `-maf 0.1` also masks
positions with raw tumour depth below 10x.

To use one aligned normal for every sample, pass `--normal_bam
normal.bam` (with `normal.bam.bai` next to it). Trimming, alignment and
duplicate marking are then skipped for the normal. `tumor_only` and
`normal_bam` cannot be combined.

## Sample map

Tab-separated, with a header. `examples/sample.map` and
`examples/sample_tumor_only.map` are templates.

```
sample_id	tumor_fastq_1	tumor_fastq_2	normal_fastq_1	normal_fastq_2
```

The normal columns are read only when each sample has its own normal to
align, that is when neither `tumor_only` nor `normal_bam` is set.

A sample sequenced over several lanes takes one row per lane, and its
lanes are joined, in row order, before trimming. A row adds a tumour lane,
and a normal lane when its normal columns are filled; leave them empty on
the other rows, or repeat the same normal pair, which is used once:

```
sample_id	tumor_fastq_1	tumor_fastq_2	normal_fastq_1	normal_fastq_2
S1	S1_L001_R1.fastq.gz	S1_L001_R2.fastq.gz	N1_R1.fastq.gz	N1_R2.fastq.gz
S1	S1_L002_R1.fastq.gz	S1_L002_R2.fastq.gz
```

A lane with an empty file is skipped with a warning, as the lab's
`merge_fastqs.sh` did. A sample's lanes must be all gzipped or all plain,
because they are joined by concatenation.

## Calling again from BAMs

Each duplicate-marked BAM, its index and its MarkDuplicates metrics are
copied to `${outdir}/bam/` (`--publish_bams false` turns this off). To call
again with other settings, without trimming and aligning again, point
`--from_bams` at that folder:

```bash
nextflow run nextflow/DupCaller.nf -profile crick \
    --sample_map sample.map --genome hg38 --genome_base /path/to/dupcaller_genomes \
    --from_bams results/bam --max_af 0.2 --outdir results_maf0.2
```

Only the sample map's `sample_id` column is read. Each sample's
`<sample_id>_tumor.mkdped.bam`, and `<sample_id>_normal.mkdped.bam` unless
`tumor_only` or `normal_bam` is set, must be in the folder. This replaces
the lab's `call_variants.nf`.

## Reference genomes

A genome's core files stay where they are at the Crick and are never
copied. `conf/genomes.config`, laid out like nf-core's `igenomes.config`,
gives each genome's FASTA, germline VCF, noise masks and calling regions.
What DupCaller adds to a genome is built once, by `build_reference.nf`,
into a folder of its own:

```bash
nextflow run nextflow/build_reference.nf -profile crick \
    --genome hg38 --genome_base /path/to/dupcaller_genomes
```

| Built into `<genome_base>/<genome>/` | Step | Notes |
|---|---|---|
| `<fasta>.repeats.tsv` | `REPEATS`: `PERF -m 2 -M 10 -u 2` | Skipped when the genome entry gives `repeat_tsv` |
| `<fasta>.{ref,tn,hp,str,dbs}.h5` | `INDEX_REFERENCE`: `DupCaller.py index` | The DupCaller index; about 18 GB for GRCh38 |
| `<fasta>.{0123,amb,ann,bwt.2bit.64,pac}` | `BWA_MEM2_INDEX`: `bwa-mem2 index` | About 28 bytes of RAM per base, about 90 GB for GRCh38 |
| `<fasta>.fai` | `FAIDX` | Only when the FASTA has no `.fai` beside it |

A file already in the folder is not built again, so a rerun after a
failure picks up where it stopped. Every run then reads the genome with:

```bash
nextflow run nextflow/DupCaller.nf -profile crick \
    --genome hg38 --genome_base /path/to/dupcaller_genomes ...
```

The FASTA, germline VCF and masks are read from their Crick locations and
the indexes from the folder; Nextflow stages them side by side, which is
where DupCaller and bwa-mem2 look. A `germline_vcf`, `noise_mask` or
`regions` given on the command line replaces the genome's. A run whose
genome has not been built stops at launch and prints the
`build_reference.nf` command.

Another genome, or the same genome from other files, is one more block in
`conf/genomes.config` or in a `-c` file:

```groovy
params {
    genomes {
        hg38_noalt {
            fasta        = '/path/to/GRCh38_no_alt.fa'
            germline_vcf = '/path/to/af-only-gnomad.hg38.vcf.gz'
            noise_mask   = ['/path/to/NOISE_GRCh38.wgns.bed.gz']
            regions      = 'chr1 chr2 ... chrX'
        }
    }
}
```

An index built by another DupCaller version may not be read: DupCaller
1.2.12 refuses one without its current index version, which an index
from 1.1.x lacks. Build the folder with the DupCaller image the pipeline
runs.

Without `--genome`, `--reference ref.fa` reads your own FASTA with its
`.fai`, bwa-mem2 index and DupCaller index all beside it.

## Protocol settings

`--protocol` sets the barcode pattern and end trimming, replacing
`barcode_pattern`, `trim_template` and `trim_read`:

| Protocol | `protocol` | barcode pattern | `trim_template` | `trim_read` |
|---|---|---|---|---|
| RE-NanoSeq | `renanoseq` | `NNNXXXX` | 7 | 7 |
| NanoSeq v2 | `nanoseqv2` | `NNNXX` | 7 | 7 |
| UDSeq | `udseq` | `NNNNNNNNXXXXXXX` | 15 | 8 |
| META-CS | `metacs` | `B` then 19 `X` (Tn5 barcode from `barcode_list`, then the 19 nt ME) | 7 | 7 |

The UDSeq pattern is the Crick lab's (an 8-base barcode, then 7 skipped
bases), not upstream's `NNNNNNNN`; the main README's
[Protocol Settings](../README.md#protocol-settings) table is upstream's.
The lab's earlier UDSeq runs used this pattern with end trimming of 7 and 7;
UDSeq libraries need 15 and 8, which `udseq` sets.
Names are matched case-insensitively, ignoring `-` and `_`. A pattern with
`B` needs `--barcode_list`, the allowed barcodes one per line. Without
`protocol`, the three parameters are used as set (defaults `NNNXXXX`, 7
and 7). A wrong pattern breaks duplex family formation without an error,
so check it against your library.

## Parameters

`nextflow.config` lists every parameter with its default and a comment.
The three booleans (`tumor_only`, `estimate_dilute`,
`keep_trimmed_fastq`) are declared with a type in `DupCaller.nf`, so
`--tumor_only false` on the command line is read as a boolean. A parameter that is not declared is refused at launch, which also
catches a misspelt name.

Every optional VCF or BED (`germline_vcf`, `noise_mask`, `target_bed`,
`indel_bed`, `gene_bed`) must be bgzip-compressed with its tabix index
(`.gz.tbi`) next to it. `noise_mask` takes one path, or a list in a `-c` or
`-params-file` file.

## Resources

`conf/base.config` sizes each process for whole-genome duplex libraries.
`threads` (default 16) is the CPU count for `BWA_MEM2_SORT` and
`CALL_VARIANTS`.

| Process | CPUs | Memory | Time |
|---|---|---|---|
| `TRIM_BARCODES` | 1 | 4 GB | 12 h |
| `BWA_MEM2_SORT` | `threads` | 32 GB + 1 GB per CPU | 24 h |
| `MARK_DUPLICATES` | 1 | 64 GB | 24 h |
| `CALL_VARIANTS` | `threads` | 6 GB per CPU | 96 h |
| `ESTIMATE_BURDEN` | 1 | 8 GB | 4 h |

A task killed for memory or time (exit 137, 143, 104, 134 or 139) is
retried twice with the request scaled by the attempt number.
`max_cpus`, `max_memory` and `max_time` cap every request; the `crick`
profile sets them for NEMO.

Other profiles: `docker`, `singularity`, and `arm` (with `docker`, on Apple
Silicon, where every image runs under amd64 emulation).

## Disk and `-resume`

Intermediate files are cleared once the next step has written its output:
`BWA_MEM2_SORT` clears the trimmed FASTQs and `MARK_DUPLICATES` clears the
sorted BAM. A cleared file is a sparse file of the original size, so it
takes no disk space and `-resume` still reuses the finished tasks. Disk use
per library peaks during alignment (trimmed FASTQs, about five times the
input fq.gz, plus the sorted BAM) and afterwards is mostly the
duplicate-marked BAM, which stays in the work directory and, unless
`--publish_bams false`, a copy of it in `${outdir}/bam/`. If a step must
rerun after its input was cleared, it stops with an error naming the
upstream task whose work directory to delete.

`--keep_trimmed_fastq` keeps the trimmed FASTQs. Changing it on `-resume`
does not rerun alignment.

## Containers

| Process | Image |
|---|---|
| `TRIM_BARCODES`, `CALL_VARIANTS`, `ESTIMATE_BURDEN`, `INDEX_REFERENCE` | `params.dupcaller_container`, default `ghcr.io/alexandrovlab/dupcaller:1.2.12` |
| `BWA_MEM2_SORT` | `community.wave.seqera.io/library/bwa-mem2_htslib_samtools` (bwa-mem2 2.3, samtools 1.22.1) |
| `BWA_MEM2_INDEX` | `quay.io/biocontainers/bwa-mem2` 2.3 |
| `MARK_DUPLICATES` | `quay.io/biocontainers/gatk4` 4.3.0.0 |

The third-party images are pinned by digest in their modules.
`--dupcaller_container` runs another DupCaller image, for example a local
`.sif` or one built from a branch.

## Output

Each sample's results are copied to `${outdir}/${sample_id}/`, the
directory `DupCaller.py call -o` and `estimate -i` produce on their own
(VCFs, coverage BEDs, error-profile tables, burden and signature outputs),
plus `${outdir}/${sample_id}_estimate_params.log`. The outputs are
described in [docs/](../docs/index.md). The duplicate-marked BAMs and their
metrics go to `${outdir}/bam/` ([Calling again from BAMs](#calling-again-from-bams)).

hg38's calling regions are chr1 to chr22, chrX and chrY, as the lab calls
them; upstream DupCaller's default leaves out chrY.

## Tests

| Command | Checks | Needs |
|---|---|---|
| `just test-wiring` | Each run mode reaches the processes it should, and bad parameter combinations are refused (`tests/nextflow/`, stub runs) | Nextflow |
| `just test-mock` | The whole pipeline on a 3.2 kb synthetic contig in 13 configurations; the default, `metacs`, `lanes` and `recall` cases compared with `tests/mock_pipeline/expected_matched_normal/` | Nextflow, Docker or Singularity |
| `just compare-to <ref>` | Every published file of every case, at `<ref>` and here | As `just test-mock` |

A stub run executes each module's `stub:` block on empty files. It checks
the wiring, and says nothing about a process's real command; the mock test
does that. `DUPCALLER_IMAGE=IMAGE` runs the mock test with that DupCaller
image and checks that every DupCaller task ran in it. CI
(`.github/workflows/test.yml`) runs both, the mock test on an image built
from the commit.
