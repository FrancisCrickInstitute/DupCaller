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
| Index the reference (optional) | `INDEX_REFERENCE`, `BWA_MEM2_INDEX` | `DupCaller.py index`, `bwa-mem2 index` |

Each process is a module under `modules/local/` and names its own
container. `DupCaller.nf` checks the parameters and connects the modules.

## Requirements

- Nextflow 26.04.0 or later (`manifest.nextflowVersion`; CI runs 26.04.6).
- Singularity or Docker. Every process runs in a container.
- A reference with its indexes (see [Reference genomes](#reference-genomes)).

## Running on NEMO

```bash
nextflow run nextflow/DupCaller.nf \
    -profile crick,nemo \
    --sample_map sample.map \
    --genome hg38 --genome_base /path/to/genomes \
    --protocol renanoseq \
    --outdir results -resume
```

`crick` is nf-core's institutional profile for NEMO: it sets the Slurm
executor, Singularity and the resource limits. `nemo`
(`conf/compute/nemo.config`) holds what is specific to this pipeline.
`nextflow.config` is read automatically from beside `DupCaller.nf`, so no
`-c` file is needed.

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

## Reference genomes

The pipeline never downloads a reference. It reads one of:

- `--genome hg38` or `--genome mm39` with `--genome_base DIR`: the files in
  `DIR/<genome>/`, including that genome's germline VCF and noise masks.
  A `germline_vcf` or `noise_mask` set explicitly replaces the genome's.
  For mm39, `--genome_strain STRAIN` (for example `C57BL_6NJ`) selects the
  germline VCF; without it calling runs with no germline VCF.
- `--reference ref.fa`: your own FASTA.

In both cases the FASTA has next to it its samtools index (`.fai`), its
bwa-mem2 index (`.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`) and its
DupCaller index (`.ref.h5`, `.tn.h5`, `.hp.h5`, `.str.h5`, `.dbs.h5`).
`--skip_bwa_index false` and `--skip_index false` (with `--repeat_tsv`)
build the last two inside the run instead; `bwa-mem2 index` needs about 28
bytes of RAM per reference base, about 90 GB for GRCh38.

`bin/fetch_genome.sh` fills a genome folder once, from a node with internet
access:

```bash
nextflow/bin/fetch_genome.sh hg38 /path/to/genomes
nextflow/bin/fetch_genome.sh mm39 /path/to/genomes C57BL_6NJ
```

It downloads the files listed in `assets/premade/<genome>.tsv` (GDC, GATK,
UCSC and Google Drive), checks each against the md5 in that file where one
is listed, unpacks the FASTA and runs `samtools faidx` and `bwa-mem2
index`. It needs `curl`, `md5sum`, `samtools`, `bwa-mem2` and `gdown` on
`PATH`, and about 50 GB of disk and 90 GB of RAM for hg38. A rerun skips
files that are already present. Google Drive can refuse a file for up to a
day once its download quota is reached; rerun later.

## Protocol settings

`--protocol` sets the barcode pattern and end trimming, replacing
`barcode_pattern`, `trim_template` and `trim_read`:

| Protocol | `protocol` | barcode pattern | `trim_template` | `trim_read` |
|---|---|---|---|---|
| RE-NanoSeq | `renanoseq` | `NNNXXXX` | 7 | 7 |
| NanoSeq v2 | `nanoseqv2` | `NNNXX` | 7 | 7 |
| UDSeq | `udseq` | `NNNNNNNN` | 15 | 8 |
| META-CS | `metacs` | `B` then 19 `X` (Tn5 barcode from `barcode_list`, then the 19 nt ME) | 7 | 7 |

Names are matched case-insensitively, ignoring `-` and `_`. A pattern with
`B` needs `--barcode_list`, the allowed barcodes one per line. Without
`protocol`, the three parameters are used as set (defaults `NNNXXXX`, 7
and 7). A wrong pattern breaks duplex family formation without an error,
so check it against your library.

## Parameters

`nextflow.config` lists every parameter with its default and a comment.
The five booleans (`skip_index`, `skip_bwa_index`, `tumor_only`,
`estimate_dilute`, `keep_trimmed_fastq`) are declared with a type in
`DupCaller.nf`, so `--skip_index false` on the command line is read as a
boolean. A parameter that is not declared is refused at launch, which also
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
duplicate-marked BAM, which stays in the work directory. If a step must
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
described in [docs/](../docs/index.md). BAMs stay in the work directory.

## Tests

| Command | Checks | Needs |
|---|---|---|
| `just test-wiring` | Each run mode reaches the processes it should, and bad parameter combinations are refused (`tests/nextflow/`, stub runs) | Nextflow |
| `just test-mock` | The whole pipeline on a 3.2 kb synthetic contig in 11 configurations, the default one compared with `tests/mock_pipeline/expected_matched_normal/` | Nextflow, Docker or Singularity |

A stub run executes each module's `stub:` block on empty files. It checks
the wiring, and says nothing about a process's real command; the mock test
does that. `DUPCALLER_IMAGE=IMAGE` runs the mock test with that DupCaller
image and checks that every DupCaller task ran in it. CI
(`.github/workflows/test.yml`) runs both, the mock test on an image built
from the commit.
