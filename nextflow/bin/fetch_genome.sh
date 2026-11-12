#!/usr/bin/env bash
# Downloads a pre-made genome into GENOME_BASE/<genome>/, the folder the
# pipeline reads with --genome <genome> --genome_base GENOME_BASE. Run it
# once, on a node with internet access; a rerun skips what is already there.
#
# Usage: fetch_genome.sh hg38|mm39 GENOME_BASE [STRAIN]
#   STRAIN (mm39 only, e.g. C57BL_6NJ) also fetches that strain's germline VCF.
#
# Needs curl, md5sum, samtools and bwa-mem2 on PATH, and gdown for the files
# hosted on Google Drive. bwa-mem2 index needs ~90 GB of RAM for hg38.
set -euo pipefail

genome="${1:?usage: fetch_genome.sh hg38|mm39 GENOME_BASE [STRAIN]}"
base="${2:?usage: fetch_genome.sh hg38|mm39 GENOME_BASE [STRAIN]}"
strain="${3:-}"
manifest="$(cd "$(dirname "${BASH_SOURCE[0]}")/../assets/premade" && pwd)/${genome}.tsv"
[ -f "$manifest" ] || { echo "no manifest for ${genome}: ${manifest}" >&2; exit 1; }

mkdir -p "${base}/${genome}"
cd "${base}/${genome}"

# Columns: file, url, md5 (empty where the source publishes none).
tail -n +2 "$manifest" | while IFS=$'\t' read -r file url md5; do
    # mm39 lists one germline VCF per mouse strain.
    case "$file" in mgp.*) [[ "$file" == "mgp.${strain}_reh.vcf.gz"* ]] || continue ;; esac
    [ -e "$file" ] && continue
    echo "fetching ${file}" >&2
    case "$url" in
        *drive.google.com*) gdown --no-cookies -O "${file}.part" "$url" ;;
        *)                  curl -fL -o "${file}.part" "$url" ;;
    esac
    if [ -n "$md5" ]; then echo "${md5}  ${file}.part" | md5sum -c - >&2; fi
    mv "${file}.part" "$file"
done

# The first manifest row is the FASTA archive.
archive="$(sed -n 2p "$manifest" | cut -f1)"
fasta="${archive%.gz}"; fasta="${fasta%.tar}"
if [ ! -e "$fasta" ]; then
    case "$archive" in
        *.tar.gz) tar -xzf "$archive" "$fasta" ;;
        *)        gzip -dc "$archive" > "${fasta}.part" && mv "${fasta}.part" "$fasta" ;;
    esac
fi
[ -e "${fasta}.fai" ]  || samtools faidx "$fasta"
[ -e "${fasta}.0123" ] || bwa-mem2 index "$fasta"
echo "${genome} is ready in $(pwd)" >&2
