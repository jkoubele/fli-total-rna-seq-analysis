#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=align
#SBATCH --error=/data/public/jkoubele/cluster_logs/align.log

file_name_read_1=${1:-"no003-1_OA3_R1.atria.fastq.gz"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
fi

file_name_stem="${file_name_read_1%%.*}"
trimmed_stem="${file_name_stem::-3}"
file_name_read_2="${trimmed_stem}_R2.atria.fastq.gz"

mkdir --parents $cell_file_prefix/jkoubele/FLI_total_RNA/BAM/"${trimmed_stem}"

$cell_file_prefix/jkoubele/STAR_2.7.11a/Linux_x86_64_static/STAR \
--runThreadN 16 \
--genomeDir $cell_file_prefix/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38/STAR_generated_genome \
--readFilesIn $cell_file_prefix/jkoubele/FLI_total_RNA/FASTQ_trimmed/"${file_name_read_2}" $cell_file_prefix/jkoubele/FLI_total_RNA/FASTQ_trimmed/"${file_name_read_1}" \
--readFilesCommand zcat \
--outFileNamePrefix $cell_file_prefix/jkoubele/FLI_total_RNA/BAM/"${trimmed_stem}"/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes GX GN \
--genomeLoad LoadAndKeep \
--limitBAMsortRAM 30000000000
