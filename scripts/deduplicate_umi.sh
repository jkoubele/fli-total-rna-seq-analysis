#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=deduplicate_umis
#SBATCH --output=/data/public/jkoubele/cluster_logs/%j_%x.log
#SBATCH --error=/data/public/jkoubele/cluster_errors/%j_%x.err
#SBATCH --mem=50G
#SBATCH --partition=all
#SBATCH --ntasks=12

sample_name=${1:-"no003-1_OA3"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker image prune -f
  docker load -i $cell_file_prefix/jkoubele/docker_images/umi_tools.tar
fi

mkdir --parents $cell_file_prefix/jkoubele/FLI_total_RNA/tmp/"$sample_name"

docker run --rm -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM/"$sample_name":/bam_folder \
 -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM_dedup_with_logs:/BAM_dedup \
 -v $cell_file_prefix/jkoubele/FLI_total_RNA/tmp/"$sample_name":/tmp --security-opt seccomp=unconfined umi_tools /bin/sh -c "mkdir --parents /BAM_dedup/${sample_name}; \
umi_tools dedup -I /bam_folder/Aligned.sortedByCoord.out.bam --paired --per-gene --gene-tag GN --chimeric-pairs discard --unpaired-reads discard --temp-dir /tmp -S /BAM_dedup/${sample_name}/deduplicated.bam; chmod 777 -R /BAM_dedup; rm -r /tmp"