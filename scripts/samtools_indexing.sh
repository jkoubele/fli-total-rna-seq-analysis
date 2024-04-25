#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=samtools_indexing
#SBATCH --error=/data/public/jkoubele/cluster_logs/samtools_indexing.log

bam_file_path=${1:-"/cellfile/datapublic/jkoubele/FLI_total_RNA/BAM/no004-0_OD1/Aligned.sortedByCoord.out.bam"}
run_on_cluster=${2:-false}

folder=$(dirname "$bam_file_path")
file_name=$(basename "$bam_file_path")

if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/samtools.tar
fi

docker run -v "$folder":/bam_folder staphb/samtools /bin/sh -c "samtools index /bam_folder/$file_name ;chmod 777 -R /bam_folder"


