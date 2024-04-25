#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=deduplicate_umis
#SBATCH --error=/data/public/jkoubele/cluster_logs/deduplicate_umis.log

sample_name=${1:-"no003-1_OA3"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/umi_tools.tar
fi

docker run -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM/"$sample_name":/bam_folder \
 -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM_dedup:/BAM_dedup --security-opt seccomp=unconfined umi_tools /bin/sh -c "umi_tools dedup \
-I /bam_folder/Aligned.sortedByCoord.out.bam --paired -S /BAM_dedup/${sample_name}/deduplicated.bam; chmod 777 -R /BAM_dedup"