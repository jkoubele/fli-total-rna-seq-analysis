#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=count_features_htseq
#SBATCH --output=/data/public/jkoubele/cluster_logs/%j_%x.log
#SBATCH --error=/data/public/jkoubele/cluster_errors/%j_%x.err
#SBATCH --mem=16G
#SBATCH --partition=all
#SBATCH --ntasks=6

sample_name=${1:-"no003-1_OA3"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/htseq.tar
fi

docker run --rm -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM_dedup/"$sample_name":/bam_folder \
-v $cell_file_prefix/jkoubele/FLI_total_RNA/gene_counts_htseq:/gene_counts \
-v $cell_file_prefix/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38:/GRCm38 \
--security-opt seccomp=unconfined htseq /bin/sh -c "htseq-count -f bam -r pos /bam_folder/deduplicated.bam /GRCm38/gencode.vM10.primary_assembly.annotation.gtf > /gene_counts/${sample_name}.txt; chmod 777 -R /gene_counts"