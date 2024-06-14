#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=compute_coverage
#SBATCH --output=/data/public/jkoubele/cluster_logs/%j_%x.log
#SBATCH --error=/data/public/jkoubele/cluster_errors/%j_%x.err
#SBATCH --partition=all
#SBATCH --ntasks=16
#SBATCH --mem=24G

sample_name=${1:-"no003-1_OA3"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker image prune -f
  docker load -i $cell_file_prefix/jkoubele/docker_images/scientific_python.tar
fi

docker run --rm -v $cell_file_prefix/jkoubele/FLI_total_RNA/BAM_dedup_exact_position/"$sample_name":/input_folder \
  -v $cell_file_prefix/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts:/scripts \
  -v $cell_file_prefix/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38:/genome_folder \
  scientific_python /bin/sh -c ". ~/.bashrc &&\
    python /scripts/extract_unique_pairs_by_strand.py --bam_file_input /input_folder/deduplicated.bam && \
    bedtools genomecov -bga -split -i /input_folder/unique_forward.bed.gz -g genome_folder/GRCm38.primary_assembly.genome.fa.fai > /input_folder/coverage_unique_forward.bedGraph && \
    bedtools genomecov -bga -split -i /input_folder/unique_reverse.bed.gz -g genome_folder/GRCm38.primary_assembly.genome.fa.fai > /input_folder/coverage_unique_reverse.bedGraph && \
    chmod 777 -R /input_folder"
