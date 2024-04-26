#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=detect_adaptors
#SBATCH --output=/data/public/jkoubele/cluster_logs/%j_%x.log
#SBATCH --error=/data/public/jkoubele/cluster_errors/%j_%x.err
#SBATCH --partition=all
#SBATCH --ntasks=3

file_name_read_1=${1:-"no003-1_OA3_R1.fastq.gz"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/atria.tar
fi

file_name_stem="${file_name_read_1%%.*}"
trimmed_stem="${file_name_stem::-3}"
file_name_read_2="${trimmed_stem}_R2.fastq.gz"

docker run --rm -v $cell_file_prefix/jkoubele/FLI_total_RNA/FASTQ_with_UMI:/FASTQ_with_UMI \
-v $cell_file_prefix/jkoubele/FLI_total_RNA/detected_adapters:/detected_adapters \
--security-opt seccomp=unconfined \
atria /bin/sh -c "atria-4.0.3/bin/atria --read1 /FASTQ_with_UMI/${file_name_read_1} --read2 /FASTQ_with_UMI/${file_name_read_2} --detect-adapter;
cp atria_adapter_detect_summary* /detected_adapters/${trimmed_stem}.txt; chmod 777 -R /detected_adapters"
