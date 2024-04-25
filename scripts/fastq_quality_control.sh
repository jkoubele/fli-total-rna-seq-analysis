#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=FASTQ_QC
#SBATCH --error=/data/public/jkoubele/cluster_logs/fastq_qc.log

cell_file_prefix="/cellfile/datapublic"
file=${1:-"no003-1_OA3_R1.fastq.gz"}
run_on_cluster=${2:-false}
input_folder=${3:-"20240219_866_YC"}
output_folder=${4:-"QC"}

if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/fastqc.tar
fi

docker run -v $cell_file_prefix/jkoubele/FLI_total_RNA/"$input_folder":/data_folder \
 -v $cell_file_prefix/jkoubele/FLI_total_RNA/"$output_folder":/QC biocontainers/fastqc:v0.11.9_cv8 /bin/sh -c "fastqc -t 3 -o /QC /data_folder/$file; chmod 777 -R /QC"