#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=aggregate_adaptors
#SBATCH --error=/data/public/jkoubele/cluster_logs/aggregate_adaptors.log

run_on_cluster=${1:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/scientific_python.tar
fi

docker run -v $cell_file_prefix/jkoubele/FLI_total_RNA/detected_adapters:/detected_adapters \
-v $cell_file_prefix/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis:/fli-total-rna-seq-analysis \
--security-opt seccomp=unconfined \
scientific_python /bin/sh -c "python fli-total-rna-seq-analysis/extract_best_adapters.py --detected_adapters_folder /detected_adapters; chmod 777 -R /detected_adapters"
