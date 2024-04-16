#!/bin/bash

#SBATCH --account=jkoubele
#SBATCH --job-name=extract_umis
#SBATCH --error=/data/public/jkoubele/cluster_logs/extract_umis.log

file_name_read_1=${1:-"no003-1_OA3_R1.fastq.gz"}
run_on_cluster=${2:-false}

cell_file_prefix="/cellfile/datapublic"
if $run_on_cluster; then
  cell_file_prefix="/data/public"
  docker load -i $cell_file_prefix/jkoubele/docker_images/umi_tools.tar
fi

file_name_stem="${file_name_read_1%%.*}"
trimmed_stem="${file_name_stem::-2}"

file_name_read_2="${trimmed_stem}R2.fastq.gz"
file_name_umi="${trimmed_stem}UMI.fastq.gz"

#file_read_1=$cell_file_prefix/jkoubele/FLI_total_RNA/20240219_866_YC/${file_name_read_1}
#file_read_2="${cell_file_prefix}/jkoubele/FLI_total_RNA/20240219_866_YC/${file_name_read_2}"
#file_umi="${cell_file_prefix}/jkoubele/FLI_total_RNA/20240219_866_YC/${trimmed_stem}UMI.fastq.gz"

#${cell_file_prefix}/jkoubele/anaconda3/bin/umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=$file_umi --read2-in=$file_read_1 \
#--stdout=${cell_file_prefix}/jkoubele/FLI_total_RNA/FASTQ_with_UMI/"$file_name_read_1"  --read2-stdout
#${cell_file_prefix}/jkoubele/anaconda3/bin/umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=$file_umi --read2-in="$file_read_2" \
#--stdout=${cell_file_prefix}/jkoubele/FLI_total_RNA/FASTQ_with_UMI/"$file_name_read_2" --read2-stdout

docker run -v $cell_file_prefix/jkoubele/FLI_total_RNA/20240219_866_YC:/data_folder \
 -v $cell_file_prefix/jkoubele/FLI_total_RNA/FASTQ_with_UMI:/FASTQ_with_UMI --security-opt seccomp=unconfined umi_tools /bin/sh -c "umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=/data_folder/$file_name_umi --read2-in=/data_folder/$file_name_read_1 \
--stdout=/FASTQ_with_UMI/$file_name_read_1  --read2-stdout; umi_tools extract --bc-pattern=NNNNNNNNNNN --stdin=/data_folder/$file_name_umi --read2-in=/data_folder/$file_name_read_2 \
--stdout=/FASTQ_with_UMI/$file_name_read_2  --read2-stdout; chmod 777 -R /FASTQ_with_UMI"