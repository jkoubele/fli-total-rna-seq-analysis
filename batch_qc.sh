for file in /cellfile/datapublic/jkoubele/FLI_total_RNA/20240219_866_YC/*fastq.gz; do
  file_name=$(basename "$file")
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/fastq_quality_control.sh $file_name true
done