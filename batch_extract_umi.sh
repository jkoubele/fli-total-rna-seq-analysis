for file in /data/public/jkoubele/FLI_total_RNA/20240219_866_YC/*_R1.fastq.gz; do
  file_name=$(basename "$file")
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/extract_umi.sh "$file_name" true
done