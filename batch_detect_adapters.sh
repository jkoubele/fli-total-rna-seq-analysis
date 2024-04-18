for file_read_1 in /data/public/jkoubele/FLI_total_RNA/FASTQ_with_UMI/*_R1.fastq.gz; do
  file_name=$(basename "$file_read_1")
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/detect_adapters.sh "$file_name" true
done
