for file in /data/public/jkoubele/FLI_total_RNA/FASTQ_trimmed/*fastq.gz; do
  file_name=$(basename "$file")
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/fastq_quality_control.sh "$file_name" true FASTQ_trimmed QC_after_trimming
done