for file in /data/public/jkoubele/FLI_total_RNA/FASTQ_trimmed/*R1.atria.fastq.gz; do
  file_name=$(basename "$file")
  echo "$file_name"
  sbatch /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/align.sh "$file_name" true
done