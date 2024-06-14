for folder in /data/public/jkoubele/FLI_total_RNA/BAM_dedup_exact_position/*; do
  sample_name=$(basename "$folder")
  echo "$sample_name"
  sbatch /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/compute_coverage.sh "$sample_name" true
done