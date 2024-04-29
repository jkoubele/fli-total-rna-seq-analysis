for folder in /data/public/jkoubele/FLI_total_RNA/BAM_dedup/*; do
  sample_name=$(basename "$folder")
  echo "$sample_name"
  sbatch -x beyer-n01,beyer-n02,beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/count_features.sh "$sample_name" true
done