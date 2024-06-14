for folder in /data/public/jkoubele/FLI_total_RNA/BAM/*; do
  sample_name=$(basename "$folder")
  echo "$sample_name"
  sbatch -x beyer-n02 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/deduplicate_umi.sh "$sample_name" true
done