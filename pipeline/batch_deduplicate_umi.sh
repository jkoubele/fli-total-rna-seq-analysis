for folder in /data/public/jkoubele/FLI_total_RNA/BAM/*; do
  sample_name=$(basename "$folder")
  if [ $sample_name = "no003-1_OA3" ]; then
   continue
  fi
  echo "$sample_name"
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/deduplicate_umi.sh "$sample_name" true
done