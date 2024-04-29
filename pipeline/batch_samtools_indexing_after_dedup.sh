for file in /data/public/jkoubele/FLI_total_RNA/BAM_dedup/*/deduplicated.bam; do
  sbatch -x beyer-n03,beyer-n02 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/samtools_indexing.sh "$file" true
done