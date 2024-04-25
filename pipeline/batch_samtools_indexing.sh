for file in /data/public/jkoubele/FLI_total_RNA/BAM/*/Aligned.sortedByCoord.out.bam; do
  sbatch -x beyer-n03 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/samtools_indexing.sh "$file" true
done