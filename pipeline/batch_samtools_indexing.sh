for file in /data/public/jkoubele/FLI_total_RNA/BAM/*/Aligned.sortedByCoord.out.bam; do
  sbatch -x beyer-n05,beyer-n01,beyer-n03,beyer-n02 /data/public/jkoubele/FLI_total_RNA/fli-total-rna-seq-analysis/scripts/samtools_indexing.sh "$file" true
done