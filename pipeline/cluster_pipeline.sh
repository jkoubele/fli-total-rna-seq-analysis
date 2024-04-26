# First QC
sh ./batch_qc.sh

# Adds UMIs to read names in R1 and R2
sh ./batch_extract_umi.sh

# Extract adapters by Atria
sh ./batch_detect_adapters.sh

# Check if all adapters in R1/R2 agree and write them to separate file
sbatch -x beyer-n03 ./aggregate_adapters.sh true

# Trimming
sh ./batch_trim_adapters.sh

# QC after trimming
sh ./batch_qc_after_trimming.sh

# Alignment
sh ./batch_align.sh

# Indexing of BAM files
sh ./batch_samtools_indexing.sh

# Deduplicate UMIs
sh ./batch_deduplicate_umi.sh