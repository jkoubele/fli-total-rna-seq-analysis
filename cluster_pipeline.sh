# First QC
sh ./batch_qc.sh

# Adds UMIs to read names in R1 and R2
sh ./batch_extract_umi.sh

# Extract adapters by Atria
sh ./batch_detect_adapters.sh

# Check if all adapters in R1/R2 agree and write them to separate file
sbatch -x beyer-n03 ./aggregate_adapters.sh true
