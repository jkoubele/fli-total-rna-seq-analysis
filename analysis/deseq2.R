library("DESeq2")

cts <- as.matrix(read.delim("/cellfile/datapublic/jkoubele/FLI_total_RNA/aggregate_counts/aggregate_feature_counts.tsv", 
                            row.names="Geneid"))

# Dropping gene-level metadata added by featureCounts, keeping only the counts
cts = cts[,startsWith(colnames(cts), "no")]
class(cts) <- "numeric"
coldata <- read.delim("/cellfile/datapublic/jkoubele/FLI_total_RNA/aggregate_counts/annotation.tsv", 
                      row.names="sample_name")
coldata$group <-factor(coldata$group)

# Renaming samples, as R didn't allow the '-' character in column names of cts and replaced it with '.'
rownames(coldata) <- chartr("-", ".",rownames(coldata))
stopifnot(all(rownames(coldata) == colnames(cts)))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ group)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

contrasts_of_interest = list(c("DR", "AL"),
                             c("NAD", "AL"),
                             c("DR_NAD", "AL"),
                             c("DR_NAD", "DR"),
                             c("DR_NAD", "NAD"))

for (contrast in contrasts_of_interest){
  print(contrast)
  res <- results(dds, contrast=c("group",contrast[1],contrast[2]))
  summary(res)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), 
            file="/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/DR_vs_AL.csv")
  write.csv(as.data.frame(resOrdered), 
            file=glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/{contrast[1]}_vs_contrast[2].csv"))
}

# res <- results(dds, contrast=c("group","DR","AL"))
# 
# resLFC <- lfcShrink(dds, coef="group_DR_vs_AL", type="apeglm")
# 
# resOrdered <- res[order(res$pvalue),]
# 
# summary(res)
# plotMA(res, ylim=c(-2,2))


