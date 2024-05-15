library("DESeq2")
library('ggplot2')
library('regionReport')
library("vsn")
library("pheatmap")
library("RColorBrewer")

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
                             c("DR_NAD", "NAD"),
                             c("DR", "NAD"))

for (contrast in contrasts_of_interest){
  print(contrast)
  res <- results(dds, contrast=c("group",contrast[1],contrast[2]))
  summary(res)
  resOrdered <- res[order(res$pvalue),]
  write.csv(as.data.frame(resOrdered), 
            file=glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/gene_DE_tables/{contrast[1]}_vs_{contrast[2]}.csv"))

  png(glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/MA_plots/{contrast[1]}_vs_{contrast[2]}.png"),
      width = 800, height = 600)
  plotMA(res, ylim = c(-3, 3))
  dev.off()
  
}

# Variance-stabilizing transformations
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

msd_ntd <- meanSdPlot(assay(ntd))
msd_vsd <- meanSdPlot(assay(vsd))
msd_rld <- meanSdPlot(assay(rld))

png(glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/transforms_mean_sd_plots/log_transform.png"),
    width = 800, height = 600)
msd_ntd$gg + ggtitle("Log transform")
dev.off()

png(glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/transforms_mean_sd_plots/vsd_transform.png"),
    width = 800, height = 600)
msd_vsd$gg + ggtitle("VSD transform")
dev.off()

png(glue::glue("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/transforms_mean_sd_plots/rld_transform.png"),
    width = 800, height = 600)
msd_rld$gg + ggtitle("Regularized Logarithm transform")
dev.off()

# Heatmaps of counts
# heatmap_indices <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# heatmap_annotation_df <- as.data.frame(colData(dds)[,"group"])
# counts_for_heatmap <- assay(ntd)[heatmap_indices,]
# rownames(heatmap_annotation_df) <-colnames(counts_for_heatmap)
# pheatmap(counts_for_heatmap, cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=heatmap_annotation_df)

sample_distance_vsd <- dist(t(assay(vsd)))
sample_distance_rld <- dist(t(assay(rld)))

sample_distance <- sample_distance_vsd

sampleDistMatrix <- as.matrix(sample_distance)
rownames(sampleDistMatrix) <- rld$group #paste(rld$group, vsd$type, sep="-")
colnames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

graphics.off()
png("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/sample_to_sample_comparison/distance_heatmap.png",
    width = 800, height = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         col=colors)
dev.off()

png("/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/sample_to_sample_comparison/PCA_rld.png",
    width = 800, height = 600)
pca_data = plotPCA(rld, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA (after Regularized Logarithm transform)") +
  geom_text(hjust=0, vjust=0, show.legend = FALSE) + scale_x_continuous(expand = c(0.2, 0.2))
dev.off()



