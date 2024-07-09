library(stats)
library(GenomicRanges)
library(rtracklayer)
library(rjson)

parent_folder <- "/cellfile/datapublic/jkoubele/FLI_total_RNA/BAM_dedup_exact_position/"
sj_selected <- read.table('/cellfile/datapublic/jkoubele/FLI_total_RNA/SJ_annotation/sj_selected.tsv')
output_folder <- "/cellfile/datapublic/jkoubele/FLI_total_RNA/intronic_slopes/"

padding <- 20

for(sample_name in list.files(parent_folder)){
  print(sample_name)
  folder_path <- paste0("/cellfile/datapublic/jkoubele/FLI_total_RNA/BAM_dedup_exact_position/", sample_name, '/')
  
  
  bed_graph_forward <- import.bedGraph(paste0(folder_path, 'coverage_unique_forward.bedGraph'))
  bed_graph_reverse <- import.bedGraph(paste0(folder_path, 'coverage_unique_reverse.bedGraph'))
  json_data <- fromJSON(file=paste0(folder_path, 'unique_read_count.json'))
  
  library_size <- json_data$valid_reads_forward + json_data$valid_reads_reverse
  strand_coverages <- list('+' = coverage(bed_graph_forward, weight = bed_graph_forward$score / library_size * 1e6), 
                           '-' = coverage(bed_graph_reverse, weight = bed_graph_reverse$score / library_size * 1e6))
  
  
  slope <- c()
  intercept <- c()
  p_value_slope <- c()
  p_value_intercept <- c()
  avg_coverage <- c()
  
  for (i in 1:nrow(sj_selected)){
    if (i%%2500==0){
      print(i)
    }
    
    row <- sj_selected[i,]
    intron_start <- row[['start']]
    intron_end <- row[['end']]
    intron_strand <- row[['strand']]
    intron_chromosome <- row[['chromosome']]
    
    intron_coverage_rle <- window(strand_coverages[[intron_strand]][[intron_chromosome]], 
                                  start=intron_start+padding,
                                  end=intron_end-padding)
    intron_coverage <- as.numeric(intron_coverage_rle)
    if (mean(intron_coverage)==0){
      slope <- c(slope, NA)
      intercept <- c(intercept, NA)
      p_value_slope <- c(p_value_slope, NA)
      p_value_intercept <- c(p_value_intercept, NA)
      avg_coverage <- c(avg_coverage, mean(intron_coverage))
      next
    }

    #intron_coverage_normalized <- intron_coverage / mean(intron_coverage)
    
    x <-if (intron_strand=='+') 1:length(intron_coverage) else length(intron_coverage):1
    
    model <- lm(intron_coverage~x)
    model_summary <- summary(model)
    
    slope <- c(slope, model_summary$coefficients["x",'Estimate'])
    intercept <- c(intercept, model_summary$coefficients["(Intercept)",'Estimate'])
    
    p_value_slope <- c(p_value_slope, model_summary$coefficients["x",'Pr(>|t|)'])
    p_value_intercept <- c(p_value_intercept, model_summary$coefficients["(Intercept)",'Pr(>|t|)'])
    avg_coverage <- c(avg_coverage, mean(intron_coverage))
  }
  
  result_df <- data.frame(slope = slope,
                          intercept = intercept,
                          p_value_slope = p_value_slope,
                          p_value_intercept = p_value_intercept,
                          avg_coverage = avg_coverage)
  
  concat_df <- cbind(sj_selected, result_df)
  write.table(concat_df, paste0(output_folder, sample_name, '.tsv'), sep = "\t", row.names = F)
}

