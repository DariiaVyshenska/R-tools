getwd()
setwd("D:/GitHub/cc_16s/ANALYSIS-BY-STAGE")
#install.packages("metap")
library(metap)
# function that correlates frequency to hypoxia, STAGE as parameter
# and creates a table per analysis

analyze_by_cohort <- function(stage){
  # filtering pre-created tables by stage: metadata (md) and frequency (freq)
  # tables
  if(length(stage) == 1){
    md_coh12 <- hyp_sam_md_sorted[hyp_sam_md_sorted$FIGO_stage == stage , 
                                  c(1, 5, 9)]
  } else{
    md_coh12 <- hyp_sam_md_sorted[hyp_sam_md_sorted$FIGO_stage == stage[1] | 
                                    hyp_sam_md_sorted$FIGO_stage == stage[2], 
                                  c(1, 5, 9)]
  }
  md_coh1 <- md_coh12[md_coh12$Cohort == 1, ]
  md_coh2 <- md_coh12[md_coh12$Cohort == 2, ]
  
  freq_coh1 <- freq_table_hyp_sorted[ ,(colnames(freq_table_hyp_sorted) 
                                        %in% md_coh1$PatientID)]
  freq_coh2 <- freq_table_hyp_sorted[ ,(colnames(freq_table_hyp_sorted) 
                                        %in% md_coh2$PatientID)]
  
  # build empty output table:
  
  corr_val_table <- setNames(data.frame(matrix(ncol = 8, 
                            nrow = length(row.names(freq_table_hyp_sorted)))), 
                            c("feature", "median-coh1", "rho-coh1", 
                              "p-val-coh1", "median-coh2","rho-coh2", 
                              "p-val-coh2", "fish_pval"))
  
  # calculate correlations, medians, and fdr
  for(i in 1:length(row.names(freq_coh1))){
    print(i)
    median_coh1 <- median(as.numeric(freq_coh1[i,]))
    median_coh2 <- median(as.numeric(freq_coh2[i,]))
    
    cor_result1 <- cor.test(as.numeric(freq_coh1[i,]), md_coh1[,3], 
                            method = c("spearman"), exact = F)
    cor_result2 <- cor.test(as.numeric(freq_coh2[i,]), md_coh2[,3], 
                            method = c("spearman"), exact = F)
    feature <- row.names(freq_coh1[i,])
    rho1 <- cor_result1$estimate
    pval1 <- cor_result1$p.value
    
    rho2 <- cor_result2$estimate
    pval2 <- cor_result2$p.value
    
    if (is.na(pval1) | is.na(pval2)){
      fish_pval <- "NA"
    } else{
      fish_results <- sumlog(c(pval1, pval2))
      fish_pval <- fish_results$p 
    }
    
    corr_val_table[i,] <- c(feature, median_coh1, 
                            rho1, pval1, median_coh2, rho2, pval2, fish_pval)
  }
  
  # metaanalysis: if correlation has the same direction across cohorts - 1, 
  # if not - 0
  corr_val_table$corr_dir <- as.numeric(corr_val_table$`rho-coh1` > 0 & 
                                          corr_val_table$`rho-coh2`> 0 |
                                          corr_val_table$`rho-coh1`< 0 &
                                          corr_val_table$`rho-coh2`< 0)
  corr_val_table$`fish<0.3` <- as.numeric(corr_val_table$fish_pval < 0.3)
  
  # updating column names for the output table
  if(length(stage) == 1){
    col_temp <- paste("stage", stage, sep = "")
  } else{
    stage_paste <- paste(stage[1], stage[2], sep = "")
    col_temp <- paste("stage", stage_paste, sep = "")
  }
  
  new_col_name <- paste (colnames(corr_val_table[-1]), col_temp, sep = "_")
  colnames(corr_val_table) <- c ("feature", new_col_name)
  
  if(length(stage) == 1){
    file_path <- paste("corr_tax_table", "_stage", stage,".csv", sep = "")
  }else{
    file_path <- paste("corr_tax_table", "_stage", stage_paste,".csv", sep = "")
  }
  write.csv(corr_val_table, file= file_path, row.names = F)
  return(corr_val_table)
}

#############
# importing data
count_table234 <- read.table(file = './normalized_table234.csv', 
                             sep = ',', header = TRUE, check.names = F)

taxonomy234 <- read.csv(file = './taxonomy234.csv', sep = ",", 
                        header = T, check.names = F)
taxonomy234 <- taxonomy234[,-3]
colnames(taxonomy234) <- c('feature', 'taxonomy')

metadata <- read.csv(file = './metadata.csv', sep = ",", 
                        header = T, check.names = F)

# create two tables - metadata and frequency table 
# for hypoxia available samples
hyp_sam_metadata <- metadata[which(metadata$HypFraction != "NA"),]
hyp_sam_md_sorted <- hyp_sam_metadata[ order(hyp_sam_metadata$PatientID), , drop = F]

freq_table_hyp <- count_table234[, c (1, which(colnames(count_table234) %in% hyp_sam_md_sorted$PatientID))]

freq_table_hyp2 <- as.data.frame(freq_table_hyp[,-1])
rownames(freq_table_hyp2) <- freq_table_hyp[,1]
freq_table_hyp <- freq_table_hyp2
freq_table_hyp_sorted <- freq_table_hyp[ , order(names(freq_table_hyp))]

# two tables to work with: hyp_sam_md_sorted and freq_table_hyp_sorted
# checking tables:
#hyp_sam_md_sorted$PatientID == colnames(freq_table_hyp_sorted)

# analysis of frequency tables: correlation of each bacteria_feature 
# with hypoxia
stage2_resutls <- analyze_by_cohort(2)
stage3_results <- analyze_by_cohort(3)
stage23_results <- analyze_by_cohort(c(2,3))

# merge all analysis tables into one big table 
all_results0 <- merge(taxonomy234, stage2_resutls, by.x = "feature", 
                      by.y = "feature")
all_results1 <- merge(all_results0, stage3_results, by.x = "feature", 
                      by.y = "feature")
all_results <- merge(all_results1, stage23_results, by.x = "feature", 
                     by.y = "feature")

# exporting merged final table
write.csv(all_results, file= "./merged_all_prelim_analysis.csv", row.names = F)
