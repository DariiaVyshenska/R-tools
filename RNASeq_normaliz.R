#this code takes count table in .csv format (first column - gene IDs, all others - counts for samples)
#it generates another table with normalized counts (log2(xi/sum(x))*10^6+1)

getwd()
setwd("D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd")


table_full <- read.csv("./CC_repl_sum_count.csv")
#colnames(input_table[1])
input_table <- table_full [-1:-5,]


gene_hits <- numeric(length(colnames(input_table))-1)
gene_hits10 <- numeric(length(colnames(input_table))-1)


normalized_table <- input_table[1]
for(i in 2:length(colnames(input_table))){
  print(colnames(input_table[i]))
  normalized_table[i] <- log2((input_table[i]/sum(input_table[i])*10^6+1))
  gene_hits[i-1] <- sum(ifelse(input_table[i] > 0, 1, 0))
  gene_hits10[i-1] <- sum(ifelse(input_table[i] > 9, 1, 0))
  
}

#here I write to the dataframe and then into the file how many genes were detected - 
#per sample - min, max, mean and median

stat_df <- data.frame(c("gene_hits","gene_hits_min10"))
stat_df[2] <-c(min(gene_hits), min(gene_hits10)) 
stat_df[3] <-c(max(gene_hits), max(gene_hits10)) 
stat_df[4] <-c(mean(gene_hits), mean(gene_hits10)) 
stat_df[5] <-c(median(gene_hits), median(gene_hits10)) 
names(stat_df) <- c("X", "min", "max", "mean", "median")

write.csv(normalized_table, "RNASeq_counts_normalized.csv", row.names = F) 
write.csv(stat_df, "stat_gene_counts.csv", row.names = F) 

