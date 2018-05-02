format.fasta <- function(df) {
  sequence_df <- setDT(df, keep.rownames = TRUE)[]
  colnames(sequence_df) <- c("OTU_ID", "Sequence")
  sequence_df <- as.data.frame(sequence_df)
  sequence_df[,2] <- as.character(sequence_df[,2])
  sequence_df
}
normaliz_16s <- function(biom_tax_df){
  input_table <- biom_tax_df
  singleton_test <- apply(input_table[,-(1:2)], 1, sum) > 1
  input_table_S <- input_table[singleton_test,]
  
  #normalizing the table
  normalized_table <- input_table_S[1:2]
  for(i in 3:length(colnames(input_table_S))){
    #print(colnames(input_table_S[i]))
    normalized_table[i] <- as.data.frame(input_table_S[i]/sum(input_table_S[i]))
  }
  return(normalized_table)
}
subtabling <- function(df, cumul_cutoff){
  d_biom_n_sorted <- df
  d_biom_n_sorted$raw_sum <- apply(d_biom_n_sorted[-(1:2)], 1, sum)
  d_biom_n_sorted <- d_biom_n_sorted[order(-d_biom_n_sorted$raw_sum),]
  
  total_freq_sum <- sum(d_biom_n_sorted$raw_sum)
  d_biom_n_sorted$total_freq <- (d_biom_n_sorted$raw_sum/total_freq_sum) *100
  d_biom_n_sorted$cum_freq <- cumsum(d_biom_n_sorted$total_freq)
  
  d_biom_n_cut <- d_biom_n_sorted[c(d_biom_n_sorted$cum_freq <= cumul_cutoff),]
  drops <- c("raw_sum","total_freq","cum_freq")
  d_biom_n_cut_final <- d_biom_n_cut[ , !(names(d_biom_n_cut) %in% drops)]
  return(d_biom_n_cut_final)
}
mason_similarity <- function(df){
  sequence_df <- df
  colnames(sequence_df) <- c("OTU", "Sequence")
  
  ###################################################################################
  # RUNNING MASON SCRIPT ON THE TABLE & EXPORTING THE MARTIX##################################################
  ###################################################################################
  
  # Create a vector of all sequences.  Sequence position in vector corresponds with row number from sequence_df
  sequence <- c()
  for(i in 1:nrow(sequence_df)){  
    sequence <- c(sequence, sequence_df[i,2])       # Filling vector with sequences 
  }
  nrow(sequence_df) == length(sequence)             # Should return TRUE
  sequence <- unlist(sequence)
  sequence_df[]<-lapply(sequence_df, as.character)
  
  # Pre allocating matrix to fill with data.  Columns is -1 to remove the diagonal.
  perc_mat <- matrix(data = NA, nrow = length(sequence), 
                     ncol = (length(sequence)-1))
  
  
  # Applying correct row and column names for comparison matrix based off of the row and column names of sequence_df
  rownames(sequence_df) <- sequence_df$OTU
  sequence_df$OTU <- NULL
  row.names(perc_mat) <- row.names(sequence_df)
  col_names <- as.vector(row.names(sequence_df))
  col_names <- col_names[-1] 
  colnames(perc_mat) <- col_names
  
  
  # Run sequence comparisons.
  x <- 1          # sequence in position 1 of sequences vector
  y <- 2          # sequence in position 2 of sequences vector
  while(x < length(sequence)){
    
    seq1 <- DNAString(sequence[x])              # Convert sequences in vector to DNAString
    seq2 <- DNAString(sequence[y])
    palign1 <- pairwiseAlignment(seq1, seq2)    # Perform pairwise alignment on DNAStrings
    
    perc_mat[x,(y-1)] <- pid(palign1, type = "PID1")   # Convert pairwise alignment score to percentage and store in matrix
    
    # Once sequence in position 1 has been compared to all other sequences, move to sequence in position 2 
    # to be compared to all other sequences.
    
    if( y == length(sequence)){   
      x <- x + 1  
      y <- x + 1 
    }
    
    # If sequence in position 1 hasn't been compared to all other sequences, keep moving through the vector
    # until it's been compared to all others.
    
    else {
      y <- y + 1
      print(y)           # to check to make sure code is still running. 
    }
    
  }
  return(perc_mat)
}
#########################

setwd("H:/Vyshenska/Rwd")
library(seqinr)
library(Biostrings)
library(data.table)

#############################
#loading files
#Q2 files from Dasha
d_fasta.file <- "./fasta_files/coh1_deblur_sequences.fasta"
d_biom.file <- "./biom_files/d_coh1_deblur_biom.txt"

#Q1 files from Khiem
k_fasta.file <- "./fasta_files/97_otus.fasta"
k_biom.file <- "./biom_files/k_coh1_biom.txt"
k_taxonomy.file <-"./taxonomy/k_97_otu_taxonomy.txt"

#output file names
d_full_normalized <- "d_biom_tax_n.csv"
k_full_normalized <- "k_biom_tax_n.csv"
d_biom_cut <- "d_biom_tax_n_cut.csv"
output_matrix_name <- "deblur_matching_martix.csv"
full_matched_pairs_filename <- "matched_pairs_fullTable.csv"
collapsed_matched_pairs <- "collapsed_matched_pairs.csv"

#conditions
cumul_cutoff <- 95   #% of cumulative OTU total frequencies to pass

############################
#formatting files
d_fasta <- as.data.frame(do.call(rbind, (read.fasta(d_fasta.file,
                                                        as.string = TRUE,
                                                        forceDNAtolower = FALSE))))
d_fasta <- format.fasta(d_fasta)


k_fasta <- as.data.frame(do.call(rbind, (read.fasta(k_fasta.file,
                                                    as.string = TRUE,
                                                    forceDNAtolower = FALSE))))
k_fasta <- format.fasta(k_fasta)

#d_biom file already contains taxonomy as the first column. pay attention to the 
#header too! remove any spaces or # as thay will mess up with the import of the table!
d_biom <- read.table (d_biom.file , header = T, check.names = F, sep = "\t")
colnames(d_biom)[1] <- c("Taxonomy_d")
#k_biom does not contain taxonomy, so I add it here after importing biom export file
k_biom <- read.table (k_biom.file , header = T, check.names = F)
k_taxonomy <- read.table (k_taxonomy.file, header = T, sep = "\t")
colnames(k_taxonomy) <- c("OTU_ID", "Taxonomy_k")
k_biom_tax <- merge(k_taxonomy,k_biom , by ="OTU_ID")


##########################

#shortening biom_tax tables to contain 99% accumulative by frequency OTUs only

# D table:removing singletons, normalize table; output full normalized table; selecting OTUs based on total (across samples) cumulative frequency cutoff, adding sequence and outputing this short table
d_biom_n <- normaliz_16s(d_biom)
write.csv(d_biom_n, d_full_normalized, row.names = F) 
d_biom_n_cut_f <- subtabling(d_biom_n, cumul_cutoff)
  #adding sequence to selected otus into the table, exporting it
d_biom_n_cut_f_seq <- merge(d_fasta, d_biom_n_cut_f)
write.csv(d_biom_n_cut_f_seq, "d_biom_n_cut_f_seq.csv", row.names = F)

# K table: removing singletons, normalize table; output full normalized table; selecting OTUs based on total (across samples) cumulative frequency cutoff, adding sequence and outputing this short table
  #adding taxonomy to K table biom###########################################

k_biom_n <- normaliz_16s(k_biom_tax)
write.csv(k_biom_n, k_full_normalized, row.names = F) 
k_biom_n_cut_f <- subtabling(k_biom_n, cumul_cutoff)
  #adding sequence to selected otus into the table, exporting it
k_biom_n_cut_f_seq <- merge (k_fasta, k_biom_n_cut_f)
write.csv(k_biom_n_cut_f_seq, "k_biom_n_cut_f_seq.csv", row.names = F)


#preparing files for merging to feed into sequence comparison function
d_otu_toMerge <- d_biom_n_cut_f_seq[,1:2]
k_otu_toMerge <- k_biom_n_cut_f_seq[,1:2]
bind_df <- rbind(d_otu_toMerge,k_otu_toMerge)

#accomodate Mason function and run it on subsampled rbind table; export matrix as csv
perc_mat_output <- mason_similarity(bind_df)

#perc_mat_output <- read.csv("dada2_matching_martix.csv", check.names = F, row.names = 1)

# Insert file path and name for output
write.csv(perc_mat_output, file= output_matrix_name)

####################################################################################################################
#perc_mat_output <- perc_mat ##remove this line for the next analysis as perc_mat now is local for mason function! and your real output went to perc_mat_output!
####################################################################################################################


#cutting the score matrix to contain only meaningful scores
d_merged_length <- length(d_otu_toMerge$OTU_ID)
score_matrix <- perc_mat_output[1:d_merged_length,(d_merged_length):length(colnames(perc_mat_output))]
max_score_vec <- apply(score_matrix, 1, max) #maximal score vector where I look into which Khiem OTUs match to any given D OTU

bul_matrix <- matrix(nrow = d_merged_length, ncol = length(k_otu_toMerge$OTU_ID))
for (i in 1:length(max_score_vec)){
  bul_matrix[i,] <- score_matrix[i,] == max_score_vec[i]
}
k_OTU_counts <- apply(bul_matrix, 1, sum)
length(k_OTU_counts)


#building a table of all matched pairs based on search for hits for each of D OTU in K OTUs (only one d_OTU can be matched with several d_OTUs but not vice versa)

cnames_sm <- colnames(score_matrix)
rnames_sm <- rownames(score_matrix)

#score_matrix <- as.matrix(score_matrix)
matched_pairs_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("dOTU_ID", "kOTU_ID", "match_score", "dOTU_match_#"))

#

for(i in 1:d_merged_length){
  
  #print(rep(rnames_sm[i], length(cnames_sm[which(bul_matrix[i,])])))
  #print(cnames_sm[which(bul_matrix[i,])])
  #print(score_matrix[i,bul_matrix[i,]])
  d_temp <- c(rep(rnames_sm[i], length(cnames_sm[which(bul_matrix[i,])])))
  #print(d_temp)
  df <- as.data.frame(d_temp)
  #print(df)
  k_temp <- c(cnames_sm[which(bul_matrix[i,])])
  #print(k_temp)
  df$k_temp <- cbind(k_temp)
  #print(df)
  df$score <- score_matrix[i,bul_matrix[i,]]
  #print(df)
  df$match_num <- k_OTU_counts[i]
  #print(df)
  matched_pairs_df <- rbind(matched_pairs_df, df)
  #print(matched_pairs_df)
  rm(df)

}
#adding taxonomy, sequence and abandance to the table

d_m_temp <- d_biom_n_cut_f_seq
d_m_temp$aban_sum <- apply(d_m_temp[-(1:3)], 1, sum)
d_m_temp$aban_allSamp <- d_m_temp$aban_sum/sum(d_m_temp$aban_sum)
d_m_temp2 <- d_m_temp[c("OTU_ID", "aban_allSamp", "Taxonomy_d", "Sequence")]
colnames(d_m_temp2) <- c("OTU_ID", "aban_sum_d", "Taxonomy_d", "Sequence_d")

k_m_temp <- k_biom_n_cut_f_seq
k_m_temp$aban_sum <- apply(k_m_temp[-(1:3)], 1, sum)
k_m_temp$aban_allSamp <- k_m_temp$aban_sum/sum(k_m_temp$aban_sum)
k_m_temp2 <- k_m_temp[c("OTU_ID", "aban_allSamp", "Taxonomy_k", "Sequence")]
colnames(k_m_temp2) <- c("OTU_ID", "aban_sum_k", "Taxonomy_k", "Sequence_k")

colnames(matched_pairs_df) <- c("dOTU_ID", "kOTU_ID", "match_score", "dOTU_match_#")

full_matched_pairs_df <- merge(matched_pairs_df, d_m_temp2, by.x = "dOTU_ID", by.y = "OTU_ID")
full_matched_pairs_df <- merge(full_matched_pairs_df, k_m_temp2, by.x = "kOTU_ID", by.y = "OTU_ID")
full_matched_pairs_df <- full_matched_pairs_df[, c(2,1,3:10)]
write.csv(full_matched_pairs_df, file= full_matched_pairs_filename, row.names = F)

#sum(full_matched_pairs_df$aban_sum_k[full_matched_pairs_df$`dOTU_match_#`>1])


#collapsing match table and k_OTU table based on unique d_OTU IDs.
#in the final table for the d_OTU IDs that matched with several k_OTU IDs:
# for k_abandance - the sum of abandances fro all matched k_OTUs for a given d_OTU
# for k_OTU ID and Taxonomy_k - the OTU data that has the max abandance among these that
# matched the given d_OTU

sorted_matched_pairs <- full_matched_pairs_df[order(full_matched_pairs_df[,4], full_matched_pairs_df[,1]),]
coll_mp <- setNames(data.frame(matrix(ncol = length(sorted_matched_pairs), nrow = 0)), colnames(sorted_matched_pairs))
temp_coll_mp <- setNames(data.frame(matrix(ncol = length(sorted_matched_pairs), nrow = 0)), colnames(sorted_matched_pairs))

for (i in 1:length(sorted_matched_pairs$`dOTU_match_#`)){
  if(sorted_matched_pairs$`dOTU_match_#`[i] == 1){
    coll_mp <- rbind(coll_mp, sorted_matched_pairs[i,])
} else {
  temp_coll_mp <- rbind(temp_coll_mp, sorted_matched_pairs[i,])
}}

temp <- setNames(data.frame(matrix(ncol = length(sorted_matched_pairs), nrow = 0)), colnames(sorted_matched_pairs))
u_dOTU_vec <- unique(temp_coll_mp$dOTU_ID)
k_aband1 <- temp_coll_mp$aban_sum_k
k_tax_vec <- temp_coll_mp$Taxonomy_k
k_otu_vec <- temp_coll_mp$kOTU_ID
count1 <- 1
for(i in 1:length(u_dOTU_vec)){
  #print(temp_coll_mp$`dOTU_match_#`[count1])
  temp <- rbind(temp, temp_coll_mp[count1,])
  temp$aban_sum_k[i] <-  sum(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)])
  #print(max(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)]))
  #print(which(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)] == max(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)])))
  temp$Taxonomy_k[i] <-  k_tax_vec[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)][which(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)] == max(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)]))]
  #print(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)])
  temp$kOTU_ID[i] <- k_otu_vec[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)][which(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)] == max(k_aband1[count1:(count1+temp_coll_mp$`dOTU_match_#`[count1] - 1)]))]
  count1 <- count1+temp_coll_mp$`dOTU_match_#`[count1]
}
coll_mp <- rbind(coll_mp, temp)
#sum(length(unique(coll_mp$dOTU_ID)))
#sum(coll_mp$`dOTU_match_#`)
write.csv(coll_mp, file= collapsed_matched_pairs, row.names = F)

#now collapsing k_OTU count table based on matched pairs table

c_smp <- sorted_matched_pairs
c_smp$cOTU_ID = paste(c_smp$dOTU_ID, c_smp$kOTU_ID, sep="_")

csmp1 <- setNames(data.frame(matrix(ncol = length(c_smp), nrow = 0)), colnames(c_smp))
csmp2 <- setNames(data.frame(matrix(ncol = length(c_smp), nrow = 0)), colnames(c_smp))

for (i in 1:length(c_smp$`dOTU_match_#`)){
  if(c_smp$`dOTU_match_#`[i] == 1){
    csmp1 <- rbind(csmp1, c_smp[i,])
  } else {
    csmp2 <- rbind(csmp2, c_smp[i,])
  }}

col_kOTU_biom <- setNames(data.frame(matrix(ncol = length(k_biom_n_cut_f_seq), nrow = 0)), colnames(k_biom_n_cut_f_seq))
for (i in 1:length(csmp1$dOTU_ID)){
  col_kOTU_biom <- rbind(col_kOTU_biom, k_biom_n_cut_f_seq[which(k_biom_n_cut_f_seq$OTU_ID == csmp1$kOTU_ID[i]),])
}
col_kOTU_biom$cOTU_ID <- csmp1$cOTU_ID
#write.csv(col_kOTU_biom, file= "col_kOTU_biom.csv", row.names = F)

#now we work with csmp2 part of the table
cu_vec <- unique(csmp2$dOTU_ID)
c_temp <- setNames(data.frame(matrix(ncol = length(k_biom_n_cut_f_seq), nrow = length(cu_vec))), colnames(k_biom_n_cut_f_seq))

for (i in 1:length(cu_vec)){
  working_vec <- which(csmp2$dOTU_ID == cu_vec[i])
  working_df <- csmp2[working_vec,]
  id_populate <- working_df$kOTU_ID[which(working_df$aban_sum_k == max(working_df$aban_sum_k))]

  working_df2 <- k_biom_n_cut_f_seq[(k_biom_n_cut_f_seq$OTU_ID %in% working_df$kOTU_ID),]
  c_temp[i, 1:2] <- working_df2[which(working_df2$OTU_ID == id_populate),1:2]
  c_temp[i,3] <- as.character(working_df2[which(working_df2$OTU_ID == id_populate),3])
  c_temp[i, 4:length(colnames(c_temp))] <- apply(working_df2[-(1:3)], 2, sum)
  rm (working_vec, working_df, id_populate, working_df2)}
c_temp$cOTU_ID <- paste(cu_vec, c_temp$OTU_ID, sep = "_") 
coll_kOTU_table <- rbind(col_kOTU_biom,c_temp)
write.csv(coll_kOTU_table, file= "coll_kOTU_table.csv", row.names = F)

###################################################
#next, I work with two tables with counts (k and d OTU tables) to calculate
# spearman correlations and p-values
#FIRST - MAKING THE TABLES TO LOOK SAME

# my computer had to be rebooted so I lost environment, that is why right now
# I am importing three files that I need for correlation calculations.
#however, if you run the code - you will have these files in your environment
# for khiem's tables - use coll_kOTU_table
# for dariia's tables - use d_biom_n_cut_f_seq 

#d_biom_n_cut_f_seq <- read.csv("d_biom_n_cut_f_seq.csv", header=T, check.names = F)
#collapsed_matched_pairs <- read.csv("collapsed_matched_pairs.csv", header=T, check.names = F)
#coll_kOTU_table <- read.csv("coll_kOTU_table.csv", header = T, check.names = F)
coll_mp$cOTU_ID <- paste(coll_mp$dOTU_ID, coll_mp$kOTU_ID, sep = "_")
temp_IDs <- coll_mp[c("dOTU_ID", "cOTU_ID")]
coll_dOTU_table <- merge(d_biom_n_cut_f_seq, temp_IDs, by.x = "OTU_ID", by.y = "dOTU_ID")

d_table <- coll_dOTU_table[c(length(coll_dOTU_table), 4:(length(coll_dOTU_table) - 1))]
a <- colnames(d_table[2:(length(d_table))])
a <- as.numeric(gsub('d', '', a))
colnames(d_table) <- as.character(c(colnames(d_table[1]), a))
d_table <- d_table[order(d_table$cOTU_ID),]

k_table <- coll_kOTU_table[c(length(coll_kOTU_table), 4:(length(coll_kOTU_table) - 1))]
b <- colnames(k_table[2:(length(k_table))])
b <- as.numeric(gsub('k', '', b))
colnames(k_table) <- as.character(c(colnames(k_table[1]), b))
k_table <- k_table[order(k_table$cOTU_ID),]
similarity_df <- coll_mp[,c("cOTU_ID", "match_score")]
similarity_df <- similarity_df[order(similarity_df$cOTU_ID),]
similarity_df$above_95 <- similarity_df$match_score >= 95
similarity_df$below_95 <- similarity_df$match_score < 95
#sum(k_table$cOTU_ID == d_table$cOTU_ID && d_table$cOTU_ID == similarity_df$cOTU_ID)
#sum(d_table$cOTU_ID == similarity_df$cOTU_ID)
#write.csv(d_table, file= "d_table.csv", row.names = F)
#write.csv(k_table, file= "k_table.csv", row.names = F)
#write.csv(similarity_df, file= "similarity_df.csv", row.names = F)

#SOCOND - creating two output tables and calculating correlations

samples <- colnames(d_table[-1])
rho_table <- setNames(data.frame(c("full","above95","below95")), c("gene_set"))
pval_table<- setNames(data.frame(c("full","above95","below95")), c("gene_set"))

for(i in 1:length(samples)){
print(samples[i])

cor_result_full <- cor.test(d_table[,c(samples[i])], k_table[,c(samples[i])], method = c("spearman"), exact = F)
rho_full <- cor_result_full$estimate
pval_full <- cor_result_full$p.value

cor_result_above95 <- cor.test(d_table[similarity_df$above_95,c(samples[i])], k_table[similarity_df$above_95,c(samples[i])], method = c("spearman"), exact = F)
rho_above95 <- cor_result_above95$estimate
pval_above95 <- cor_result_above95$p.value

cor_result_below95 <- cor.test(d_table[similarity_df$below_95,c(samples[i])], k_table[similarity_df$below_95,c(samples[i])], method = c("spearman"), exact = F)
rho_below95 <- cor_result_below95$estimate
pval_below95 <- cor_result_below95$p.value


rho_table[,c(samples[i])] <- c(rho_full, rho_above95, rho_below95)
pval_table[,c(samples[i])] <- c(pval_full, pval_above95, pval_below95)
rm(rho_full, rho_above95, rho_below95, pval_full, pval_above95, pval_below95, cor_result_full, cor_result_above95, cor_result_below95)

##

pdf_name <- paste(samples[i],"pdf", sep = ".")
pdf(pdf_name)
plot(d_table[,c(samples[i])], k_table[,c(samples[i])], xlab = "QIIME2", ylab = "QIIME1", main = "All samples")
abline(lm(d_table[,c(samples[i])]~ k_table[,c(samples[i])]))

plot(d_table[similarity_df$above_95,c(samples[i])], k_table[similarity_df$above_95,c(samples[i])], xlab = "QIIME2", ylab = "QIIME1", main = "Samples: similarity >= 95%")
abline(lm(d_table[similarity_df$above_95,c(samples[i])] ~ k_table[similarity_df$above_95,c(samples[i])]))

plot(d_table[similarity_df$below_95,c(samples[i])], k_table[similarity_df$below_95,c(samples[i])], xlab = "QIIME2", ylab = "QIIME1", main = "Samples: similarity < 95%")
abline(lm(d_table[similarity_df$below_95,c(samples[i])] ~ k_table[similarity_df$below_95,c(samples[i])]))
dev.off()
rm(pdf_name)

}
write.csv(rho_table, file= "rho_table.csv", row.names = F)
write.csv(pval_table, file= "pval_table.csv", row.names = F)

