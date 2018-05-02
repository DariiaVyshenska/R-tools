setwd("D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd/Q1vsQ2")
input_file_name <- "d_coh1_dada2_biom.txt"


input_table <- read.table(input_file_name, header=T, check.names = F)

#removing singletons
singleton_test <- apply(input_table[,-1], 1, sum) > 1
input_table_S <- input_table[singleton_test,]

#normalizing the table
normalized_table <- input_table_S[1]
for(i in 2:length(colnames(input_table_S))){
  #print(colnames(input_table_S[i]))
  normalized_table[i] <- (input_table_S[i]/sum(input_table_S[i]))
}

#writing table into the file
write.csv(normalized_table, "normalized_table.csv", row.names = F) 


