# This function takes list of tables with info about DEGs and pvalue
# threshold, it spits out dataframe with the name of the KD-target in first 
# column and number of genes that pass this pvalue threshold for that gene
# as a second column

degs_df <- function(tables, pval){
  deg_vec <- vector(length =  length(tables))
  for(i in 1:length(tables)){
    print(names(tables[i]))
    deg_vec[i] <- sum(tables[[i]][["Parametric p-value"]] <= pval)
  }
  df <- data.frame("KDtarget" = names(tables), "DEGs" = deg_vec)
  return(df)
}