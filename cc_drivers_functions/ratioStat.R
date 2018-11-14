# this funciton takes matrix of ratios and creates two tables:
# table with all combination for each KD target and ratios for these combinations
# table with statistics - min, max, and max/min ratios for each KD-target
# function returns only statistics, but exports both tables as .csv files
# you also need to provid pval under which the ratios matrix was created
# this pval will be used in file names only.

ratioStat <- function(matrixA, pval){
  kd_tar_sum <- ncol(matrixA)
  min_vec <- vector(length = kd_tar_sum)
  max_vec <- vector(length = kd_tar_sum)
  maxBymin_vec <- vector(length = kd_tar_sum)
  names_vec <- vector(length = kd_tar_sum)
  
  kd_tar_sum <- ncol(matrixA)
  gene_name <- colnames(matrixA)[1]
  print(gene_name)
  all_comps <- paste(gene_name, colnames(matrixA)[1:kd_tar_sum], sep = "vs")
  names_apply <- all_comps[!all_comps == paste(gene_name, gene_name, sep = "vs")]
  ratios <- matrixA[1,c(2:kd_tar_sum)]
  
  min_vec[1] <- min(ratios)
  max_vec[1] <- max(ratios)
  maxBymin_vec[1] <- max(ratios)/min(ratios)
  names_vec[1] <- gene_name
  
  df <- data.frame(names_apply, ratios)
  names(df) <- c(gene_name, paste(gene_name, "ratio", sep = "_"))
  
  rm(all_comps, names_apply, gene_name, ratios)
  a <- 1
  b <- 3
  for(i in 2:(kd_tar_sum - 1)){
    gene_name <- colnames(matrixA)[i]
    cat("Processing ", gene_name, "\n")
    all_comps <- paste(gene_name, colnames(matrixA)[1:kd_tar_sum], sep = "vs")
    names_apply <- all_comps[!all_comps == paste(gene_name, gene_name, sep = "vs")]
    
    ratios <- c(matrixA[c(1:a),i],matrixA[i, c(b:kd_tar_sum)])
    df_temp <- data.frame(names_apply, ratios)
    names(df_temp) <- c(gene_name, paste(gene_name, "ratio", sep = "_"))
    
    min_vec[i] <- min(ratios)
    max_vec[i] <- max(ratios)
    maxBymin_vec[i] <- max(ratios)/min(ratios)
    names_vec[i] <- gene_name
    
    df <- cbind(df, df_temp)
    
    a <- a+1
    b <- b+1
    rm(all_comps, names_apply, ratios, df_temp, gene_name)
  }
  
  gene_name <- colnames(matrixA)[kd_tar_sum]
  print(gene_name)
  all_comps <- paste(gene_name, colnames(matrixA)[1:kd_tar_sum], sep = "vs")
  names_apply <- all_comps[!all_comps == paste(gene_name, gene_name, sep = "vs")]
  ratios <- matrixA[1:(kd_tar_sum - 1), kd_tar_sum]
  df_temp <- data.frame(names_apply, ratios)
  names(df_temp) <- c(gene_name, paste(gene_name, "ratio", sep = "_"))
  df <- cbind(df, df_temp)
  
  min_vec[kd_tar_sum] <- min(ratios)
  max_vec[kd_tar_sum] <- max(ratios)
  maxBymin_vec[kd_tar_sum] <- max(ratios)/min(ratios)
  names_vec[kd_tar_sum] <- gene_name
  stat_df <- data.frame(KDtarget = names_vec, min_ratio = min_vec, 
                        max_ratio = max_vec, "max_ratio/min_ratio" = maxBymin_vec,
                        check.names = F)
  
  
  df_fname <- paste("ratio_tablesAllKDpval", pval, ".csv", sep = "")
  stat_fname <- paste("statOfratiospval", pval, ".csv", sep = "")
  write.csv(df, df_fname, row.names = F)
  write.csv(stat_df, stat_fname, row.names = F)
  return(stat_df)
}