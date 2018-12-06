library(metap)
# Function: extracts rows based on key word (value variable) - use "Coefficient"
#or "pvalue"
extract_table <- function(input_table, value){
  new_table <- input_table[, (grepl(value, colnames(input_table)))]
  return(new_table)
}

# Function: filters table either based on p-value ("p_val") threshold or on 
# corr coef directionality ("coef_dir"; adds coefition direction column where 1
# is positive, -1 is negative correlation). 
# If any NA is detected - row is excluded

subtable_table <- function(test_table, test, p_val_threshold = NULL){  
  if(test == "coef_dir"){
    test_value <- "Coefficient"
  } else if (test == "p_val"){
    test_value <- "pvalue"
  }
  subtable <- extract_table(test_table, test_value)
  # checking for NA's among the data to exclude these rows in the future
  na_vec <- apply(subtable, 1, function(x) any(is.na(x)))

  if(test == "coef_dir"){
    subtable_vec <- apply(subtable, 1, function(x) length(unique(sign(x))) == 1)
    if (sum(na_vec) == 0){
      cat("No NA in correlations!\n")
      return_table <- test_table[subtable_vec, ]
      return_table$cor_direction <- sign(subtable[subtable_vec,1])
    } else if (sum(na_vec) > 0){
      cat("NAs were found in correlations!\n")
      return_table <- test_table[subtable_vec & !(na_vec), ]
    }
  } else if(test == "p_val"){
    cat("This part of Script needs to be written!\n")
    # take min of all pvals and then apply threshold to the minimal
  }
  return(return_table)
}  


# Function: calculates Fisher's combined probability test p-value and
# attaches this p-value as extra column "fish_pval". If parameter "fdr" = T
# than fdr for fisher will be calcualted and added as "fisher_fdr" column.
# Returns table
fisher_calc <- function(test_table, fdr = F){
  # extracting rows with either pvalues or corr coefficients
  subtable <- extract_table(test_table, "pvalue")
  # replace any zeros
  subtable[subtable==0.0000000]<-0.00000001
  # calculate fisher
  test_table$fish_pval <- apply(subtable, 1, function(x) if(any(is.na(x))){NA}else(sumlog(x)$p))
  if (fdr == T){test_table$fish_fdr <- p.adjust(test_table$fish_pval, 
                                                method = "fdr")}
  return(test_table)
}
