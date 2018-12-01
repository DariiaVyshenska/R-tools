# Function: extracts rows based on key word (value variable) - use "Coefficient"
#or "pvalue"
extract_table <- function(input_table, value){
  new_table <- input_table[, (grepl(value, colnames(input_table)))]
  return(new_table)
}

# Function: filters table either based on p-value ("p_val") threshold or on 
# corr coef directionality ("coef_dir"). If any NA is detected - row is excluded

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
    } else if (sum(na_vec) > 0){
      cat("NAs were found in correlations!\n")
      return_table <- test_table[subtable_vec & !(na_vec), ]
    }
  } else if(test == "p_val"){
    print(" This part of Script needs to be written!")
    # take min of all pvals and then apply threshold to the minimal
  }
}  
  
