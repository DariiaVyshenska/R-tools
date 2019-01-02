library(metap)
# Function: extracts rows based on key word (value variable) - use "Coefficient"
#or "pvalue"
extract_table <- function(input_table, value){
  new_table <- input_table[, (grepl(value, colnames(input_table)))]
  return(new_table)
}

# Function: filters table either based on p-value ("p_val") threshold or on 
# corr coef directionality; adds coefition direction column where 1
# is positive, -1 is negative correlation ("coef_dir"); adds median of correlation
# coefitient("corcoef_median").

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
      return_table$corcoef_median <- apply(subtable[subtable_vec,],1,function(x) median(x))
    } else if (sum(na_vec) > 0){
      cat("NAs were found in correlations!\n")
      return_table <- test_table[subtable_vec & !(na_vec), ]
      return_table$cor_direction <- sign(subtable[subtable_vec & !(na_vec),1])
      return_table$corcoef_median <- apply(subtable[subtable_vec & !(na_vec),],1,function(x) median(x))
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


# Function: imports and merges separate correlation analysis results files 
# into one big table (plus exporting it as .csv file as a backup; a must for
# large analysis files)
# Input: a path to the folder with ONLY correlation analysis files (Richard's
# script output; one file = one analysis) - example 
# input_path = "D:/GitHub/R_wd/corr_raw_input/"
# Output: single merged table with merged analysis, analysis headers renamed
# according to the analysis # (from analysis 1 for all to analysis 1, 2, 3 ect.)
# Function also exports the merged table as .csv (mandatory)
merge_raw_corr <- function(input_path){
  #getting the list of files in the specified folder
  files <- list.files(input_path)
  file_path <- paste(input_path, files[1], sep = "")
  #creating the main table, importing the first file with analysis
  merged <- read.table(file_path, header = T, check.names = F, sep = ",")
  # looping through the rest of the files, importing them and merging with
  # our main merging table
  for(i in 2:length(files)){
    file_path <- NA
    file_path <- paste(input_path, files[i], sep = "")
    tableToAdd <- read.table(file_path, header = T, check.names = F, sep = ",")
    
    analysis <- paste("Analys", i, sep = " ")
    header <- names(tableToAdd)
    newheader <- gsub("Analys 1", analysis, header)
    names(tableToAdd) <- newheader
    merged <- merge(merged, tableToAdd, by.x = "pairName", 
                    by.y = "pairName")
  }
  # exporting intact but already merged analysis files as a single table
  write.csv(merged, "merged_analysisTable.csv", row.names = F)
  
  return(merged)
}