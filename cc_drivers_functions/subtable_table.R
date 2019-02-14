library(metap)
# Function: extracts rows based on key word (value variable) - use "Coefficient"
#or "pvalue"
extract_table <- function(input_table, value){
  new_table <- input_table[, (grepl(value, colnames(input_table)))]
  return(new_table)
}

# Function: filters table either based on p-value ("p") threshold or on 
# corr coef directionality; adds coefition direction column where 1
# is positive, -1 is negative correlation ("c"); adds median of correlation
# coefitient("corcoef_median") - in the mode of "coef_dir" (encountering any NA leads
# to row exclusion).
# in the mode of "p_val" it only filters the table by individual pvalues

subtable_table <- function(test_table, test, p_val_threshold = NULL){  
  if(test == "c"){
    test_value <- "Coefficient"
  } else if (test == "p"){
    test_value <- "pvalue"
  }
  subtable <- extract_table(test_table, test_value)
  # checking for NA's among the data to exclude these rows in the future
  na_vec <- apply(subtable, 1, function(x) any(is.na(x)))

  if(test_value == "Coefficient"){
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
  } else if(test_value == "pvalue"){
    # selects only rows where maxium pvalue less or equal to pvalue threshold
    min_pval_v <- apply(subtable, 1, max)
    subtable_vec <- as.numeric(min_pval_v) <= p_val_threshold
    return_table <- test_table[subtable_vec, ]
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
  test_table$fish_pval <- apply(subtable, 1, function(x) if(any(is.na(x))){NA}else(sumlog(as.numeric(x))$p))
  if (fdr == T){test_table$fish_fdr <- p.adjust(as.numeric(test_table$fish_pval), 
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
  merged <- read.table(file_path, header = T, check.names = F, sep = ",", 
                       stringsAsFactors = F)
  # looping through the rest of the files, importing them and merging with
  # our main merging table
  for(i in 2:length(files)){
    file_path <- NA
    file_path <- paste(input_path, files[i], sep = "")
    tableToAdd <- read.table(file_path, header = T, check.names = F, sep = ",",
                             stringsAsFactors = F)
    
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

# FUNCTION: calculates PUC
# input: two tables 
# a) fc_pucInput is a data frame, first column is a unique id,
# second column is fold change sign (1 if gene is up, -1 if gene is down). this
# table should not contain any other sumbols, otherwise you will need to change
# the funciton; 
# b) edge_pucInput is a data frame, first column must contain the pairName, 
# second - correlation direction (1 or -1 only!)
# output: vector of puc (1 if pass, 0 if no pass, NA if one or both of nodes in 
# the edge did not have corresponding fold change value in table #1), the order
# is respectful to edges in initial table. can be used dicectly to append as
# a column to initial data frame with edges
puc <- function(fc_pucInput, edge_pucInput){
  #preparing edges
  genes_L <- strsplit(edge_pucInput[,1], split = "<==>", fixed = T)
  edge_pucInput$n1 <- unlist(lapply(genes_L, `[[`, 1))
  edge_pucInput$n2 <- unlist(lapply(genes_L, `[[`, 2))
  edges <- edge_pucInput[,-1]
  
  #prep hash table for fold change references
  fcenv<-new.env()
  for(i in seq(nrow(fc_pucInput))){
    fcenv[[ fc_pucInput[i,1] ]]<- fc_pucInput[i,2]
  }
  
  # getting vector with PUC: 1 if pass, 0 if no pass, NA if one of the nodes
  # is not represented by fold change
  puc_vector <- vector()
  for(i in 1:nrow(edges)){
    fc_n1 <- fcenv[[as.character(edges$n1[i])]]
    fc_n2 <- fcenv[[as.character(edges$n2[i])]]
    cd <- edges$cor_direction[i]
    
    if(!is.null(fc_n1) & !is.null(fc_n2)){  
      if(fc_n1 == fc_n2 & cd == 1 | fc_n1 != fc_n2 & cd == -1){
        puc_vector[i] <- 1
      } else if(fc_n1 == fc_n2 & cd == -1 | fc_n1 != fc_n2 & cd == 1){
        puc_vector[i] <- 0
      } else {
        print("Something is whrong with PUC calc loop #2")
        cat("fcn1: ", fc_n1, "\n")
        cat("fcn2: ", fc_n2, "\n")
        cat("cd: ", cd, "\n")}
    } else if(is.null(fc_n1) | is.null(fc_n2)){
      puc_vector[i] <- NA
    } else {
      print("Something is whrong with PUC calc loop #1")
      cat("fcn1: ", fc_n1, "\n")
      cat("fcn2: ", fc_n2, "\n")
      cat("cd: ", cd, "\n")}
  }
  return(puc_vector)
}

# FUNCTION: splits edge name into two columns with each node names
# input: edge table with edge name column "pairName"
# output: same table with extra two columns: "n1" with node one name, "n2" -
# node two name
split_edge_names <- function(table){
  genes_L <- strsplit(table[,which(names(test_table) == "pairName")], 
                      split = "<==>", fixed = T)
  table$n2 <- unlist(lapply(genes_L, `[[`, 2))
  table$n1 <- unlist(lapply(genes_L, `[[`, 1))
  return(table)
}