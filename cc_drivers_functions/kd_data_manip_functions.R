require(gtools)
# FUNCTION: importing data
# function takes a path to a folder that contains tab delimited files
# files are imported as data.frames and all of them put into a list 
# each table in the list has a name of the file name from wich it was imported
import_df_asList <- function(path){
  list_of_tables <- list()
  for (file in list.files(path)){
    cat("Importing file: ", file, "\n", sep = "")
    file_path <- paste(input_folder, file, sep = "/")
    table <- read.table(file = file_path, sep = '\t', header = TRUE, 
                        check.names = F, stringsAsFactors = F) 
    file_name <- gsub(".txt", "", file)
    list_of_tables[[file_name]] <- table
  }
  return(list_of_tables)
}

# FUNCTION: takes list of tables with info about DEGs and pvalue
# threshold, it spits out dataframe with the name of the KD-target in first 
# column and number of genes that pass this pvalue threshold for that gene
# as a second column
degs_df <- function(tables, pval){
  deg_vec <- vector(length =  length(tables))
  for(i in 1:length(tables)){
    print(names(tables[i]))
    deg_vec[i] <- sum(as.numeric(tables[[i]][["Parametric p-value"]]) <= pval)
  }
  df <- data.frame("KDtarget" = names(tables), "DEGs" = deg_vec)
  return(df)
}

# FUNCTION: takes list of tables with unfiltered degs, outputs a vector
# of union of all degs for all KDs under certain pvalue threshold
# also exports degs for this pval for each kd file and union as a file
degs_union_vec <- function(tables, pval){
  # selecting gene ID's that pass pval threshold into a list
  deg_list <- list()
  for(i in 1:length(tables)){
    gene_name <- names(tables[i])
    sel_vec <- as.numeric(tables[[i]][["Parametric p-value"]]) <= pval
    id_vec <- tables[[i]][["UniqueID"]]
    deg_list[[gene_name]] <- id_vec[sel_vec]
  }
  
  # exporting the list containing all degs under this pval into a file
  file_nameAllDegs <- paste("degs_per_KD_pval", pval, ".txt", sep = "")
  sink(file_nameAllDegs)
  print(deg_list)
  sink()
  
  # getting the union of all degs under this pval and exporting it into a file
  # as well as returing it as the function output
  file_nameDegsUnion <- paste("degs_union_pval", pval, ".txt", sep = "")
  uni_id <- unique(unlist(deg_list))
  write.table(uni_id, file_nameDegsUnion, col.names = F, row.names = F)
  return(uni_id)
}

