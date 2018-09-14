getwd()
setwd("D:/GitHub/R_wd/")

#======================
# FUNCTIONS
#======================


# FUNCTION: transforming a list into a data.frame for exprot:
# input - list of vectors with names;
# output - data.frame with names of vectors in column 1, collapsed comma separated
# values from each vector in column 2
transform_list_out <- function(input_list, annotations){

  out_list <- data.frame(KD = character(), 
                         Genes=character(), 
                         stringsAsFactors=FALSE) 
  for(i in 1:length(input_list)){
    symbols <- unique(annotations$geneID[annotations$EnsemblID %in% input_list[[i]]])
    col_genes <- paste(symbols, collapse = ",")
    out_list[i,] <- c(names(input_list[i]), col_genes)
  }
  return(out_list)
}


# FUNCTION: importing data
# function takes a path to a folder that contains tab delimited files
# files are imported as data.frames and all of them put into a list 
# each table in the list has a name of the file name from wich it was imported
import_df_asList <- function(path){
  list_of_tables <- list()
  for (file in list.files(path)){
    #file <- "TPX2.txt"
    cat("Importing file: ", file, "\n", sep = "")
    file_path <- paste(input_folder, file, sep = "/")
    table <- read.table(file = file_path, 
                        sep = '\t', header = TRUE, check.names = F, stringsAsFactors = F) 
    gene_name <- gsub(".txt", "", file)
    list_of_tables[[gene_name]] <- table
  }
  return(list_of_tables)
}


# FUNCTION: filtering list of tables for PUC complience
# input - list of tables that contain column "PUC-complience?" (0 or 1)
# and filters this list to contain either all genes (puc = NA), 
# puc-complient genes (puc = T), puc not-complient genes (puc = F)
# output - list of tables
check_puc <- function(all_tables_test, puc_threshold = NA){
  if (is.na(puc_threshold)){
    cat("No filtering for PUC is done!\n")
    return(all_tables)
  } else if (puc_threshold == T){
    puc_val <- 1
  } else if (puc_threshold == F){
    puc_val <- 0
  } else if (cat("'puc' argument must be T/F/NA!\n"))
  cat("Filters PUC\n")
  new_tables_list <- list()
  for (table in 1:length(all_tables)){
    table_name <- names(all_tables[table])
    print(table_name)
    condition <- all_tables[[table]][["PUC-complient?"]] == puc_val
    new_tables_list[[table_name]] <- all_tables[[table]][eval(condition),]
  }
  return(new_tables_list)
}


# FUNCTION: filtering tables, creating list of genes and list of their ratios
# input - list of imported tables
# output - matrix with shared gene numbers/normalized ratios
#
# all_tables_lis - list of tables from BRB
# puc = NA/T/F
# pval_threshold = numeric OR F
# fdr_threshold = numeric OR F
# shared_direct = "same", "diff" or "all"
# normalized = F OR T
# export = T/F
# annotation_file_path = file to path containing two columns - EnsemblID and geneID (tab del)
get_matrix <- function (all_tables_list, puc = NA, pval_threshold, 
                        fdr_threshold, shared_direct, 
                        normalized, export = F, annotation_file_path = NA) {
  # filtering for PUC if nesessary
  all_tables_list <- check_puc(all_tables = all_tables_list, puc_threshold = puc)
  
  genes_list <- list()
  ratio_list <- list()
  
  # Extracting sublists based on p-value and FDR thresholds
  for(table in 1:length(all_tables_list)){
    sel_vec_signif <- vector()
    gene_name <- NA
    cat("Working with: ", names(all_tables_list[table]), "\n")
    if(is.numeric(pval_threshold) & is.numeric(fdr_threshold)){
      cat("Filtering by p-value and FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold &+
        all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (is.numeric(pval_threshold) & fdr_threshold == F){
      cat("Filtering only by p-value.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < pval_threshold
    } else if (pval_threshold == F & is.numeric(fdr_threshold)){
      cat("Filtering only by FDR.\n\n")
      sel_vec_signif <- all_tables_list[[table]][["FDR"]] < fdr_threshold
    } else if (pval_threshold == F & fdr_threshold == F){
      sel_vec_signif <- all_tables_list[[table]][["Parametric p-value"]] < 1.1
    } else { cat ("FDR and p-value MUST be either 'F' or numeric!\n")}
    
    gene_name <- names(all_tables_list[table])
    genes_list[[gene_name]] <- all_tables_list[[table]][["UniqueID"]][sel_vec_signif]
    ratio_list[[gene_name]] <- all_tables_list[[table]][["KD/C-all"]][sel_vec_signif]
  }
  
  cat(lengths(genes_list), "\n")
  
  ####
  # Writing out all genes that pass the threshold in an .csv file
  if(export == T) {
    if(!is.na(annotation_file_path)){
      annotation_file <- annotation_file_path
      annotations <- read.table(annotation_file_path, stringsAsFactors = F, 
                                check.names = F, header = T)
    } else {cat("You must provide annotation file!\n")
      break
      }
    

    genes_list_tr <- transform_list_out(genes_list, annotations)
    fullList_file_name <- paste("fullGeneListFDR-", 
                                fdr_threshold, "pval-", pval_threshold,".csv", sep="")
    write.csv(genes_list_tr, fullList_file_name, row.names = F)
  }
  #####
  
  # vector with total numbers of genes controlled by each KD - for normalization
  all_numbers <- lengths(genes_list)
  # collecting KD genes names
  nm_vec <- names(genes_list)
  # creating matrix to populate with numbers of shared genes
  shared_numbers <- matrix(0, nrow = length(genes_list), ncol = length(genes_list))
  rownames(shared_numbers) <- nm_vec
  colnames(shared_numbers) <- nm_vec
  
  # creating a list to populate with the shared genes names - for export in .csv file
  shared_genes <- list()
  
  # populating matrix and list based on directionality and normalization parameters
  for (i in 1:(length(genes_list) - 1)){
#i <- 1
    cat("I'm working with: ", nm_vec[i], "\n", sep = "")
    for (j in (i+1):length(genes_list)){
#j <- 2
      cat("comparing with ", names(genes_list[j]), "\n")
      colhead <- paste(names(genes_list[i]),"vs", names(genes_list[j]), sep = "")
      shared_only <- genes_list[[i]][genes_list[[i]] %in% genes_list[[j]]]
      length(shared_only)
      
      if(shared_direct == "all"){
        if (normalized == F){shared_numbers[i,j] <- length(shared_only)
        } else if (normalized == T){
          shared_numbers[i,j] <- length(shared_only)/mean(c(all_numbers[i], all_numbers[j]))
        } else {cat("ERROR: parameter 'normalized' must be T or F!\n")}
        if(length(shared_only) == 0){
          shared_genes[[colhead]] <- NA
        } else {shared_genes[[colhead]] <- shared_only}
      } else {
        shared_same_dir <- c()
        
        for(gene in shared_only){
          ratio_i <- ratio_list[[i]][genes_list[[i]] == gene]
          ratio_j <- ratio_list[[j]][genes_list[[j]] == gene]
          print(gene)
          print(ratio_i)
          print(ratio_j)

          if(shared_direct == "same"){
            condition <- expression(ratio_i > 1 & ratio_j > 1 | ratio_i < 1 & ratio_j < 1)
          } else if (shared_direct == "diff"){
            condition <- expression(ratio_i < 1 & ratio_j > 1 | ratio_i > 1 & ratio_j < 1)
          }
          #print(condition)
          #print(eval(condition))
          if(eval(condition) == T){
            shared_same_dir <- c(shared_same_dir, gene)
          } else if (eval(condition) == F){print("Different direction!")
            } else {print("ERROR in conditions!")}
        }
        
        #print(shared_same_dir)
        if (normalized == F){shared_numbers[i,j] <- length(shared_same_dir)
        } else if (normalized == T){
          shared_numbers[i,j] <- length(shared_same_dir)/mean(c(all_numbers[i], all_numbers[j]))
        } else {cat("ERROR: parameter 'normalized' must be T or F!\n")}
        if(length(shared_same_dir) == 0){
          shared_genes[[colhead]] <- NA
        } else {shared_genes[[colhead]] <- shared_same_dir}
      }
    } 
  }
  
  if(export == T) {
    # exporting in .csv files matrix with numbers/ratios of shared genes
    # and list of shared genes for each KD
  
    shared_g_df <- transform_list_out(shared_genes, annotations)
    matrix_file_name <- paste("shared-numb-fdr", 
                              fdr_threshold, "KDdir", toupper(shared_direct),".csv", sep="")
    gene_list_file <- paste("shared-genes-fdr", 
                            fdr_threshold, "KDdir", toupper(shared_direct),".csv", sep = "")
    write.csv(shared_numbers, matrix_file_name)
    write.csv(shared_g_df, gene_list_file, row.names = F)
  }
  return(shared_numbers)
}
#======================================================================
# ALL genes data analysis
#======================================================================

# importing all tables as a list (analysis for all genes detected by RNASeq)
input_folder <- "./tables-all-ge-analysis"
all_tables <- import_df_asList(input_folder)

# getting list of matrixes generated with different p-value thresholds
matrix_list <- list()
i <- 1
for (p in seq(0.0005, 0.05, by = 0.0005)){
  matrix_list[[i]] <- get_matrix(all_tables_list = all_tables, 
                                                pval_threshold = p, fdr_threshold = F, 
                                                shared_direct = "same", normalized = T)
  i <- i + 1
}


# extracting ratios for each comparison and plotting it in the file
#########################
# creating data frame to populate
kd_names <- row.names(matrix_list[[1]])
comparis_names <- vector()
for(n in kd_names){
  for(n2 in kd_names){
    comparis_names <- append(comparis_names, paste(n, "vs", n2, sep = ""))
  }
}
df_ratios <- data.frame(matrix(ncol = 100, nrow = 81))
row.names(df_ratios) <- comparis_names
colnames(df_ratios) <- seq(0.0005, 0.05, by = 0.0005)

# populating data frame with ratios
count <- 1
for(matrix_ in matrix_list){
  df_ratios[, count] <- as.vector(t(matrix_))
  count <- count +1
}

# building plots and exporting them to the file
pdf(file='ratio_plots.pdf')
for(r in 1:length(row.names(df_ratios))){
  if (sum(df_ratios[r,]) != 0){
    cat("Working with: ",row.names(df_ratios[r,]), "\n")
    plot(x = colnames(df_ratios), y = df_ratios[r, ], xlab = "p-value threshold",
         ylab = "ratio (overlapped/mean_total#)", main = row.names(df_ratios[r,]))
  }
}
dev.off()
#########################################################################################

#======================================================================
# Nature Communication genes data anaysis
#======================================================================
# importing all tables as a list (analysis for genes detected by RNASeq AND present
# in Nature communication signature

input_folder <- "./test_files"
all_tables_test <- import_df_asList(input_folder)

test_out <- get_matrix(all_tables_list = all_tables_test, puc = T,
                               pval_threshold = 0.01, fdr_threshold = F, 
                               shared_direct = "all", normalized = F, export = T, 
                               annotation_file_path = "./human_annotation_GTF.txt")
