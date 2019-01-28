# Function: partially filters merged raw correlation table and calculates
# key network statistics (no PUC is implemented yet!)
# input: raw single merged correltation table (can and should contain several
# cohorts for metaanalysis; MUST NOT contain groups that should be analysied 
# separately)
# output: .csv file with table with all the data, filtered by correlation 
# directionality and individual p-values (for each cohort); .csv file with 
# table with key network statistics (function also returns this function)


calc_stats <- function(merged_table, indivPval = 0.3, fishPval = 0.05, 
                       fdr_threshold = 0.1){
  require(combinat)
  source("D:/GitHub/R_scripts/cc_drivers_functions/subtable_table.R")
  ##############################################################################
  # SUPPLEMENTARY SMALL FUNCITONS
  # function that gives a matrix of all possible combinations without repetition
  expand.grid.unique <- function(x, y, include.equals=FALSE)
  {
    x <- unique(x)
    
    y <- unique(y)
    
    g <- function(i)
    {
      z <- setdiff(y, x[seq_len(i-include.equals)])
      
      if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    
    do.call(rbind, lapply(seq_along(x), g))
  }
  #
  unique_nodes <- function(intable){
    genes_L <- strsplit(intable[,1], split = "<==>", fixed = T)
    nodes_number <- length(unique(c(lapply(genes_L, `[[`, 1),lapply(genes_L, `[[`, 2))))
    return(nodes_number)
  }
  ################################################################################
  m_table_cd <- subtable_table(merged_table, "coef_dir")
  m_table_indpval <- subtable_table(m_table_cd, "p_val", 
                                    p_val_threshold = indivPval)
  m_table_cd_indpval_fish <- fisher_calc(m_table_indpval, fdr = T)
  
  filt_table_name <- paste("analysisTable_CD_IP", indivPval, ".csv", sep = "")
  write.csv(m_table_cd_indpval_fish, filt_table_name, 
            row.names = F)
  
  col_num <- ncol(m_table_cd)
  short_table <- m_table_cd[,c(1, (col_num-3):col_num)]
  
  # calculate initial number of nodes we started with
  starting_nodes <- unique_nodes(merged_table)
  
  # calculate number of theoretical edges
  theor_edges <- nrow(expand.grid.unique(c(1:starting_nodes),c(1:starting_nodes), 
                                         include.equals = F))
  # calculate number of nodes and edges that pass correlation directionality
  # filter
  cd_nodes <- unique_nodes(m_table_cd)
  cd_edges <- nrow(m_table_cd[as.numeric(m_table_cd$cor_direction) != 0,])
  
  # calculate number of nodes and edges that pass correlation directionality 
  # filter + pass individual pvalue filter
  cd_ip_nodes <- unique_nodes(m_table_indpval)
  cd_ip_edges <- nrow(m_table_indpval[m_table_indpval$cor_direction != 0,])
  
  # calculate number of nodes and edges that pass correlation directionality 
  # filter + pass individual pvalue filter + fisher p-value filter
  m_table_cd_indpval_F <- m_table_cd_indpval_fish[as.numeric(m_table_cd_indpval_fish$fish_pval) <= fishPval,]
  cd_ip_f_nodes <- unique_nodes(m_table_cd_indpval_F)
  cd_ip_f_edges <- nrow(m_table_cd_indpval_F[as.numeric(m_table_cd_indpval_F$cor_direction) != 0,])
  
  # calculate number of nodes and edges that pass correlation directionality 
  # filter + pass individual pvalue filter + fisher p-value filter + 
  # fisher FDR filter
  m_table_cd_indpval_F_FDR <- m_table_cd_indpval_F[as.numeric(m_table_cd_indpval_F$fish_fdr) <= fdr_threshold,]
  cd_ip_f_fdr_nodes <- unique_nodes(m_table_cd_indpval_F_FDR)
  cd_ip_f_fdr_edges <- nrow(m_table_cd_indpval_F_FDR[as.numeric(m_table_cd_indpval_F_FDR$cor_direction) != 0,])
  
  # ratio of filtered edges to filtered nodes
  ratio_edgTonod <- as.numeric(cd_ip_f_fdr_edges)/as.numeric(cd_ip_f_fdr_nodes)
  
  # for correlated pairs that pass:
  # cor dir, individ pval, fisher, fisher fdr filters we calculate:
  # number of positive correlations
  corr_vec <- m_table_cd_indpval_F_FDR$corcoef_median
  
  pos_cor <- sum(sign(as.numeric(corr_vec)) == 1)
  # number of negative correlations
  neg_cor <- sum(sign(as.numeric(corr_vec)) == -1)
  # ratio of number of positive over number of negative correlations
  ratio_posToneg <- as.numeric(pos_cor)/as.numeric(neg_cor)
  # maximum absolute median positive correlation
  max_pos <- abs(max(corr_vec[sign(as.numeric(corr_vec)) == 1]))
  # maximum absolute median negative correlation
  max_neg <- abs(min(corr_vec[sign(as.numeric(corr_vec)) == -1]))
  
  
  ######
  # Creating table with all statistics & exporting it into a file
  feat_val <- c(starting_nodes, theor_edges, cd_nodes, cd_edges, cd_ip_nodes, cd_ip_edges, 
                cd_ip_f_nodes, cd_ip_f_edges, cd_ip_f_fdr_nodes, cd_ip_f_fdr_edges, 
                ratio_edgTonod, pos_cor, neg_cor, ratio_posToneg, max_pos, max_neg)
  features <- c("# of starting nodes", "# of theoretical edges",  "# of cor dir nodes",
                "# of cor dir edges",  "# of cor dir, ip nodes",  "# of cor dir, ip edges",
                "# of cor dir, ip, fish nodes",  "# of cor dir, ip, fish edges",
                "# of cor dir, ip, fish, fish_fdr nodes",
                "# of cor dir, ip, fish,fish_fdr edges",  "ratio of # edges to # nodes",
                "# of positive correlations",  "# of negative correlations",
                "ratio of positive to negative corr",  "max absolute positive corr",
                "max absolute negative corr"
  )
  
  stat_table <- data.frame(FEATURES = features, FEATURE_VALUE = feat_val)
  
  stat_fileName <- paste("StatNW_CD_ip", indivPval, "_fishPval", fishPval, "_FDR",
                         fdr_threshold, ".csv", sep = "")
  write.csv(stat_table, stat_fileName, row.names = F)
  return(stat_table)
}
