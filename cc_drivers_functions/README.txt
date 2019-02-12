kd_data_manip_functions
#contains several functions for KD project
#
#import_df_asList
# FUNCTION: importing data
# function takes a path to a folder that contains tab delimited files
# files are imported as data.frames and all of them put into a list 
# each table in the list has a name of the file name from wich it was imported
#
#degs_df.R
# This function takes list of tables with info about DEGs and pvalue
# threshold, it spits out dataframe with the name of the KD-target in first 
# column and number of genes that pass this pvalue threshold for that gene
# as a second column
#
#degs_union_vec
# FUNCTION: takes list of tables with unfiltered degs, outputs a vector
# of union of all degs for all KDs under certain pvalue threshold
# also exports degs for this pval for each kd file and union as a file
#
#fc_table_extr
# FUNCTION: takes list of tables with unfiltered degs, outputs a data frame
# with DEGs' fold changes (gene considered to be DEg under specified pvalue 
# threshold, otherwise it will be NA in a table) for each KD 

ratioStat.R
# this funciton takes matrix of ratios and creates two tables:
# table with all combination for each KD target and ratios for these combinations
# table with statistics - min, max, and max/min ratios for each KD-target
# function returns only statistics, but exports both tables as .csv files
# you also need to provid pval under which the ratios matrix was created
# this pval will be used in file names only.



subtable_table.R
# contains three functions: extract_table and subtable_table, fisher_calc. Subtable_table uses extract_table
#
# extract_table
# Function: extracts rows based on key word (value variable) - use "Coefficient"
# or "pvalue"
# subtable_table
# Function: filters table either based on p-value ("p_val") threshold or on 
# corr coef directionality; adds coefition direction column where 1
# is positive, -1 is negative correlation ("coef_dir"); adds median of correlation
# coefitient("corcoef_median") - in the mode of "coef_dir" (encountering any NA leads
# to row exclusion).
# in the mode of "p_val" it only filters the table by individual pvalues
#
# fisher_calc
# Function: calculates Fisher's combined probability test p-value and
# attaches this p-value as extra column "fish_pval". If parameter "fdr" = T
# than fdr for fisher will be calcualted and added as "fisher_fdr" column.
# Returns table
#
# merge_raw_corr
# Function: imports and merges separate correlation analysis results files 
# into one big table (plus exporting it as .csv file as a backup; a must for
# large analysis files)
# Input: a path to the folder with ONLY correlation analysis files (Richard's
# script output; one file = one analysis) - example 
# input_path = "D:/GitHub/R_wd/corr_raw_input/"
# Output: single merged table with merged analysis, analysis headers renamed
# according to the analysis # (from analysis 1 for all to analysis 1, 2, 3 ect.)
# Function also exports the merged table as .csv (mandatory)
#
#puc
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

nw_calc_stats.R
# calc_stats
# Function: partially filters merged raw correlation table and calculates
# key network statistics (no PUC is implemented yet!)
# input: raw single merged correltation table (can and should contain several
# cohorts for metaanalysis; MUST NOT contain groups that should be analysied 
# separately)
# output: .csv file with table with all the data, filtered by correlation 
# directionality and individual p-values (for each cohort); .csv file with 
# table with key network statistics (function also returns this function)
#
# puc_stat
# Function: filters correlation directionality filtered and fisher+fdr 
# calculated table based on individual p-value, fisher, and fisher fdr threshold
# plus gives the statistics on edges according to PUC pass/no pass/not applicable
# and on nodes - ups and downs, fold change consistency regulation across drivers
#
# input: individual p-value(ip), fisher combined p-value(fp), fisher fdr (f_fdr),
# and table of the network that contains columns:
# pairName - gene pair name
# n1 - gene 1 name from the gene pair
# n2 - gene 2 name from the gene pair
# set of columns for each cohort that must be included in metaanalysis (from 
# Richard's script output)
# cor_direction - corelation direction for the gene pair (1/-1)
# fish_pval - for the pair(edge)
# fish_fdr - for the pair(edge)
# PUC - for the pair(edge), can be 1/0/NA
# n1_fc_consist and n2_fc_consist - fold change consistency for gene 1 and 
# gene 2 respectfully, can be 1 (positive cor), -1 (negative cor), 0 (inconsis-
# tent across drivers)
# n1_regulating_drivers and n2_regulating_drivers - number of drivers by which
# gene 1 or gene 2 respectfully regulated (any number >= 1)
#
# output: 
# file with full but filtered with ip, fp, f_fdr table;
# file with staistics about the filtered network
# returns: statistics data frame
