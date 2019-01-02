ratioStat.R
# this funciton takes matrix of ratios and creates two tables:
# table with all combination for each KD target and ratios for these combinations
# table with statistics - min, max, and max/min ratios for each KD-target
# function returns only statistics, but exports both tables as .csv files
# you also need to provid pval under which the ratios matrix was created
# this pval will be used in file names only.

degs_df.R
# This function takes list of tables with info about DEGs and pvalue
# threshold, it spits out dataframe with the name of the KD-target in first 
# column and number of genes that pass this pvalue threshold for that gene
# as a second column

subtable_table.R
# contains three functions: extract_table and subtable_table, fisher_calc. Subtable_table uses extract_table
# extract_table
# Function: extracts rows based on key word (value variable) - use "Coefficient"
# or "pvalue"
# subtable_table
# Function: filters table either based on p-value ("p_val") threshold or on 
# corr coef directionality; adds coefition direction column where 1
# is positive, -1 is negative correlation ("coef_dir"); adds median of correlation
# coefitient("corcoef_median"). 
# If any NA is detected - row is excluded
# fisher_calc
# Function: calculates Fisher's combined probability test p-value and
# attaches this p-value as extra column "fish_pval". If parameter "fdr" = T
# than fdr for fisher will be calcualted and added as "fisher_fdr" column.
# Returns table
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