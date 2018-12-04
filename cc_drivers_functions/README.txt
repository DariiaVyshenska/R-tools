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
# corr coef directionality ("coef_dir"). If any NA's are detected - the row 
# will be excluded. uses function extract_table
# fisher_calc
# Function: calculates Fisher's combined probability test p-value and
# attaches this p-value as extra column "fish_pval". Returns table with attached fisher column
# uses function subtable_table
