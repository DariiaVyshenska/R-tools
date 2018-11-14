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

