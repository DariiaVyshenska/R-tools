setwd("H:/Vyshenska/Rwd")

# this function takes .gtf file as input, gives .txt comma separated file as an output
# output file contains repetitions. needs additional cleaning to contain only unique occurrences  
processFile <-  function(filepath, outputfilename) {
  
  # open file to write output to and the file to analyse
  sink(outputfilename)
  con = file(filepath, "r")
  
  # looping through the lines in the genome file
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # gets gene Ids and gene names, cleans them from not needed symbols
    line_split <- unlist(strsplit(line, split = ";"))
    gene_id_raw <- line_split[grep("gene_id", line_split)]
    gene_name_raw <- line_split[grep("gene_name", line_split)]
    gene_id <- gsub(' gene_id \"|\"', '',gene_id_raw)
    gene_name <- gsub(' gene_name \"|\"', '',gene_name_raw)
    new_line <- paste(gene_id, gene_name, sep = ",")
    cat(new_line, "\n")
  }
  
  close(con)
  sink()

}
##########

# get human annotation file
processFile("./genes.gtf", "human_annotations.txt")
