#This code uses count txt files (from 4 replicates per sample) for 20 samples
#it creates and exports 2 tables - one with all counts per all replicates, second - sum table of counts for all samples (summation
#of counts per gene across all replicates, for each sample)


getwd()
setwd("D:/MORGUN LAB/11_TEMP_ANALYSIS/R_wd")

#here code reads all the files in the working directory and it will use it for murging
sam.num <- 20
rep.num <- 4

files.list <-  c(list.files(path="."))
#print(files.list)
#length(files.list)
library(gtools)
#now I sort the file names so that it will be nicely sorted
files.list <- mixedsort(files.list)
#create a table to merge all other tables to it. I aslo attach the header to this table - table1
table1 <-  read.table(files.list[1])
header.list <- c("geneID", files.list[1])
names(table1) <- header.list
#merging the tables
d <- length(files.list)
for (replicate in files.list[2:d]){
  temp.table <- data.frame()
  temp.table <- read.table(replicate)
  names(temp.table) <- c("geneID",replicate)
  table1 <- merge(table1, temp.table, all = T)
}

####################################################

#checking column names to also check if merging was done right on replicates. 
#This particular code is written for the table where first column is geneID, and all others are replicates by 4.
#you expecting nothing to be printed if replicates are placed together. for add. checking - you can delite one "#" and see printing TRUE for all replicates placed together

col.names.check <- colnames(table1)
col.names.check.no1name <- c(col.names.check[2:length(col.names.check)])
library(stringr)
list.col.check <- str_split(col.names.check.no1name, pattern = "_L00")
for(l in c(seq(1,length(list.col.check),by=rep.num))){
  #print (l)
  #print (list.col.check[[(l)]][1])
  #temp.vector <- vector(mode="numeric", length=0)
  for (x in c(0:3)){
    if (list.col.check[[(l)]][1] == list.col.check[[(l+x)]][1]) {
      #print ("TRUE")
    } else {
      print ("FALSE for")
      print (list.col.check[[(l)]][1])
      print (l)
    }
    #temp.vector <- append(temp.vector, list.col.check[[(l+x)]][1])
    #print (temp.vector)
  }}

##########writing the full table with all replicates into excel file

write.csv(table1, "merged_raw.csv", row.names = F)


##########################
#generating table with sum of all replicates per sample (4 repl-s per sample)

mylist <- list() #create an empty list
i <- 1 # variable for counting in the loop

for (m in c(seq(2,length(colnames(table1)),by=rep.num))) {
  vec <- numeric(nrow(table1)) #preallocate a numeric vector
  k <- 0
  k <- m+3
  vec <- rowSums(table1[m:k])
  mylist[[i]] <- vec #push each sum of 4 replicates in to the list
  i <- i+1 
}
df <- do.call("rbind",mylist) #combine all vectors into a matrix
df_transposed <- t(df) #transpose martix to have geneIDs as rows and samples as columns

matrix.col.names <- numeric(sam.num) #creating vector for col. names
sp <- 1#variable for for loop
for(l in c(seq(1,length(list.col.check),by=rep.num))){
  matrix.col.names[sp] <- list.col.check[[l]][1]
  sp <- sp+1
  #print (l)
  #print (list.col.check[[(l)]][1])
  #temp.vector <- vector(mode="numeric", length=0)
} #for loop to generate list of unique sample names

matrix.row.names <- table1$geneID #extracting row names

rownames(df_transposed) <- matrix.row.names #assigning row names to the matrix
colnames(df_transposed) <- matrix.col.names #assigning column names to the matrix

#writing the matrix into excel file
write.csv(df_transposed, "repl_sum_count.csv") 

#######

#counts_mx <- df_transposed[-1:-5,]

