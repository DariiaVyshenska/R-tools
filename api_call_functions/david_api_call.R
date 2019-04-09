## this portion is required if you need to install RDAVIDWebService
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RDAVIDWebService", version = "3.8")
# 
# browseVignettes("RDAVIDWebService")
## Be aware that you may need to update/(re)install your Java

# Funcion:
# makes a single call to david for a given vector of genes (test_list).
# requires to provide file names for cluster file (clusterFileName) and
# annotation table file (annotFileName). Also requires registered email of
# the user.
# uses ENSEMBL gene IDs and pre-set annotation categories:
# "GOTERM_BP_3", "GOTERM_BP_4", "GOTERM_BP_5", "GOTERM_BP_DIRECT",
# "GOTERM_CC_3", "GOTERM_CC_4", "GOTERM_CC_5", "GOTERM_CC_DIRECT",
# "GOTERM_MF_3", "GOTERM_MF_4", "GOTERM_MF_5", "GOTERM_MF_DIRECT",
# "COG_ONTOLOGY", "UP_KEYWORDS", "UP_SEQ_FEATURE", "OMIM_DISEASE",
# "BBID", "BIOCARTA", "KEGG_PATHWAY", "REACTOME_PATHWAY"


david_call <- function(test_list, clusterFileName, annotFileName, email){
  library("RDAVIDWebService")
  david_url <- "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
  david<-DAVIDWebService(email=email, url=david_url)
  setAnnotationCategories(david, c("GOTERM_BP_3", "GOTERM_BP_4", "GOTERM_BP_5", "GOTERM_BP_DIRECT",
                                   "GOTERM_CC_3", "GOTERM_CC_4", "GOTERM_CC_5", "GOTERM_CC_DIRECT",
                                   "GOTERM_MF_3", "GOTERM_MF_4", "GOTERM_MF_5", "GOTERM_MF_DIRECT",
                                   "COG_ONTOLOGY", "UP_KEYWORDS", "UP_SEQ_FEATURE", "OMIM_DISEASE",
                                   "BBID", "BIOCARTA", "KEGG_PATHWAY", "REACTOME_PATHWAY"))
  
  result<-addList(david, test_list, idType="ENSEMBL_GENE_ID", listName="test", listType="Gene")
  termCluster<-getClusterReport(david, type="Term")
  getClusterReportFile(david, type="Term",fileName=clusterFileName)
  getFunctionalAnnotationChartFile(david, annotFileName)
}