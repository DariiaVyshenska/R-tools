This folder contains modules and functions for calling different database APIs using R.

david_api_call.R
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