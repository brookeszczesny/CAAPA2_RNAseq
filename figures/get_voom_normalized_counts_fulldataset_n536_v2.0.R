

# Extract and write out voom normalized count data (log2CPM) to a csv file with sids and counts for a list of user input genes

# takes a user input filename (output_filename) and user input list of genes (target_genes) and outputs "ouput_filename_voom_normalized_counts_Fulldataset_n536.csv"
# to the current working directory
# genes must be input by user in single line space separated list in the format (Ensembl_ID-genename1 Ensembl_ID-genename2 Ensemble_ID-genename3 etc)
# genes not tested in the full dataset will not appear in the output file

# These voom normalized counts are extracted from the final full group n536 Case/Control Differential Expression analysis 

# To run on skyline 
# module load r
# module load r-dplyr
# R
# source("get_voom_normalized_counts_fulldataset_n536_v2.0.R")
# follow the user input prompts
# -------------------------------------------------------------------------------------------------

library(dplyr)
#library(tibble)


normalized_counts_path <- "/data/lad_gphs/analysis/CAAPA/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_voom_norm_counts.csv"

# Input list of genes example
#ENSG00000096060-FKBP5 ENSG00000115414-FN1 ENSG00000113580-NR3C1 ENSG00000080824-HSP90AA1 ENSG00000096384-HSP90AB1 ENSG00000081913-PHLPP1 ENSG00000110958-PTGES3 ENSG00000170606-HSPA4 ENG00000082175-PGR ENSG00000004478-FKBP4 ENSG00000142208-AKT1 ENSG00000198793-MTOR   

#output_filename <- readline(prompt="Enter a name for your output file:")
cat("Enter a name for your output file:")
output_filename <- readLines(con="stdin", n=1)

target_genes <- readline(prompt="Enter a space separated list target genes :")
target_genes <- strsplit(target_genes,"\\s+")[[1]]

# Normalized counts as input
normalized_counts=read.csv(normalized_counts_path,header=T,row.names=1)
names(normalized_counts) <- sub("^X", "", names(normalized_counts))

# Only select target genes from normalized counts dataframe
data_keep_rows <- target_genes
normalized_counts_subset <- normalized_counts[rownames(normalized_counts) %in% data_keep_rows, ]

# Transpose normalized counts to combine with phenotype file. Now, first row is sample ID and column is gene name.
normalized_tr<-t(normalized_counts_subset)
output_path <- paste(c(output_filename[[1]],"_voom_normalized_counts_Fulldataset_n536.csv"),collapse="")
write.csv(normalized_tr,file=output_path ,row.names=TRUE) # Write out counts for file






	
