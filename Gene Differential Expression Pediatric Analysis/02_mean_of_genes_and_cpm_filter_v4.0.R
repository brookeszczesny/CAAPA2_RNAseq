
rm(list=ls())

library(dplyr)
library(tidyr)
library(sva)
library(DESeq2)
library(tibble)


### Set File Paths
# Active asthma case definition used

#subset_name <- "full_n536_"
#subset_name <- "Adult_n376_"
subset_name<- "Ped_n160_"

#output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/"
#output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/adult/"
output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/ped/"

setwd(output_path)
sink(paste(c(output_path,subset_name,"02_mean_of_genes_and_cpm_filter_Out.txt"),collapse=""))

#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Adult_n376_withPC.csv"
pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Ped_n160_withPC.csv"

cpm_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/merged_read_counts/60664_genes/merged_CPM_CoCo_counts.txt"
raw_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/merged_read_counts/60664_genes/merged_raw_CoCo_counts.txt"

#Set cpm threshold; percentage corresponding to freq of smaller group in cases/controls
#cpm_threshold_value <- 0.4720 #full
#cpm_threshold_value <- 0.4628 #Adult
cpm_threshold_value <- 0.4938 #Pediatric

# Set mean threshold; used to select for genes with Deseq2 normalized count mean greater than or equal to 20
mean_thresh <- 20

#-------------------------------------------------------------------------------------------------

### Read in pheno and merged raw counts from CoCo
# note that this does not use the pcs
pheno = read.csv(pheno_withPC_path, header=TRUE, row.names=1)
raw = read.table(raw_path,header=T,row.names=1,sep=",")
names(raw) <- sub("^X", "", names(raw))

### Read in genes merged CPM counts
# Ensembl gene name rows, sample ID columns
cpm = read.table(cpm_path,header=T,row.names=1,sep=",")
names(cpm) <- sub("^X", "", names(cpm))

### Subset CPM counts to pheno samples (Adult or Ped).
# cpm_final contains sample ID columns, Ensembl gene name rows
final_samples=row.names(pheno)
cpm_final= cpm[,c(t(final_samples))]
write.csv(cpm_final,paste(c(output_path,subset_name,"cpm_counts.csv"),collapse=""),quote=F,row.names=T)

### Filter 1: CPM > 0 for cpm_threshold_value of samples### Run filter
# set n equal to the number of columns (sample IDs) in cpm_final
# create new column called Threshold; Applying across the row (for each gene)
# find how many times the CPM is > 0; if the CPM is greater than 0 for cpm_threshold_value or more samples return TRUE
# write a csv with Ensemble gene names rows and Threshold value TRUE/False Column)
n <- dim(cpm_final)[2]
cpm_final$Threshold = ifelse(apply(cpm_final, 1, function(x) length(x[x>0])>= cpm_threshold_value * n), TRUE, FALSE)
cpm_thresh<-cpm_final %>%
   select(Threshold)
   
write.csv(cpm_thresh,paste(c(output_path,subset_name,"_threshold", cpm_threshold_value,"filtered_cpm_counts.csv"),collapse=""),row.names=T,quote=F)

### Subset raw counts to pheno samples.
final_samples=row.names(pheno)
raw_final= raw[,c(t(final_samples))]
write.csv(raw_final,paste(c(output_path,subset_name,"_raw_counts.csv"),collapse=""),quote=F,row.names=T)

### DESeq2
colData=read.csv(pheno_withPC_path,header=T,row.names=1) #fix this to be full pheno  with pc subset
countData=read.csv(paste(c(output_path,subset_name,"_raw_counts.csv"),collapse=""),header=T,row.names=1)
names(countData) <- sub("^X", "", names(countData))
countData[1:5,1:5]
countData <- round(countData,0)
countData[1:5,1:5]

head(colData)
colData$asthmaCaseControl <- 0
colData$asthmaCaseControl[colData$Asthma == "Case"] <- 2
colData$asthmaCaseControl[colData$Asthma == "Control"] <- 1
colData$asthmaCaseControl <- as.factor(colData$asthmaCaseControl)

all.equal(colnames(countData), rownames(colData))

### DESeq2 design + normalized counts

cds <- DESeqDataSetFromMatrix(countData = countData, colData = colData,design = ~ asthmaCaseControl+as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+as.factor(site)+PC1+PC2)
dds <- estimateSizeFactors(cds)
normalized_counts <- counts(dds, normalized=TRUE)

write.csv(normalized_counts,paste(c(output_path,subset_name,"deseq2_normalized_counts.csv"),collapse=""),quote=F,row.names=T)

### Get mean of normalized counts for each gene
normalized_counts<-read.csv(paste(c(output_path,subset_name,"deseq2_normalized_counts.csv"),collapse=""),header=T,row.names=1)
names(normalized_counts) <- sub("^X", "", names(normalized_counts))
normalized_counts$mean <- rowMeans(normalized_counts)
normalized_counts_mean <- normalized_counts %>%
   select(mean)
   
write.csv(normalized_counts_mean,paste(c(output_path,subset_name,"deseq2_normalized_counts_mean.csv"),collapse=""),quote=F,row.names=T)

###Select for genes passing cpm threshold and mean threshold
cpm_thresh_pass<- cpm_thresh %>%
	rownames_to_column("gene_ID")%>%
	filter(Threshold=="TRUE")
	
normalized_counts_mean_pass<- normalized_counts_mean %>%
	rownames_to_column("gene_ID")%>%
	filter(mean>=mean_thresh)

passing_filtered_genes <- merge(cpm_thresh_pass,normalized_counts_mean_pass,by="gene_ID")
n_passing_filtered_genes <-nrow(passing_filtered_genes)
print(paste(c(n_passing_filtered_genes, " genes pass the filter!"),collapse=""))

write.csv(passing_filtered_genes,paste(c(output_path,subset_name,"_", n_passing_filtered_genes,"_filtered_genes.csv"),collapse=""),quote=F,row.names=F)

sessionInfo()
savehistory(file=paste(c(output_path,subset_name,"02_mean_of_genes_and_cpm_filter.Rhistory"),collapse=""))
sink()
closeAllConnections()