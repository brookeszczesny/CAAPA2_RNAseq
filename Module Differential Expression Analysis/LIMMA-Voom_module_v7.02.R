
rm(list=ls())

library(edgeR)
library(limma)
library(sva)
library(qvalue)
library(dplyr)
library(tibble)
library(stringr)

### Set File Paths

# Active asthma Case Definition used
subset_name <- "Run4_module"

output_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/"

setwd(output_path)
#sink(paste(c(output_path,subset_name,"LIMMA-voom_Out.txt"),collapse=""))

# using full dataset with full module expression not excluding samples
pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"
#pheno_withPC_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/wgcna_pheno.csv"

# module expression
module_voom_counts_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/Run4\ _module_expression.csv"



#-------------------------------------------------------------------------------------------------------
## Read in raw counts and  pheno with PCs

colData=read.csv(pheno_withPC_path,header=T,row.names=1)
colData <- colData[,1:14]
head(colData)
dim(colData)

module_voom_counts=read.csv(module_voom_counts_path,header=T,row.names=1)
names(module_voom_counts) <- sub("^X", "", names(module_voom_counts))  
module_voom_counts[1:5,1:5]
#module_voom_counts <- round(module_voom_counts,0)
#module_voom_counts[1:5,1:5]
countData <- module_voom_counts[,-1]# drop module labels column
countData[1:5,1:5]
dim(countData)[2] == dim(colData)[1] # check sample number matches for pheno and count data


colData$asthmaCaseControl <- 0
colData$asthmaCaseControl[colData$Asthma == "Case"] <- 2
colData$asthmaCaseControl[colData$Asthma == "Control"] <- 1
colData$asthmaCaseControl <- as.factor(colData$asthmaCaseControl)
head(colData)

all.equal(colnames(countData), rownames(colData))


mod1 <- model.matrix(~ colData$asthmaCaseControl+as.factor(colData$Library_Prep_Batch)+colData$Age+as.factor(colData$Gender)+colData$Agilent_RINe+colData$GC+as.factor(colData$site)+colData$PC1+colData$PC2)

#### Run LIMMA
fit <- lmFit(countData, mod1);
fit.eb <- eBayes(fit)
save(fit.eb, file=paste(c(subset_name, "LIMMA-voom_module_out",".rda"),collapse=""))

# make results table
res <- topTable(fit.eb, coef="colData$asthmaCaseControl2",confint=TRUE,number=Inf,adjust.method="BH")

# too few p-values to calculate q so we will use BH correction

results_name <- paste(c(output_path,subset_name,"_module_DE_results.csv"),collapse="")
write.table(res,results_name,sep=",",row.names=T,col.names=NA,quote=F)

#----------------------------------------------------------------------------------
