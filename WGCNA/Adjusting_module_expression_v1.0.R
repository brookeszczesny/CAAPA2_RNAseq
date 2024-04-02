# To correct module expression values for age, site,sex, batch,sex, gc, PCs, etc

#Step1: run a linear regression for each module: Module ~ site + age + sex + PC1 + PC2 + other covariates I am forgetting (do not include asthma)
#Step 2: save residuals from the model1
#Step 3: add the intercept from the model back to the residuals to rescale the values onto the original.
rm(list=ls())

library(tibble)
library(dplyr)
library(tidyr)
library(limma)


output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/module_expression/16FEB2023/"
setwd(output_path)

#module expression 
# mean of batch corrected counts for genes assigned to each module

module_expression_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/Run4\ _module_expression.csv"

pheno_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"

#-----------------------------------------------------------------------------------------------------------------------------------

# Module expression as input (only normalized for library size)
module_expression=read.csv(module_expression_path,header=T,row.names=1)
names(module_expression) <- sub("^X", "", names(module_expression))

#module_key <- module_expression%>% rownames_to_column("module_colors")%>% select(1:2) # module number color key
module_expression <- module_expression[,-1]# drop module labels column
module_expression[1:5,1:5]

# Phenotype file
colData=read.csv(pheno_path,header=T,row.names=1)


mod <- model.matrix(~ as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+as.factor(site)+PC1+PC2,data=colData)
fit <- lmFit(module_expression, mod);

intercepts <- as.data.frame(fit$coefficients[,1])%>% rownames_to_column("module_colors")

residuals_matrix <- as.data.frame(residuals(fit,module_expression)) %>% rownames_to_column("module_colors")
residuals_matrix[1:5,1:5]

int_res <- merge(intercepts, residuals_matrix,by="module_colors", all=TRUE) %>% rename("intercept"="fit$coefficients[, 1]")

module_expression_new <- int_res[,3:538]+ int_res$intercept
module_expression_new <- cbind( int_res[, 1],module_expression_new)%>% rename("module_colors"="int_res[, 1]")
module_expression_new[1:5,1:5]

write.table(module_expression_new,"module_expression_adj_covariates.csv",sep=",",row.names=F)
