
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
#subset_name <- "full_n536_"
#subset_name <- "Adult_n376_"
#subset_name<- "Ped_n160_"

#subset_name <- "Nigeria_n66_"
#subset_name <- "Baltimore_n76_"
#subset_name <- "Barbados_n81_"
subset_name <- "Brazil_n80_"
#subset_name <- "Chicago_n84_"
#subset_name <- "Denver_n87_"
#subset_name <-"Washington_DC_n62_"

#output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/"
#output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/adult/"
#output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/ped/"
output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/"

setwd(output_path)
sink(paste(c(output_path,subset_name,"LIMMA-voom_Out.txt"),collapse=""))

# confirm if analysis is by site
#by_site <- FALSE
by_site <- TRUE

#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Adult_n376_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Ped_n160_withPC.csv"

#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Nigeria_n66_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Baltimore_n76_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Barbados_n81_withPC.csv"
pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Brazil_n80_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Chicago_n84_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Denver_n87_withPC.csv"
#pheno_withPC_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Washington_DC_n62_withPC.csv"

# Active asthma
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/adult/Adult_n376__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/ped/Ped_n160__raw_counts.csv"

#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Nigeria_n66__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Baltimore_n76__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Barbados_n81__raw_counts.csv"
raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Brazil_n80__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Chicago_n84__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Denver_n87__raw_counts.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/Washingtond_DC_n62__raw_counts.csv"

passed_genes_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_filtered_genes.csv"
#passed_genes_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/adult/Adult_n376__21789_filtered_genes.csv"
#passed_genes_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/adult_vs_ped/ped/Ped_n160__21887_filtered_genes.csv"


#-------------------------------------------------------------------------------------------------------

### Further subset raw counts to genes that pass the filter
raw_counts = read.csv(raw_counts_path,header=T,row.names=1)
names(raw_counts) <- sub("^X", "", names(raw_counts))
passed_genes = read.csv(passed_genes_path,header=T,row.names=1)

filtered_raw_counts <- raw_counts %>%
	rownames_to_column("gene_ID")%>%
	filter(gene_ID %in% rownames(passed_genes))

n_passing_filtered_genes <- dim(filtered_raw_counts)[1]
filtered_raw_counts_path <- paste(c(output_path,subset_name,"_", n_passing_filtered_genes,"_filtered_raw_counts.csv"),collapse="")
write.csv(filtered_raw_counts,filtered_raw_counts_path,quote=F,row.names=F)

### Read in raw counts and  pheno with PCs
colData=read.csv(pheno_withPC_path,header=T,row.names=1)
countData=read.csv(filtered_raw_counts_path,header=T,row.names=1)
n_passing_filtered_genes <- dim(countData)[1]
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

### Add asthma case/control as group to create a DGE oblect for edgeR

DGE <- DGEList(counts=countData,group=colData$asthmaCaseControl)

#### Normalize using edgeR

y <- calcNormFactors(DGE);

#### voom + SVA

# if analysis is not by site, include site as a factor
# if analysis is by site, do not include site as factor
if (by_site==FALSE) {
	mod <- model.matrix(~ asthmaCaseControl+as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+as.factor(site)+PC1+PC2,data=colData)
	mod0 <-  model.matrix(~ as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+as.factor(site)+PC1+PC2,data=colData)
} else {
	mod <- model.matrix(~ asthmaCaseControl+as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+PC1+PC2,data=colData)
	mod0 <-  model.matrix(~ as.factor(Library_Prep_Batch)+Age+as.factor(Gender)+Agilent_RINe+GC+PC1+PC2,data=colData)
}

v <- voom(counts=y, design=mod,plot=TRUE)

voom_norm_counts <- v$E
write.table(voom_norm_counts,paste(c(output_path,subset_name,"_", n_passing_filtered_genes,"_voom_norm_counts.csv"),collapse=""),sep=",",row.names=T,quote=F)


svseq <- sva(v$E, mod, mod0);

#### Add SVs to the model and rerun voom
if (by_site==FALSE) {
	mod1 <- model.matrix(~ colData$asthmaCaseControl+as.factor(colData$Library_Prep_Batch)+colData$Age+as.factor(colData$Gender)+colData$Agilent_RINe+colData$GC+as.factor(colData$site)+colData$PC1+colData$PC2+svseq$sv)
} else {
	mod1 <- model.matrix(~ colData$asthmaCaseControl+as.factor(colData$Library_Prep_Batch)+colData$Age+as.factor(colData$Gender)+colData$Agilent_RINe+colData$GC+colData$PC1+colData$PC2+svseq$sv)
}

v <- voom(counts=y, design = mod1)

#### Run LIMMA
fit <- lmFit(v, mod1);
fit.eb <- eBayes(fit)
save(fit.eb, file=paste(c("LIMMA-voom_out",subset_name,".rda"),collapse=""))

#-------------------------------------------------------------------------------------------------------
# make results table
res <- topTable(fit.eb, coef="colData$asthmaCaseControl2",confint=TRUE,number=Inf)

res$SE <- res$logFC/res$t # std error = effect size / t-stat

p.value.asthma <- fit.eb$p.value[,"colData$asthmaCaseControl2"] # matches p in toptable; used to get correct order of q-values
q.value.asthma <- qvalue(p.value.asthma)$q # calculate q-value
pval_qval<-cbind(p.value.asthma,q.value.asthma)

res <- merge(res,pval_qval,by='row.names',all=TRUE)

subset_name_label <- substr(subset_name,1,nchar(subset_name)-1)

log2FC_label <-paste(c("log2FC.",subset_name_label),collapse="")
CI.L_label <-paste(c("CI.L.",subset_name_label),collapse="")
CI.R_label <-paste(c("CI.R.",subset_name_label),collapse="")
AveExp_label <- paste(c("AveExp.",subset_name_label),collapse="")
tstat_label <- paste(c("t.",subset_name_label),collapse="")
pvalue_label <-paste(c("p.value.",subset_name_label),collapse="")
res$adj.P.Val <- NULL # remove adjusted p-value
res$B <- NULL # remove log odds
SE_label <- paste(c("SE.",subset_name_label),collapse="")
res$p.value.asthma <- NULL # remove duplicate pvalue column
qvalue_label <- paste(c("q.value.",subset_name_label),collapse="")

# relabel columns in format value_type.subset_name; change rownames to gene_ID column
colnames(res) <- c("gene_ID",log2FC_label,CI.L_label,CI.R_label,AveExp_label,tstat_label,
pvalue_label,SE_label, qvalue_label)

results_name <- paste(c(output_path,subset_name,"_", n_passing_filtered_genes,"_results.csv"),collapse="")
write.table(res,results_name,sep=",",row.names=F,quote=F)

#### QQ plot
library(ggplot2)
library(RVA)

voom<-read.csv(results_name,row.names=1,header=T)
png(filename=paste(c(output_path,subset_name,"_", n_passing_filtered_genes,"_qqplot.png"),collapse=""),width = 768, height = 768)
p0 <- plot_qq(data=voom,p.value.flag=pvalue_label,ci=0.95) + theme_bw(base_size = 24)+annotate(
     geom = "text",
     x = -Inf,
     y = Inf,
     hjust = -0.15,
     vjust = 1 + 0.15 * 3,
     label = sprintf("λ = %.4f",median(qchisq(1-voom[[pvalue_label]],1), na.rm=T) / qchisq(0.5,1)),
     size = 12
   )
p0 
dev.off()

### p-value histogram
png(filename=paste(c(output_path,subset_name,"_",n_passing_filtered_genes,"_pvalue_hist.png"),collapse=""),width = 768, height = 768)
p1<-hist(voom[[pvalue_label]],xlab="P-Value",cex.axis=1.5,cex.lab=1.5,main="")
dev.off()


### volcano plot
# color genes with qvalue < 0.05 by #sites with pvalue < 0.05
png(filename=paste(c(output_path,subset_name,"_",n_passing_filtered_genes,"_volcano_plot.png"),collapse=""),width = 768, height = 768)
p2<-plot(res[[log2FC_label]],-log10(res[[pvalue_label]]),main=paste(c(subset_name,"Volcano Plot"),collapse=" "),xlab="log2FC",ylab="-log10(P-Value)", col = ifelse(res[[qvalue_label]]<0.05,'red','black'), pch = 20 )
p2+ legend(x="topright",legend=c("q-value < 0.05","q-value >= 0.05"), fill=c("red", "black"))
dev.off()



sessionInfo()
savehistory(file=paste(c(output_path,subset_name,"LIMMA-voom.Rhistory"),collapse=""))
#sink()
closeAllConnections()