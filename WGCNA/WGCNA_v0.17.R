
# WGCNA
rm(list=ls())
library(edgeR)
library(limma)
library(WGCNA)
library(dplyr)
library(tibble)

options(stringsAsFactors = FALSE);
#enableWGCNAThreads()

normalized_counts_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_voom_norm_counts.csv"
pheno_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"
DE_results_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/summary/combined_results.csv"
#raw_counts_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__raw_counts.csv"

#limma_data_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/LIMMA-voom_out.rda"
voom_counts_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_voom_norm_counts.csv"


output_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/"

setwd(output_path)

subset_name <- "Run4"

#parameters
FDR_threshold <- 0.15 #cuttoff to filter into WGCNA
sample_outlier_cutoff <- 58
min_module_size <- 15
max_block_size <- 80000 #default is 5000; but we have 10GB memory
network_type <- "signed"
tom_type <- "signed"
deep_split <- 3
pam_stage <- TRUE
power_beta<-12
R2_cutoff <- 0.9



#----------------------------------------------------------------
# filtering based on DE results
DE_results <- read.csv(DE_results_path,header=T)

# get genes which had q-value of 0.4 or less
FDR_filtered_genes <- DE_results %>%
	filter(q.value.full_n536 < FDR_threshold)%>%
	select(gene_ID)
	
FDR_filtered_genes <- as.vector(FDR_filtered_genes[,1]) #convert to vector

n_FDR_filtered_genes<- length(FDR_filtered_genes)
print(n_FDR_filtered_genes)


#------------------------------------------------------------
# removing batch effects from count data

# read in $E voom counts generated in DEG pipeline
voom_counts = read.csv(voom_counts_path,header=T,row.names=1)
names(voom_counts) <- sub("^X", "", names(voom_counts))
dim(voom_counts)

# select counts for genes passing the FDR threshold
voom_counts_filtered <- voom_counts %>%
	filter(row.names(voom_counts) %in% FDR_filtered_genes)
dim(voom_counts_filtered) # verify dims
voom_counts_filtered[1:5,1:5]

#Load traits 
pheno1 =read.csv(pheno_path,header=T,row.names=1)
pheno1 <- pheno1[,1:14] # select up to PC2
head(pheno1)

lib_batches <- pheno1$Library_Prep_Batch # select batches as vector
sites <- pheno1$site
sexes <- pheno1$Gender
covariates <- subset(pheno1, select=c("GC","Agilent_RINe"))# additional continuous covariates to adjust
covariates <- data.matrix(covariates)

head(pheno1)
pheno1$asthmaCaseControl <- 0
pheno1$asthmaCaseControl[pheno1$Asthma == "Case"] <- 2
pheno1$asthmaCaseControl[pheno1$Asthma == "Control"] <- 1
pheno1$asthmaCaseControl <- as.factor(pheno1$asthmaCaseControl)

all.equal(colnames(voom_counts_filtered), rownames(pheno1))

# design without library prep batch and covariates GC content and RIN
design <- model.matrix(~asthmaCaseControl+Age+as.factor(site)+PC1+PC2, data=pheno1) #as.factor(Gender)
batch_corrected_counts <- removeBatchEffect(voom_counts_filtered,batch=lib_batches,batch2=sexes,covariates=covariates,design=design)
batch_corrected_counts[1:5,1:5]
write.table(batch_corrected_counts,file="voom_counts_batch_RIN_GC_corrected.csv",sep=",",row.names=T,col.names=NA,quote=F)

# plot counts before and after correction
pdf(file=paste(subset_name,"lib_rin_gc_voom_counts_before&after_rmvbatcheffects_mdsplot.pdf"))
par(mfrow = c(2,3))
p1 <- plotMDS(voom_counts_filtered, labels=pheno1$Library_Prep_Batch, col=pheno1$Library_Prep_Batch, main="Library Prep Batch")
p2 <- plotMDS(voom_counts_filtered, labels=pheno1$Agilent_RINe, col=pheno1$Agilent_RINe,main="RIN")
p3 <- plotMDS(voom_counts_filtered, labels=pheno1$GC, col=pheno1$GC,main="GC Content")


p4 <- plotMDS(batch_corrected_counts, labels=pheno1$Library_Prep_Batch, col=pheno1$Library_Prep_Batch,main="Library Prep Batch")
p5 <- plotMDS(batch_corrected_counts, labels=pheno1$Agilent_RINe, col=pheno1$Agilent_RINe,main="RIN")
p6 <- plotMDS(batch_corrected_counts, labels=pheno1$GC, col=pheno1$GC,main="GC Content")
dev.off()

pdf(file=paste(subset_name,"age_sex_site_voom_counts_before&after_rmvbatcheffects_mdsplot.pdf"))
par(mfrow = c(2,3))
p7 <- plotMDS(voom_counts_filtered, labels=pheno1$Age, col=pheno1$Age,main="Age")
p8 <- plotMDS(voom_counts_filtered, labels=pheno1$Gender, col=unclass(factor(pheno1$Gender)),main="Sex")
p9 <- plotMDS(voom_counts_filtered, labels=pheno1$site, col=unclass(factor(pheno1$site)),main="site")

p10 <- plotMDS(batch_corrected_counts, labels=pheno1$Age, col=pheno1$Age,main="Age")
p11 <- plotMDS(batch_corrected_counts, labels=pheno1$Gender, col=unclass(factor(pheno1$Gender)),main="Sex")
p12 <- plotMDS(batch_corrected_counts, labels=pheno1$site, col=unclass(factor(pheno1$site)),main="site")

dev.off()

pdf(file=paste(subset_name,"asthma_voom_counts_before&after_rmvbatcheffects_mdsplot.pdf"))
p0 <- plotMDS(voom_counts_filtered, labels=pheno1$Asthma, col=unclass(factor(pheno1$Asthma)),main="Asthma Status")
p00 <- plotMDS(batch_corrected_counts, labels=pheno1$Asthma, col=unclass(factor(pheno1$Asthma)),main="Asthma Status")
dev.off()
#------------------------------------------------------------
# Format Data for WGCNA

#Load traits 
pheno =read.csv(pheno_path,header=T,row.names=1)
all_traits = pheno[,1:14] #drop ancestry pcs 2-32
head(all_traits)

# input batch corrected counts (gene rows, sample ID columns)
# format expression data for consensus analysis
expr_data = as.data.frame(t(batch_corrected_counts))

# cluster on Euclidean distance
sample_tree = hclust(dist(expr_data),method = "average")	


pdf(file = paste(c(subset_name, "SampleClustering.pdf"),collapse=""), width = 11, height = 8.5);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
p1<-plot(sample_tree, main = "Sample_clustering", xlab="", sub="", cex = 0.7)+
abline(h=sample_outlier_cutoff,col="red")
dev.off();


# Determine cluster under the line
# cut branches at cutHeight keeping branches with minSize of 10
clust = cutreeStatic(sample_tree, cutHeight = sample_outlier_cutoff, minSize = 10)

# 0 contains cut samples,1 contain keep samples
table(clust)
keepSamples = (clust==1)
discardSamples = (clust==0)
datExpr = expr_data[keepSamples, ]

discard_datExpr = expr_data[discardSamples, ]
discarded_samples<- pheno[rownames(discard_datExpr),]
write.table(discarded_samples,file="discarded_samples.csv",sep=",",row.names=T,col.names=NA,quote=F)

#Re-create traits table with only samples that passed clustering
keep_samples <- intersect(rownames(all_traits), rownames(datExpr))
wgcna_traits <- pheno[keep_samples,1:14]
write.table(wgcna_traits,file="wgcna_pheno.csv",sep=",",row.names=T,col.names=NA,quote=F)


# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
n_genes = ncol(datExpr)
n_samples = nrow(datExpr)


save(datExpr,wgcna_traits,n_genes,n_samples,file = "wgcna_data_input.RData")

#------------------------------------------------------------
# scale free index and mean connectivity plots for determining power

# load expression data and traits
wgcna_data = load("wgcna_data_input.RData")

# choose soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file = paste(c(n_FDR_filtered_genes,subset_name,"scale_free_fit_mean_connectivity_plots.pdf"),collapse=""), width = 11, height = 8.5)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

p1<-plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

p1<-p1+abline(h=0.90,col="red") # using a R^2 cutoff of 0.90

p1<-p1+plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(p1)
dev.off()

#-----------------------------------------------------------
#disableWGCNAThreads()

net = blockwiseModules(datExpr, 
maxBlockSize= 8000, # default is 5000
power = power_beta,
networkType = network_type,
TOMType = tom_type,
deepSplit= deep_split, # default is 2
pamStage= pam_stage,
minModuleSize = min_module_size,
numericLabels = TRUE, 
saveTOMs = TRUE,
saveTOMFileBase = "tom_file_base",
verbose = 3)

save(net, file=paste(subset_name,"_net.RData"))

table(net$colors) # top row number of modules, bottom number of genes, 0 is for genes outside of modules

pdf(file = paste(c(subset_name,"cluster_dendogram_with_module_colors.pdf"),collapse=""), width = 11, height = 8.5)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()


module_labels = net$colors
module_colors = labels2colors(net$colors)
MEs = net$MEs
gene_tree = net$dendrograms[[1]]
save(net,MEs, module_labels, module_colors, gene_tree, file = paste(subset_name,"networkconstruction_auto.RData"))


#------------------------------------------------------------
#Relating Modules to traits
net_data = load(paste(subset_name,"networkconstruction_auto.RData"))

datTraits = wgcna_traits[,c(0,1,2,5,7,10,13,14)] # drop categorical traits for heatmap


#sample clustering dendrogram relating traits
pdf(file = paste(subset_name,"sample_dendogram_with_trait_heatmap.pdf"), width = 11, height = 8.5)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(datTraits),
main = "Sample dendrogram and trait heatmap")
dev.off()



# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, module_colors)$eigengenes
MEs = orderMEs(MEs0) # reorder MEs so that similar ones are next to eachother
moduleTraitCor = cor(MEs, datTraits, use = "p");# weighted pearson correlatiom
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, n_samples);

#sizeGrWindow(10,6)
pdf(file = paste(subset_name,"module_eigengene_trait_heatmap.pdf"), width = 11, height = 8.5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()


# merge network info together (gene module numbers & colors)
network_data <- data.frame(module_labels)
network_data <- cbind(network_data,data.frame(module_colors))
network_data <- rownames_to_column(network_data,"gene_ID")
write.table(network_data,paste(subset_name,"module_assignments.csv"),sep=",",row.names=F)


#----------------------------------------------------------------------------------
# Determining Module Expression
network_data <- column_to_rownames(network_data,"gene_ID")
batch_corrected_counts_modules <- merge(network_data,batch_corrected_counts, by=c(0)) # merge module info with batch corrected counts
dim(batch_corrected_counts_modules)

sorted_counts_modules <- batch_corrected_counts_modules[order(batch_corrected_counts_modules$module_labels),] # sort my module label

module_expr <- sorted_counts_modules %>%
	group_by(module_colors) %>%
	summarise_if(is.numeric,mean,na.rm=TRUE)
	
module_expr <- as.data.frame(module_expr)
module_expr <- module_expr[!(module_expr$module_labels==0),] # drop expression value for unassigned genes

write.table(module_expr,paste(subset_name,"_module_expression.csv"),sep=",",row.names=F)

#---------------------------------------------------------------------------------------------
#Determining overlapping modules
network_data <- read.csv(paste(subset_name,"module_assignments.csv"),header=T)
dim(network_data)


color_numbers <- unique(network_data[,c("module_labels","module_colors")])
color_counts <- data.frame(table(network_data$module_colors))
names(color_counts) <- c("module_colors","gene_count")
network_summary <- merge(color_numbers,color_counts,by="module_colors")
write.table(network_summary , file = paste(subset_name,"module_number&color_gene_counts.csv"), sep = ",",row.names=F,col.names=T)

network_color_key <- table(network_data$module_colors,network_data$module_labels)
write.csv(network_color_key , file = paste(subset_name,"color_key.csv"), sep = ",",row.names=T,col.names=NA)




#---------------------------------------------------------
# Determine hubs for each module
#selecting for gene with highest connectivity looking at all genes in file

network_data = load(paste(subset_name,"networkconstruction_auto.RData"))
colorh <- net$colors

wgcna_data = load("wgcna_data_input.RData")
datExpr = expr_data[[1]]

# find hubs not including module 0 (genes not included in modules)
hubs = chooseTopHubInEachModule(datExpr$data, colorh, omitColors=0, power = power_beta, type=network_type)
write.table(data.frame(hubs), file = paste(subset_name,"hub_genes.csv"), sep = ",",row.names=T,col.names=NA)


#---------------------------------------------------------------------------------------------
# finding DEGs in modules
library(VennDiagram)


DE_results_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_results.csv"
module_results_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/11JUL2022/Run3/module_assignments.csv"

DE_res <-read.csv(DE_results_path,header=T)

# get genes from DE analysis with q < 0.05
DE_res_q005 <- DE_res %>% 
	filter(q.value.full_n536 <0.05)
DEGs_q005<- DE_res_q005$gene_ID

# get genes from DE analysis with p<0.05
DE_res_p005<- DE_res %>% 
	filter(p.value.full_n536<0.05)
DEGs_p005<- DE_res_p005$gene_ID


modules <- read.csv(module_results_path ,header=T)
n_modules <- length(unique(modules$module_labels))-1

common_gene_count <- data.frame(matrix(NA, nrow=n_modules+1, ncol=2))
for (i in 0:26){
	print(paste("module",i))
	module_genes <- modules %>%
		filter(module_labels == i) %>%
		select("gene_ID")
	module_genes <- module_genes[["gene_ID"]] # conver to vector
	n_module_genes <- length(module_genes)
	
	overlap <- calculate.overlap(x=list(DEGs_q005,module_genes)) #find DEGs in current module
	common_genes <-overlap$a3
	n_common_genes<-length(overlap$a3)
	
	overlap2 <- calculate.overlap(x=list(DEGs_p005,module_genes)) #find DEGs in current module
	common_genes2 <-overlap2$a3
	n_common_genes2 <-length(overlap2$a3)
	
	print(common_genes)
	print(n_common_genes)
	common_gene_count[i+1,1] <- i
	common_gene_count[i+1,2] <- n_common_genes
	common_gene_count[i+1,3] <- n_common_genes2

	

}
names(common_gene_count) <- c("module_labels","q005_genes","p005_genes")
common_gene_count

#q005_common_genes_table <- do.call(rbind,q005_common_genes)
write.table(common_gene_count, file = paste(subset_name,"DEG_count_per_module.csv"), sep = ",",row.names=F,col.names=T)



sessionInfo()
