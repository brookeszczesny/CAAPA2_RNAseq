

# Subset data to relevant data for figures


output_path <- "/users/bszczesn/CAAPA2/Multiomics_Figures/figure_inputs/"
setwd(output_path)

# -------------------------------------------------------------------------------------------------
# Figure 1

library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(patchwork)
library(VennDiagram)
library(gridExtra)
library(ggpubr)
library(tidyr)

normalized_counts_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_voom_norm_counts.csv"
pheno_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"
combined_res <- read.csv( "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/summary/combined_results.csv")

# Input list of genes for gene summary PDF
target_genes<-c(
"ENSG00000115414-FN1"
)

### Violin plots data prep


# Normalized counts as input
normalized_counts=read.csv(normalized_counts_path,header=T,row.names=1)
names(normalized_counts) <- sub("^X", "", names(normalized_counts))

# Phenotype file
colData=read.csv(pheno_path,header=T,row.names=1)


# Only select target genes from normalized counts dataframe
data_keep_rows <- target_genes
normalized_counts_subset <- normalized_counts[rownames(normalized_counts) %in% data_keep_rows, ]

# Transpose normalized counts to combine with phenotype file. Now, first row is sample ID and column is gene name.
normalized_tr=t(normalized_counts_subset)
write.csv(normalized_tr, "Figure1BC_FN1_voom_normalized_counts.csv",row.names=TRUE) # Write out FN1 counts for figure


# Subset phenotype file to include only columns: Site, age_group, and Asthma status
colData_site_asthma=colData %>% select(site,Asthma,age_group)

# Make sure row names of count file and pheno file are the same
all.equal(rownames(normalized_tr), rownames(colData_site_asthma))

# Combine phenotype and normalized counts
normalized_site_asthma=cbind(normalized_tr,colData_site_asthma)
write.csv(normalized_site_asthma, "Figure1BC.csv",row.names=TRUE) # Write out FN1 counts for figure

# Add a new columns for violin colors
normalized_site_asthma$site_asthma <- paste(normalized_site_asthma$site,normalized_site_asthma$Asthma)
normalized_site_asthma$site_asthma_age <- paste(normalized_site_asthma$site,normalized_site_asthma$Asthma,normalized_site_asthma$age_group)
normalized_site_asthma$age_asthma <- paste(normalized_site_asthma$age_group,normalized_site_asthma$Asthma)

# filter filter to adults and peds
normalized_site_asthma_adults <- normalized_site_asthma %>% filter(age_group == "adult")
normalized_site_asthma_peds <- normalized_site_asthma %>% filter(age_group == "ped")

# specify order of appearance in site_age_asthma violins
# site plot
normalized_site_asthma <- normalized_site_asthma %>%
  mutate(site=factor(site,levels=c("Nigeria", "Barbados","Washington DC", "Baltimore", "Chicago", "Denver", "Brazil")) )

# site_asthma plot
normalized_site_asthma <- normalized_site_asthma %>%
  mutate(site_asthma=factor(site_asthma,levels=c( "Nigeria Case","Nigeria Control" , "Barbados Case","Barbados Control","Washington DC Case", "Washington DC Control", 
	"Baltimore Case","Baltimore Control", "Chicago Case","Chicago Control", "Denver Case", "Denver Control", "Brazil Case", "Brazil Control")))

# site_asthma_age plot
normalized_site_asthma <- normalized_site_asthma %>%
  mutate(site_asthma_age=factor(site_asthma_age,levels=c( "Nigeria Case adult", "Nigeria Control adult", "Nigeria Case ped", "Nigeria Control ped",
	"Barbados Case adult","Barbados Control adult",
	"Washington DC Case adult", "Washington DC Control adult",
	"Baltimore Case adult", "Baltimore Control adult",
	"Chicago Case adult", "Chicago Control adult", "Chicago Case ped", "Chicago Control ped",
	"Denver Case ped", "Denver Control ped",
	"Brazil Case adult", "Brazil Control adult")))


   

### forest plots

# get col positions of values from combined res
log2FC_pos <- grep( "log2FC." , names(combined_res))# get column indices for log2FC
SE_pos <- grep( "SE." , names(combined_res))# get column indices for SE
CIL_pos <- grep( "CI.L." , names(combined_res))# get column indices for CIL
CIR_pos <- grep( "CI.R." , names(combined_res))# get column indices for CIR

log2FC_vals <- combined_res[,c(1,log2FC_pos)] %>% pivot_longer(cols=log2FC.full_n536:log2FC.Nigeria_n66,names_to="dataset",values_to="log2FC")%>% mutate(dataset=sub("log2FC.", "",dataset))
SE_vals <- combined_res[,c(1,SE_pos)] %>% pivot_longer(cols=SE.full_n536:SE.Nigeria_n66,names_to="dataset",values_to="SE")%>% mutate(dataset=sub("SE.", "",dataset))
CIL_vals <- combined_res[,c(1,CIL_pos)] %>% pivot_longer(cols=CI.L.full_n536:CI.L.Nigeria_n66,names_to="dataset",values_to="CI.L")%>% mutate(dataset=sub("CI.L.", "",dataset))
CIR_vals <- combined_res[,c(1,CIR_pos)] %>% pivot_longer(cols=CI.R.full_n536:CI.R.Nigeria_n66,names_to="dataset",values_to="CI.R")%>% mutate(dataset=sub("CI.R.", "",dataset))

# combine results info single df
forest_list <- list(log2FC_vals,SE_vals,CIL_vals,CIR_vals)
forestdata <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("gene_ID","dataset")), forest_list) 
#head(forestdata)

# specify order of appearance in forest plot
dataset_order <- c("full_n536","Adult_n376", "Ped_n160","Nigeria_n66", "Barbados_n81","Washington_DC_n62","Baltimore_n76","Chicago_n84","Denver_n87","Brazil_n80")

# function to parse gene_ID and return ensembl ID
get_ensembl_id <- function(gene_id) {
	gene_id <- str_trim(gene_id)
	ensembl_id <-regmatches(gene_id,regexpr("-",gene_id),invert=TRUE)
	ensembl_id<- ensembl_id[[1]][1]
	return(ensembl_id)
}

# function to parse gene_ID and return gene-name
get_gene_name <- function(gene_id) {
	gene_id <- str_trim(gene_id)
	gene_name <-regmatches(gene_id,regexpr("-",gene_id),invert=TRUE)
	gene_name <- gene_name[[1]][2]
	return(gene_name)
}

### Generate gene summary pdf

gene_pdf = ""
#i =1
for (i in 1:length(target_genes)) {
	gene_id <- names(normalized_site_asthma)[i]
	print(i)
	gene_id <- str_trim(gene_id)
	print(gene_id)
	
	gene_name <- get_gene_name(gene_id)	
	ensembl_id <- get_ensembl_id(gene_id)
	print(ensembl_id)
	print(gene_name)
	
	# rename column i with ensembl id
	names(normalized_site_asthma)[i] <- ensembl_id 
	names(normalized_site_asthma_adults)[i] <- ensembl_id 
	names(normalized_site_asthma_peds)[i] <- ensembl_id 

	#### violin plots

	# site and asthma status specific	
	p0 <- ggplot(normalized_site_asthma, aes_string(x= "site_asthma", y= ensembl_id, color="Asthma")) +geom_violin()+
	theme_classic()+
	geom_jitter(shape=16, size=0.5, aes(colour= Asthma)) +
	theme(axis.text.x = element_text(size=rel(0.8),angle = 45, vjust = 1, hjust=1))+theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+
	stat_summary(fun=median, geom="point", size=2, color="black") +
	xlab("")+
	theme(plot.title = element_text(size = rel(0.5)))+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
	
	#ggsave(filename=paste(c("site_asthma_status_violin_",gene_name,".jpeg"),collapse=""),plot=p0,path= output_path)

	#  Case-Control 
	p1 <- ggplot(normalized_site_asthma, aes_string(x= "Asthma", y= ensembl_id, color="Asthma")) +geom_violin()+
	theme_classic()+
	geom_jitter(shape=16,size=0.5,aes(colour= Asthma)) +
	theme(axis.text.x = element_text(size=rel(0.8),angle = 45, vjust = 1, hjust=1))+ theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+	
	stat_summary(fun=median, geom="point", size=2, color="black") +
	ggtitle("Full Group")+
	xlab("")+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
	
	#ggsave(filename=paste(c("asthma_status_violin_",gene_name,".jpeg"),collapse=""),plot=p1,path=output_path)

	#site specific
	p2 <- ggplot(normalized_site_asthma, aes_string(x= "site", y= ensembl_id , color="site")) +geom_violin()+
	theme_classic()+
	geom_jitter(shape=16,size=0.5,aes(colour= site)) +
	theme(axis.text.x = element_text(size=rel(0.8),angle = 45, vjust = 1, hjust=1))+theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+	
	stat_summary(fun=median, geom="point", size=2, color="black") +
	xlab("")+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
	
	#ggsave(filename=paste(c("site_violin_",gene_name,".jpeg"),collapse=""),plot=p2,path=output_path)
	
	#  Adult Case-Control 
	p3 <- ggplot(normalized_site_asthma_adults, aes_string(x= "Asthma", y= ensembl_id , color="Asthma")) +geom_violin()+
	scale_color_manual(values=c("orchid1","chartreuse"))+		
	theme_classic()+
	geom_jitter(shape=16,size=1,aes(colour= Asthma)) +
	theme(axis.text.x = element_text(size=rel(0.8),angle = 45, vjust = 1, hjust=1))+ theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+	
	stat_summary(fun=median, geom="point", size=2, color="black") +
	ggtitle("Adults Only")+
	xlab("")+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
		
	#ggsave(filename=paste(c("asthma_status_adults_violin_",gene_name,".jpeg"),collapse=""),plot=p3,path=output_path)
	
	#  Ped Case-Control 
	p4 <- ggplot(normalized_site_asthma_peds, aes_string(x= "Asthma", y= ensembl_id , color="Asthma")) +geom_violin()+
	scale_color_manual(values=c("purple1","mediumseagreen"))+		
	theme_classic()+	
	geom_jitter(shape=16,size=0.5,aes(colour= Asthma)) +
	theme(axis.text.x = element_text(size=rel(0.8),angle = 45, vjust = 1, hjust=1))+ theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+	
	stat_summary(fun=median, geom="point", size=2, color="black") +
	ggtitle("Pediatrics Only")+
	xlab("")+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
		
	#ggsave(filename=paste(c("asthma_status_peds_violin_",gene_name,".jpeg"),collapse=""),plot=p4,path=output_path)
	
	
	# site , age_group, and asthma status specific	
	p5 <- ggplot(normalized_site_asthma, aes_string(x= "site_asthma_age", y= ensembl_id , color="age_asthma")) +
	geom_violin()+
	scale_color_manual(values=c("orchid","chartreuse","purple1","mediumseagreen"))+		
	theme_classic()+
	geom_jitter(shape=16,size=0.5,aes(colour= age_asthma)) +
	theme(axis.text.x = element_text(size=rel(0.7),angle = 45, vjust = 1, hjust=1))+theme(legend.position = "none")+
	theme(axis.text.y = element_text(size=rel(0.8)))+	
	stat_summary(fun=median, geom="point", size=2, color="black") +
	xlab("")+
	ylab("log2CPM")+
	theme(axis.title.y = element_text(size=rel(0.5)))+
	theme(plot.title = element_text(size = rel(0.5)))
	
	#ggsave(filename=paste(c("site_age_asthma_status_violin_",gene_name,".jpeg"),collapse=""),plot=p5,path= output_path)
	
	
	#### forest plot
	
	# filter data to target gene
	gene_res <- forestdata %>%
		filter(gene_ID == gene_id)%>%
		slice(match(dataset_order,dataset))
		gene_res$dataset <- factor(gene_res$dataset,levels=rev(gene_res$dataset))
	print(gene_res)

	f1 <- ggplot(data=gene_res, aes(y=dataset, x=log2FC,xmin=CI.L,xmax=CI.R))+
		geom_point()+
		scale_y_discrete(breaks=dataset_order <- c("full_n536","Adult_n376", "Ped_n160","Nigeria_n66", "Barbados_n81","Washington_DC_n62","Baltimore_n76","Chicago_n84","Denver_n87","Brazil_n80"),
		labels= c("Full","Adult", "Pediatric", "Nigeria", "Barbados","WashingtonDC", "Baltimore", "Chicago", "Denver", "Brazil"))+
		geom_errorbarh(height=0.1)+
		geom_vline(xintercept=0,color="black",linetype="dashed",alpha=0.5)+
		theme_classic()+
		theme(axis.title.y = element_blank())+
		theme(axis.title.x = element_text(size=rel(.8), margin = margin(t = -100, unit = "pt")))+#axis.title.x = element_text(size=rel(1),vjust=10))+
		ggtitle("Effect Sizes/ 95% CI")
		
	#ggsave(filename=paste(c("forest_plot_",gene_name,".jpeg"),collapse=""),plot=f1,path= output_path)

	# get log2FC (beta), q-value, and p-value for the current gene
	gene_combined_res <- combined_res[combined_res$gene_ID == gene_id,] 
	full_values <- paste(c("	log2FC = ", signif(gene_combined_res$log2FC.full_n536, digits=3), "," ,
		"		p-value = ",signif(gene_combined_res$p.value.full_n536,digits=3), ",",
		"		q-value = ", signif(gene_combined_res$q.value.full_n536,digits=3)),collapse="")
		
	pdf(paste(c(output_path,gene_name,"_violin_plots_summary.pdf"), collapse=""), width=11, height=8.5)
	p12 <-p1+p2
	p120<-p12 -p0
	p34 <- p3+p4
	p345 <-p34- p5
	gene_pdf <- (p120/p345 | f1) + 
	plot_annotation(title= gene_name, subtitle = full_values, theme=theme(plot.title=element_text(size=24))) +
	plot_layout(widths= c(3,1))
	print(gene_pdf)
	dev.off()
	}

write.csv(gene_res, "Figure_1D_FN1_effect_sizes.csv", row.names=FALSE)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

# note that all genes tested in full are tested in all sites; there are 391 genes included in combined_res not tested in full
# add a column to full combined results counting number of sites with p value < 0.05
combined_res <- read.csv( "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/summary/combined_results.csv")

qval_pos <- grep( "q.value." , names(combined_res))# get column indices for qvalues
pval_pos <- grep( "p.value." , names(combined_res))# get column indices for qvalues


combined_res_tested <- combined_res %>%
	filter(!is.na(log2FC.full_n536)) # filter to only genes tested in full and sites


combined_res_pvals <- combined_res_tested[,pval_pos]
combined_res_pvals_sites <- combined_res_pvals[,4:10]

combined_res_tested$n_sites_p005 <- rowSums(combined_res_pvals_sites<0.05, na.rm=FALSE) # do not ignore NAs for sites; does not matter because all NAs are not tested in all 7 sites 


combined_res_tested$volcano_color[combined_res_tested$q.value.full_n536>0.05] ="black" # q >0.05 in full


combined_res_tested$volcano_color[combined_res_tested$q.value.full_n536<0.05 & combined_res_tested$n_sites_p005==0] ="purple" # 0 sites
combined_res_tested$volcano_color[combined_res_tested$q.value.full_n536<0.05 & combined_res_tested$n_sites_p005==1] ="blue" # 1 site
combined_res_tested$volcano_color[combined_res_tested$q.value.full_n536<0.05 & combined_res_tested$n_sites_p005>=2 & combined_res_tested$n_sites_p005<=3] ="green" # 2-3 sites
combined_res_tested$volcano_color[combined_res_tested$q.value.full_n536<0.05 & combined_res_tested$n_sites_p005>=4 & combined_res_tested$n_sites_p005<=7] ="red" # 4-7 sites

write.table(combined_res_tested, file="Figure_1A_Volcano_combined_res_with_nsitesp005.csv", sep=",", row.names=F)

# Volcano Plot genes with q < 0.05 colored by number of sites with p < 0.05
#png(filename="volcano_plot_colored_by_nsite_p<0.05.png",width = 768, height = 768)
#p1<-plot(combined_res_tested$log2FC.full_n536,-log10(combined_res_tested$p.value.full_n536),xlab="log2FC",ylab="-log10(p-Value)", col = combined_res_tested$volcano_color, pch = 20 ) #main= "Full n536 Volcano Plot Colored by Number of Sites with P-value < 0.05 ",
#legend(x="topright",legend=c("0 sites","1 site", "2-3 sites", "4-7 sites", "NS"),fill=c("purple","blue","green","red", "black"))
#dev.off()

#pdf(file="volcano_plot_colored_by_nsite_p<0.05.pdf",width = 8, height = 8)
volcano_pdf<- ggplot(combined_res_tested,aes(x=`log2FC.full_n536`,y=-log10(`p.value.full_n536`)))+
	geom_point(aes(x=`log2FC.full_n536`,y=-log10(`p.value.full_n536`),color=`volcano_color`))+
	theme_classic()+
	scale_colour_manual(name = NULL, values = setNames(c("purple","blue","green","red", "black"),c("purple","blue","green","red", "black")), labels=c("0 sites","1 site", "2-3 sites", "4-7 sites", "NS")) +
	scale_size_manual(values = c(0.2))+
	xlab("log2FC") + ylab("-log10(p-Value)")+
	theme(legend.position=c(.85, .95), legend.key.size = unit(2, "pt"),
	axis.title.x = element_text(margin = margin(t = -100, unit = "pt")))

print(volcano_pdf)
#dev.off()


# filter combined_res table to genes with q < 0.05 and count how many genes occur for each site # with p<0.05
combined_res_q005 <- combined_res_tested %>%
	filter(q.value.full_n536 <0.05)
p005_nsite_freq <- table(combined_res_q005$n_sites_p005)
print("# of sites with p<0.05 frequency table")
print(p005_nsite_freq)



layout <- "
##BBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
AABBBB
"

pdf("Figure1_v3.0.pdf", width=12, height=6)
figure1_pdf <- volcano_pdf + gene_pdf +
	plot_layout(design = layout)
print(figure1_pdf)
dev.off()  


#------------------------------------------------------------------------------------------------------------

rm(list=ls())


library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(patchwork)
library(VennDiagram)
library(gridExtra)
library(ggpubr)
library(tidyr)


output_path <- "/users/bszczesn/CAAPA2/Multiomics_Figures/supplementary_forest_plots/24JAN2023/"
setwd(output_path)

normalized_counts_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/full/full_n536__21831_voom_norm_counts.csv"

pheno_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_full_n536_withPC.csv"


combined_res <- read.csv( "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/summary/combined_results.csv")

# Input list of genes for gene summary PDF
# top 15 genes
target_genes<-c(
"ENSG00000115414-FN1", #1
"ENSG00000127948-POR", #2
"ENSG00000155659-VSIG4", #3
"ENSG00000143502-SUSD4", #4
"ENSG00000244694-PTCHD4", #5
"ENSG00000182601-HS3ST4", #6
"ENSG00000158528-PPP1R9A", #7
"ENSG00000115306-SPTBN1", #8
"ENSG00000172554-SNTG2", #9
"ENSG00000263961-RHEX", #10
"ENSG00000096060-FKBP5", #11
"ENSG00000197694-SPTAN1", #12
"ENSG00000239445-ST3GAL6-AS1", #13
"ENSG00000087460-GNAS", #14
"ENSG00000140937-CDH11" #15
) 



### Violin plots data prep


# Normalized counts as input
normalized_counts=read.csv(normalized_counts_path,header=T,row.names=1)
names(normalized_counts) <- sub("^X", "", names(normalized_counts))

# Phenotype file
colData=read.csv(pheno_path,header=T,row.names=1)


# Only select target genes from normalized counts dataframe
data_keep_rows <- target_genes
normalized_counts_subset <- normalized_counts[rownames(normalized_counts) %in% data_keep_rows, ]

# Transpose normalized counts to combine with phenotype file. Now, first row is sample ID and column is gene name.
normalized_tr=t(normalized_counts_subset)

# Subset phenotype file to include only columns: Site, age_group, and Asthma status
colData_site_asthma=colData %>% select(site,Asthma,age_group)

# Make sure row names of count file and pheno file are the same
all.equal(rownames(normalized_tr), rownames(colData_site_asthma))

# Combine phenotype and normalized counts
normalized_site_asthma=cbind(normalized_tr,colData_site_asthma)

# Add a new columns for violin colors
normalized_site_asthma$site_asthma <- paste(normalized_site_asthma$site,normalized_site_asthma$Asthma)
normalized_site_asthma$site_asthma_age <- paste(normalized_site_asthma$site,normalized_site_asthma$Asthma,normalized_site_asthma$age_group)
normalized_site_asthma$age_asthma <- paste(normalized_site_asthma$age_group,normalized_site_asthma$Asthma)

# filter filter to adults and peds
normalized_site_asthma_adults <- normalized_site_asthma %>% filter(age_group == "adult")
normalized_site_asthma_peds <- normalized_site_asthma %>% filter(age_group == "ped")


### forest plots

# get col positions of values from combined res
log2FC_pos <- grep( "log2FC." , names(combined_res))# get column indices for log2FC
SE_pos <- grep( "SE." , names(combined_res))# get column indices for SE
CIL_pos <- grep( "CI.L." , names(combined_res))# get column indices for CIL
CIR_pos <- grep( "CI.R." , names(combined_res))# get column indices for CIR

log2FC_vals <- combined_res[,c(1,log2FC_pos)] %>% pivot_longer(cols=log2FC.full_n536:log2FC.Nigeria_n66,names_to="dataset",values_to="log2FC")%>% mutate(dataset=sub("log2FC.", "",dataset))
SE_vals <- combined_res[,c(1,SE_pos)] %>% pivot_longer(cols=SE.full_n536:SE.Nigeria_n66,names_to="dataset",values_to="SE")%>% mutate(dataset=sub("SE.", "",dataset))
CIL_vals <- combined_res[,c(1,CIL_pos)] %>% pivot_longer(cols=CI.L.full_n536:CI.L.Nigeria_n66,names_to="dataset",values_to="CI.L")%>% mutate(dataset=sub("CI.L.", "",dataset))
CIR_vals <- combined_res[,c(1,CIR_pos)] %>% pivot_longer(cols=CI.R.full_n536:CI.R.Nigeria_n66,names_to="dataset",values_to="CI.R")%>% mutate(dataset=sub("CI.R.", "",dataset))

# combine results info single df
forest_list <- list(log2FC_vals,SE_vals,CIL_vals,CIR_vals)
forestdata <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c("gene_ID","dataset")), forest_list) 
#head(forestdata)
write.csv(forestdata, "forest_plot_effect_sizes.csv", row.names=FALSE)

# specify order of appearance in forest plot
dataset_order <- c("full_n536","Adult_n376", "Ped_n160","Nigeria_n66", "Barbados_n81","Washington_DC_n62","Baltimore_n76","Chicago_n84","Denver_n87","Brazil_n80")
#c("Full","Adult", "Pediatric", "Nigeria", "Barbados","WashingtonDC", "Baltimore", "Chicago", "Denver", "Brazil")


# function to parse gene_ID and return ensembl ID
get_ensembl_id <- function(gene_id) {
	gene_id <- str_trim(gene_id)
	ensembl_id <-regmatches(gene_id,regexpr("-",gene_id),invert=TRUE)
	ensembl_id<- ensembl_id[[1]][1]
	return(ensembl_id)
}

# function to parse gene_ID and return gene-name
get_gene_name <- function(gene_id) {
	gene_id <- str_trim(gene_id)
	gene_name <-regmatches(gene_id,regexpr("-",gene_id),invert=TRUE)
	gene_name <- gene_name[[1]][2]
	return(gene_name)
}



### Generate gene summary pdf
plot_list <- list()
for (i in 1:length(target_genes)) {
	gene_id <- names(normalized_site_asthma)[i]
	print(i)
	gene_id <- str_trim(gene_id)
	print(gene_id)
	
	gene_name <- get_gene_name(gene_id)	
	ensembl_id <- get_ensembl_id(gene_id)
	print(ensembl_id)
	print(gene_name)
	
	# rename column i with ensembl id
	names(normalized_site_asthma)[i] <- ensembl_id 
	names(normalized_site_asthma_adults)[i] <- ensembl_id 
	names(normalized_site_asthma_peds)[i] <- ensembl_id 

	
	#### forest plot
	
	# filter data to target gene
	gene_res <- forestdata %>%
		filter(gene_ID == gene_id)%>%
		slice(match(dataset_order,dataset))
		gene_res$dataset <- factor(gene_res$dataset,levels=rev(gene_res$dataset))
	print(gene_res)

	f1 <- ggplot(data=gene_res, aes(y=dataset, x=log2FC,xmin=CI.L,xmax=CI.R))+
		geom_point()+
		scale_y_discrete(breaks=c("full_n536","Adult_n376", "Ped_n160","Nigeria_n66", "Barbados_n81","Washington_DC_n62","Baltimore_n76","Chicago_n84","Denver_n87","Brazil_n80"),
		labels= c("Full","Adult", "Pediatric", "Nigeria", "Barbados","WashingtonDC", "Baltimore", "Chicago", "Denver", "Brazil"))+
		geom_errorbarh(height=0.1)+
		geom_vline(xintercept=0,color="black",linetype="dashed",alpha=0.5)+
		theme_classic()+
		theme(axis.title.y = element_blank())+
		#theme(axis.title.x = element_text(size=rel(1), margin = margin(t = -100, unit = "pt")))+#axis.title.x = element_text(size=rel(1),vjust=10))+
		#ggtitle(paste(gene_name,"Effect Sizes/ 95% CI"))
		ggtitle(gene_name)
		
		if (i %in% c(4,6)) {
			plot_title_color <- "#6A51A3" #blueviolet
			} else if (i %in% c(8,11)) {
			plot_title_color <- "#DD3497" #deeppink
			} else if (i %in% c(14,3,2)) {
			plot_title_color <-"#4EB3D3" #light blue
			} else {
			plot_title_color <- "black"
		}
		
	f1<- f1+theme(plot.title = element_text(colour = plot_title_color))
		
	plot_list[[i]] <- f1

	}
	
# Color code titles by function group
# Airway Remodeling purple: FN1, CDH11
# Th2 Inflammation: VSIG4, HS3ST4
# Drug response blue: PTCHD4, SPTBN1, FKBP5
# NA Black: POR, SUSD4, PPP1R9A, SNTG2, RHEX, SPTAN1, ST#GAL6-AS1, GNAS
	
	# FN1 4 ; POR 5; VSIG4 8; SUSD4 7; PTCHD4 14
	# HS3ST4 11; PPP1R9A 9;	SPTBN1 3; SNTG2 10;	RHEX 15 
	# FKBP5 2; SPTAN1 12; ST3GAL6-AS1 13; GNAS 1; CDH11 6
	

	forest2_pdf <- plot_list[[4]] + (plot_list[[5]]+theme(axis.text.y=element_blank())) + (plot_list[[8]]+theme(axis.text.y=element_blank())) + (plot_list[[7]]+theme(axis.text.y=element_blank())) + (plot_list[[14]]+theme(axis.text.y=element_blank())) +
		plot_list[[11]] + (plot_list[[9]]+theme(axis.text.y=element_blank())) + (plot_list[[3]]+theme(axis.text.y=element_blank())) + (plot_list[[10]]+theme(axis.text.y=element_blank())) + (plot_list[[15]]+theme(axis.text.y=element_blank())) +
		plot_list[[2]] + (plot_list[[12]]+theme(axis.text.y=element_blank())) + (plot_list[[13]]+theme(axis.text.y=element_blank())) + (plot_list[[1]]+theme(axis.text.y=element_blank())) + (plot_list[[6]]+theme(axis.text.y=element_blank())) +
		plot_layout(ncol=5)
		
	pdf(paste(c(output_path,"Top15genes_forest_plots_v5.0.pdf"), collapse=""), width=11, height=8.5)
	print(forest2_pdf)
	dev.off()
