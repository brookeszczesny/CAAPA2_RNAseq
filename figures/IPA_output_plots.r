
upstream_reg_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/Run2_29AUG2022_upstream_regulators.xls"
networks_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/Run2_29AUG2022_networks.xls"
pathways_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/Run2_29AUG2022_canonical_pathways.xls"

dexa_targets_path <-  "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/dexamethasone_targets_Run2_29AUG2022.xls"
dexa_targets_path2 <-  "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/dexamethasone_targets2_Run2_29AUG2022.xls"
dexa_mech_targets_path <-  "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/dexamethasone_mechnetwork_targets_Run2_29AUG2022.xls"

flut_targets_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/fluticasone_propionate_targets_Run2_29AUG2022.xls"
flut_mech_targets_path <-  "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/IPA/Run2_29AUG2022/fluticasone_mechnetwork_targets_Run2_29AUG2022.xls"



DEGs_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/summary/combined_results.csv"
output_path <- "/users/bszczesn/CAAPA2/Multiomics_Figures/Figure2/"
setwd(output_path)

library(dplyr)
library(readxl)
library(tibble)
library(VennDiagram)
library(ggplot2)
library(ggrepel)

# read in combined results file
DEG_results <- read.csv(DEGs_path)
upstream_regs <- read_excel(upstream_reg_path, skip=1)
networks <- read_excel(networks_path, skip=1)
pathways <- read_excel(pathways_path, skip=1)


# plotting upstream regulators with volcano type plot
# color genes with p-value of overlap < 0.05 red

#label top 10 upstream regs
pdf(file="upstream_regulators_plot_with_labels.pdf",width = 8, height = 8)
p1<- ggplot(upstream_regs,aes(x=`Activation z-score`,y=-log10(`p-value of overlap`)))+
	geom_point(aes(x=`Activation z-score`,y=-log10(`p-value of overlap`), color = -log10(`p-value of overlap`)>1.30103))+theme_classic()+
	geom_vline(xintercept=2,color="blue") + geom_vline(xintercept=-2,color="blue") + geom_hline(yintercept=1.3,color="blue") +
	scale_colour_manual(name = NULL, values = setNames(c('red','black'),c(T,F)), labels=c("p-value < 0.05","p-value >= 0.05")) +
	xlab("Activation z-score") + ylab("-log10(P-Value of Overlap)")+
	ggtitle("IPA Upstream Regulators") +
	geom_text_repel(data=upstream_regs[1:10,],aes(`Activation z-score`,-log10(`p-value of overlap`),label=`Upstream Regulator`)) +
	theme(legend.position="bottom")
p1
dev.off()

# label dex, flut, IL4, and TGFB1
pdf(file="upstream_regulators_plot_with_labels4_v0.03.pdf",width = 6, height = 8)
p1<- ggplot(upstream_regs,aes(x=`Activation z-score`,y=-log10(`p-value of overlap`)))+
	geom_point(aes(x=`Activation z-score`,y=-log10(`p-value of overlap`), color = -log10(`p-value of overlap`)>1.30103))+theme_classic()+
	geom_vline(xintercept=2,color="blue") + geom_vline(xintercept=-2,color="blue") + geom_hline(yintercept=1.3,color="blue") +
	scale_colour_manual(name = NULL, values = setNames(c('red','black'),c(T,F)), labels=c("p-value < 0.05","p-value >= 0.05")) +
	xlab("Activation z-score") + ylab("-log10(P-Value of Overlap)")+
	#ggtitle("IPA Upstream Regulators") +
	geom_text_repel(data=upstream_regs[c(1,2,4,6),],aes(`Activation z-score`,-log10(`p-value of overlap`),label=`Upstream Regulator`)) +
	theme(legend.position="bottom")
p1
dev.off()


# label dex, IL4
pdf(file="upstream_regulators_plot_with_labels2_v0.03.pdf",width = 6, height = 8)
p1<- ggplot(upstream_regs,aes(x=`Activation z-score`,y=-log10(`p-value of overlap`)))+
	geom_point(aes(x=`Activation z-score`,y=-log10(`p-value of overlap`), color = -log10(`p-value of overlap`)>1.30103))+theme_classic()+
	geom_vline(xintercept=2,color="blue") + geom_vline(xintercept=-2,color="blue") + geom_hline(yintercept=1.3,color="blue") +
	scale_colour_manual(name = NULL, values = setNames(c('red','black'),c(T,F)), labels=c("p-value < 0.05","p-value >= 0.05")) +
	xlab("Activation z-score") + ylab("-log10(P-Value of Overlap)")+
	#ggtitle("IPA Upstream Regulators") +
	geom_text_repel(data=upstream_regs[c(1,2),],aes(`Activation z-score`,-log10(`p-value of overlap`),label=`Upstream Regulator`)) +
	theme(legend.position="bottom")
p1
dev.off()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------



# plotting Networks by score
png(filename="networks_plot.png",width = 768, height = 768)
p2 <- plot(networks$ID, networks$Score, main="IPA networks", xlab="Network ID", ylab="Score", pch=20)
dev.off


# plotting canonical pathways with volcano type plot
png(filename="canonical_pathways_plot.png",width = 768, height = 768)
p3<-plot(pathways$"z-score", pathways$"-log(p-value)",main="IPA Canonical Pathways",
xlab="z-score",ylab="-log10(P-Value)", col = ifelse(pathways$"-log(p-value)">1.3,'red','black'), pch = 20, xlim=c(-3,3) )
legend(x="topright",legend=c("p-value < 0.05","p-value >= 0.05"), fill=c("red", "black"), cex=0.8)
abline(v = 2 , col = "blue")
abline(v = -2 , col = "blue")
abline(h = 1.3, col = "blue")
dev.off()


# plotting canonical pathways with volcano type plot
png(filename="canonical_pathways_with_labels_plot.png",width = 768, height = 768)
p3<-plot(pathways$"z-score", pathways$"-log(p-value)",main="IPA Canonical Pathways",
xlab="z-score",ylab="-log10(P-Value)", col = ifelse(pathways$"-log(p-value)">1.3,'red','black'), pch = 20 ,xlim=c(-3,3))
text(pathways$"z-score", pathways$"-log(p-value)"-0.1,labels=pathways$"Ingenuity Canonical Pathways")
legend(x="topright",legend=c("p-value < 0.05","p-value >= 0.05"), fill=c("red", "black"), cex=0.8)
abline(v = 2 , col = "blue")
abline(v = -2 , col = "blue")
abline(h = 1.3, col = "blue")
dev.off()





# for each row of DEG table 
# get the gene name from the gene_ID column and add to table

DEG_results <- DEG_results %>%
	mutate(gene_name=substring(gene_ID,17))
head(DEG_results)
	
	
#select top regulators
num_regulators <- 5

top_regs_table <- upstream_regs[1:num_regulators,c("Upstream Regulator", "Target Molecules in Dataset")]
top_regs <- top_regs_table$"Upstream Regulator"

# add NA columns for each regulator to DEG table
for (x in 1:num_regulators) { 
	current_reg_name <- top_regs_table[[x,1]]
	DEG_results[current_reg_name] <- NA

}
head(DEG_results)


	
DEG_results <- DEG_results %>%
	mutate(dexamethasone = gene_name %in% top_regs_table[1, "Target Molecules in Dataset"])

DEG_results <- DEG_results %>%
	mutate(dexamethasone = gene_name %in% top_regs_table[1, "Target Molecules in Dataset"])
	
	
#--------------------------------------------------------------------------------------------------------------------------------------
# Look at Dexamethasone and Fluticasone propionate targets	
DEG_results <- read.csv(DEGs_path)

dexa_targets1 <- read_excel(dexa_targets_path, skip=1)
dexa_targets2 <- read_excel(dexa_targets_path2, skip=1)

flut_targets <- read_excel(flut_targets_path, skip=1)

dexa_targets <- rbind(dexa_targets1, dexa_targets2) #combine all targets into single df

	
sorted_DEGs <- DEG_results[order(DEG_results$"q.value.full_n536"),]
sorted_DEGs <- sorted_DEGs %>%
	mutate(ensembl_ID = substr(gene_ID,1,15))# add ensembl ID to DEGs sorted by qvalue in full table

rownames(sorted_DEGs)<-NULL

# find position of a DEG in the results sorted by qvalue for the full set given an ensembl ID
get_DEG_pos <- function(ensembl_id) { 
	pos <- rownames(sorted_DEGs[sorted_DEGs$ensembl_ID==ensembl_id,])
	pos <- as.numeric(pos)
	return(pos)
}


flut_targets <- flut_targets %>% rowwise() %>%
	mutate(DEG_rank = get_DEG_pos(ID)) # add DEG rank to target table


dexa_targets <- dexa_targets %>% rowwise() %>%
	mutate(DEG_rank = get_DEG_pos(ID))
	
dexa_targets %>% print(n=Inf)

# hist of target position in DEG table
f1 <- hist(flut_targets$DEG_rank, main="Distribution of Target DEG rank by q-value in DE Atopic Active Asthmatics n=536",
	xlab="DEG rank by q-value")

d1 <- hist(dexa_targets$DEG_rank, main="Distribution of Target DEG rank by q-value in DE Atopic Active Asthmatics n=536",
	xlab="DEG rank by q-value", ylim = c(0,15))
	
# overlap of targets
overlap_targets <- calculate.overlap(list(dexa_targets$ID,flut_targets$ID))
a1<- length(overlap_targets$a1) # dexa
a2<- length(overlap_targets$a2) # flut
a3<- length(overlap_targets$a3) # Both

v1<- draw.pairwise.venn(area1 = a1,    # Draw pairwise venn diagram
                   area2 = a2,
                   cross.area = a3)
                   
                   
flut_dexa_shared_deg_targets <- sorted_DEGs[sorted_DEGs$ensembl_ID %in% overlap_targets$a3,]
#--------------------------------------------------------------------------------------------------------------------
# Look at Dexamethasone and Fluticasone propionate mechanistic network targets	

DEG_results <- read.csv(DEGs_path)

dexa_mech_targets <- read_excel(dexa_mech_targets_path, skip=1)

flut_mech_targets <- read_excel(flut_mech_targets_path, skip=1)

# overlap of targets
overlap_mech_targets <- calculate.overlap(list(dexa_mech_targets$Target,flut_mech_targets$Target))
m1<- length(overlap_mech_targets$a1) # dexa
m2<- length(overlap_mech_targets$a2) # flut
m3<- length(overlap_mech_targets$a3) # Both

v2<- draw.pairwise.venn(area1 = m1,    # Draw pairwise venn diagram
                   area2 = m2,
                   cross.area = m3)

flut_dexa_shared_mech_deg_targets <-as.data.frame(overlap_mech_targets$a3)

write.csv(flut_dexa_shared_mech_deg_targets, "flut_dexa_shared_mechanisticnetworks_targets.csv")
#--------------------------------------------------------------------------------------------------------------------
# Finding DEGs in upstream regulators

upstream_regs <- read_excel(upstream_reg_path, skip=1)
p005_upstream_regs <- upstream_regs %>% 
	filter(`p-value of overlap`<0.05)

dim(p005_upstream_regs)[1] # number of regulators with p-value of overlap < 0.05

	
DEGs_389 <-sorted_DEGs %>% filter(q.value.full_n536<0.05) # get only DEGs q-value<0.05 in full
	
	
overlap_regs <- calculate.overlap(list(DEGs_389$ensembl_ID, p005_upstream_regs $ID)) # find DEGs in regulators with p-value of overlap < 0.05
a01 <- length(overlap_regs$a1) # DEGs
a02 <- length(overlap_regs$a2) # regulators
a03 <- length(overlap_regs$a3) # Both

deg_reg_shared <- DEGs_389[DEGs_389$ensembl_ID %in% overlap_regs$a3,] # DEGs that are regulators
#--------------------------------------------------------------------------------------------------------------------
# Examine pathways 
pathways <- read_excel(pathways_path, skip=1)

p005_pathways <- pathways %>%
	filter(`-log(p-value)` < 1.3)
dim(p005_pathways)[1] # pathways with p-value < 0.05

predicted_p005_pathways <- p005_pathways %>%
	filter(abs(`z-score`) > 2 )
	
	
#--------------------------------------------------------------------------------------------------------------------	