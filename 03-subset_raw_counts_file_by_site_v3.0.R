
rm(list=ls())

# Active Asthma Case Definition

# input site subset names
subset_name0 <- "Baltimore_n76_"
subset_name1 <- "Barbados_n81_"
subset_name2<- "Brazil_n80_"
subset_name3 <- "Chicago_n84_"
subset_name4 <- "Denver_n87_"
subset_name5 <- "Nigeria_n66_"
subset_name6 <-"Washington_DC_n62_"

subset_name <- c(subset_name0,subset_name1,subset_name2,subset_name3,subset_name4,subset_name5,subset_name6)
n_subset <- length(subset_name)

output_path <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/differential_expression_analysis/05MAY2022/active_asthma/by_site/"
setwd(output_path)
sink(paste(c(output_path,"03_subset_raw_counts_by_site_Out.txt"),collapse=""))

# input pheno subsets by site
pheno_path0 <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Baltimore_n76_withPC.csv"
pheno_path1 <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Barbados_n81_withPC.csv"
pheno_path2<- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Brazil_n80_withPC.csv"
pheno_path3 <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Chicago_n84_withPC.csv"
pheno_path4 <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Denver_n87_withPC.csv"
pheno_path5 <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Nigeria_n66_withPC.csv"
pheno_path6 <-"/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/phenotype_data/27APR2022/asthma_active/pheno_Washington_DC_n62_withPC.csv"

pheno_path <- c(pheno_path0, pheno_path1,pheno_path2,pheno_path3,pheno_path4,pheno_path5,pheno_path6)

# read in full raw counts file
raw_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/merged_read_counts/60664_genes/merged_raw_CoCo_counts.txt"
raw = read.table(raw_path,header=T,row.names=1,sep=",")
names(raw) <- sub("^X", "", names(raw))

### Subset raw counts to pheno samples for each site
for (i in 1:n_subset) {
	print(subset_name[i])
	pheno = read.csv(pheno_path[i], header=TRUE, row.names=1)
	final_samples=row.names(pheno)
	raw_final= raw[,c(t(final_samples))]
	write.csv(raw_final,paste(c(output_path,subset_name[i],"_raw_counts.csv"),collapse=""),quote=F,row.names=T)
}

sessionInfo()
savehistory(file=paste(c(output_path,"03_subset_raw_counts_by_site.Rhistory"),collapse=""))
sink()
closeAllConnections()