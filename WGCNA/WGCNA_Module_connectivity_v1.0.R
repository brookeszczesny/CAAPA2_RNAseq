
# Get moduleEigengenes and correlations
# Get module expression and correlations

rm(list=ls())

library(edgeR)
library(limma)
library(WGCNA)
library(dplyr)
library(tibble)
library(igraph)

output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/cytoscape/module_connectivity/07DEC2022/"
setwd(output_path)

module_expr_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/Run4\ _module_expression.csv"

#load("/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/wgcna_data_input.RData")
load("/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/Run4\ networkconstruction_auto.RData")


# WGCNA dendro and heatmap of eigengene correlation
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

# Format module expression matrix (module rows, sample ID columns)
module_expr <- read.csv(module_expr_path,header=T)
names(module_expr) <- sub("^X", "", names(module_expr))
module_expr$module_labels <- NULL
module_expr <- column_to_rownames(module_expr,"module_colors") 
module_Expr <- t(module_expr)
module_expr[1:5,1:5]
module_Expr[1:5,1:5]

dim(module_expr)

# heatmap of module expression pearson pairwise correlation
module_cor <- cor(x = module_Expr, y = NULL, method="pearson", "pairwise.complete.obs")
module_cor[1:5,1:5]
heatmap(x = module_cor,symm = TRUE)

# heatmap of MEs pearson pairwise correlation
module_eigen_cor <- cor(x = MEs, y = NULL, method="pearson", "pairwise.complete.obs")
module_eigen_cor[1:5,1:5]
heatmap(x = module_eigen_cor,symm = TRUE)


# make edge table using module expression correlation values
g <- graph.adjacency(module_cor,weighted=T)
get.edgelist(g)
edges <- get.data.frame(g)

# filter out edges connecting to  self
edges <- edges %>%
	filter(from != to)

#remove duplicates
edges <- edges[order(edges$weight),]
dim(edges)
head(edges)

edges <- edges[seq(1, nrow(edges), 2), ]
dim(edges)
head(edges)	
write.csv(edges,"Module_edge_table_expressionCor.csv", row.names=F)
