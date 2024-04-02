# Get WGCNA module expression Pearson Pairwise correlations and pvalues
# specifically we are interested in the modules used for 3 axes of asthma: M4, M5, and M6

# M6 = red
# M5 = green
# M4 = yellow



library(igraph)
library(WGCNA)
library(dplyr)
library(tibble)
library(ggpubr)

output_path <- "/dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/WGCNA/14AUG2022/run4/module_expression_correlation/04AUG2023/"
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


#module expression pairwise pearson correlation with p-values

xp <- data.frame(module_Expr)

M6M5 <- ggscatter(xp, x = 'red', y='green', add = 'reg.line', conf.int = T, xlab = 'M6 Red', ylab = 'M5 green', title ='Module Expression M6 and M5')
M6M5 <- M6M5 +  stat_cor(method = 'pearson')
pdf("M6M5_module_expression_cor.pdf", width=11, height=8.5)
print(M6M5)
dev.off()


M6M4 <- ggscatter(xp, x = 'red', y='yellow', add = 'reg.line', conf.int = T, xlab = 'M6 Red', ylab = 'M4 yellow', title ='Module Expression M6 and M4')
M6M4 <- M6M4 +  stat_cor(method = 'pearson')
pdf("M6M4_module_expression_cor.pdf", width=11, height=8.5)
print(M6M4)
dev.off()

M4M5 <- ggscatter(xp, x = 'yellow', y='green', add = 'reg.line', conf.int = T,xlab = 'M4 yellow', ylab ='M5 green' , title ='Module Expression M4 and M5')
M4M5 <- M4M5+  stat_cor(method = 'pearson')
pdf("M4M5_module_expression_cor.pdf",  width=11, height=8.5)
print(M4M5)
dev.off()