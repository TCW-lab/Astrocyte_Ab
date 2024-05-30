#WGCNA
#BiocManager::install("WGCNA")
library(rlang, lib = "/projectnb/tcwlab/LabMember/mwu/R/")
library("ggplot2") # For PCA plot
library("WGCNA") # For network analysis
library("DESeq2") # For data transformation
#BiocManager::install("sva")
library("sva") # For batch/covariate correction
#install.packages('RColorBrewer', lib = "/projectnb/tcwlab/LabMember/mwu/R/")
library("GSA") # To extract gene sets
library("erer") # For list writing
library("car") # For Levene's homogeneity test
library("bestNormalize") # for data transformation
library("RColorBrewer") # For PCA color pallette

tom <- readRDS("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/TOM.rds")

df <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/normalized_expression.csv")

# Label modules using colors, grey color is reserved for unassigned genes.
module_color <- labels2colors(tom$colors, zeroIsGrey=TRUE)
names(module_color)<-names(tom$colors)

# Calculate module eigengenes, which represent the expression profile of their respective module. They are the first PCs.
ME <- tom$MEs
# Assign colors to MEs
ME_color <- moduleEigengenes(df[, -1], module_color)$eigengenes
# Order MEs by hierarchical clustering
ME <- orderMEs(ME_color)

# Plot the cluster tree of module eigengenes
plotEigengeneNetworks(ME, "Astrocytes, with KOs \nModule eigengene correlation dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))

# Merge modules based on enrichment validation:
module_merged <- mergeCloseModules(df[, -1], module_color, cutHeight=0.6, verbose=4)

# Combine old module colors with merged module colors for plotting
color_combined <- cbind(module_color, module_merged$colors)

outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/"

# Plot dendrogram with both original module color and merged colors
pdf(paste0(outdir, "Module_Dendrogram_combined_colors.pdf"), width = 10, height = 10)
plotDendroAndColors(tom$dendrograms[[1]], main=paste("Cluster Module Dendrogram with original and merged colors"), colors=color_combined, groupLabels=c("Module colors", "Merged colors"), hang=0.03, dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
dev.off()

# Plot dendrogram with only the merged module color
pdf(paste0(outdir, "Module_Dendrogram_merged_colors.pdf"), width = 10, height = 10)
plotDendroAndColors(tom$dendrograms[[1]], main=paste("Cluster Module Dendrogram with merged colors"), colors=module_merged$colors, groupLabels=c("Merged colors"), hang=0.03, dendroLabels=FALSE, addGuide=TRUE, guideHang=0.05)
dev.off()

# Plot the eigengene correlation from original module
pdf(paste0(outdir, "Module_Eigengene_Correlation_Dendrogram.pdf"), width = 10, height = 10)
plotEigengeneNetworks(ME, "Module Eigengene Correlation Dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))
dev.off()

# Plot the eigengene correlation from merged module
pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Dendrogram.pdf"), width = 10, height = 10)
plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Dendrogram", plotDendrograms=TRUE, plotHeatmaps=FALSE, plotAdjacency=FALSE, excludeGrey=TRUE, marDendro=c(1,4,4,1))
dev.off()

# Plot the eigengene heatmap
pdf(paste0(outdir, "Merged_Module_Eigengene_Correlation_Heatmap.pdf"), width = 10, height = 10)
plotEigengeneNetworks(module_merged$newMEs, "Merged Module Eigengene Correlation Heatmap", plotHeatmaps=TRUE, plotDendrograms=FALSE, plotAdjacency=FALSE, excludeGrey=FALSE, marHeatmap=c(8,10,3,3), xSymbols=names(module_merged$newMEs), ySymbols=names(module_merged$newMEs), printAdjacency = TRUE)
dev.off()

#Filtered input expression data from Feature Count
df_filtered <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/expression_filtered.csv")

library(data.table)

# Save gene id and module color infromation from merged module into a table
modules_df<-data.table(gene_id=names(module_merged$colors),module=module_merged$colors)

# Merge the ME group with the expression data on gene_id
modules_df<-merge(modules_df, df_filtered, by.x = 'gene_id', by.y = 'gene_name')

# Order by ME group
modules_df<-modules_df[order(module)]

#Convert to gene ID to EnsembleID
#Save the gene ids to csv for online conversion using Biotools
#write.csv(modules_df$gene_id, "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/genes_to_convert.csv", row.names = FALSE, quote = FALSE)

gene_symbol <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/gene_names_converted.csv")

#Keep the gene ids if conversion to gene symbols failed
empty_rows <- gene_symbol$gene_symbol == ""
gene_symbol$gene_symbol[empty_rows] <- gene_symbol$gene_id[empty_rows]

modules_df$gene_symbol <- gene_symbol$gene_symbol

#write.csv(modules_df, "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/module_genes.csv", row.names = FALSE, quote = FALSE)