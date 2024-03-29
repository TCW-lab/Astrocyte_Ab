#!/bin/usr/Rscript

 library(rlang, lib = "/projectnb/tcwlab/LabMember/mwu/R")
 library(WGCNA)
 library(flashClust)
 library(curl)

 deseq2_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
 gsea_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/"
 if(!file.exists(gsea_root)) {
   dir.create(gsea_root, mode="0755", recursive=TRUE)
 }

 deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCtrl,Uptake_33vs44_Ab,Uptake_33vs44_Ctrl,Uptake_44_AbvsCtrl"

 deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

 gmt_dir <- "/projectnb/tcwlab/MSigDB/"

 sft_list <- c(14, 14, 14, 14, 16, 18, 18, 14, 14, 12)

 wgcna_dendro <- function(x) {
   comparison <- deseq2_conf_list[x]
   dir <- paste0(deseq2_root, comparison)
   subdir <- list.dirs(dir)[2]
   df <- read.csv(paste0(subdir, "/", "deseq2_normalized_counts.csv"))

   #Transpose the data
   t_df <- as.data.frame(t(df[, -1]))
   colnames(t_df) <- rownames(df)

   outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/WGCNA"

   softPower <- 14
   adjacency <- adjacency(t_df, power = softPower)

   #Create the Topological Dissimilarity Matrix
   TOM.dissimilarity <- 1-TOMsimilarity(adjacency)

   #Cluster the dissimilarity matrix
   geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

   file_dir <- paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/WGCNA/", comparison, "/")

   #Plot dendrogram
   pdf(paste0(file_dir, comparison, "_dendrogram.pdf"), height = 10, width = 10)
   print(plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
              labels = FALSE, hang = 0.04))
   dev.off()

   #To identify modules
   Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
   #Assign colors to each module
   ModuleColors <- labels2colors(Modules)

   #Plot dendrogram
   pdf(paste0(file_dir, comparison, "_module_dendrogram.pdf"), height = 10, width = 10)
   print(plotDendroAndColors(geneTree, ModuleColors,"Module",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05,
                             main = "Gene dendrogram and module colors"))
   dev.off()

   #Identify eigengenes
   MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
   MEs <- MElist$eigengenes 

   write.csv(MEs, paste0(file_dir, comparison, "_eigengenes.csv"))
 }


 for (i in 1) {
   wgcna_dendro(i)
 } 
