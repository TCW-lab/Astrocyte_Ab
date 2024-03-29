#!/usr/bin/Rscript

library(rlang, lib = "/projectnb/tcwlab/LabMember/mwu/R")
library(WGCNA)
library(flashClust)
library(curl)
library(ggplot2)

deseq2_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
gsea_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/"
if(!file.exists(gsea_root)) {
      dir.create(gsea_root, mode="0755", recursive=TRUE)
    }

deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCtrl,Uptake_33vs44_Ab,Uptake_33vs44_Ctrl,Uptake_44_AbvsCtrl"

deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

gmt_dir <- "/projectnb/tcwlab/MSigDB/"

plot_WGCNA_threshold <- function(x) {
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "deseq2_normalized_counts.csv"))
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/WGCNA"
  
  #Transpose the data
  t_df <- as.data.frame(t(df[, -1]))
  colnames(t_df) <- rownames(df)
  
  #Quality QC, default settings are used here
  gsg <-goodSamplesGenes(t_df)
  
  file_dir <- paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/WGCNA/", comparison, "/")
  if(!file.exists(file_dir)) {
    dir.create(file_dir, mode="0755", recursive=TRUE)
  }
  
  write.csv(summary(gsg), paste0(file_dir, comparison, "_QC_summary.csv"))
  
  #Cluster samples based on distance
  sampleTree <- hclust(dist(t_df), method = "average")
  
  #Save as pdf
  pdf(paste0(file_dir, comparison, "_sample_cluster.pdf"), height = 10, width = 10)
  print(plot(sampleTree, main = paste0(comparison, " Sample Clusters"), sub="", xlab=""))
  dev.off()
  
  #Calculates networks from multiple soft threshold values
  spt <- pickSoftThreshold(t_df)
  
  #Plot the R^2 as a function of soft threshold
  sft_r2 <- ggplot(spt$fitIndices, aes(x = Power, y = SFT.R.sq, label = Power)) + 
    geom_hline(yintercept = 0.8, colour = "red") +
    geom_text(colour = "red") +
    labs(title = paste0(comparison, " Scale independence"), x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit,signed R^2")
  
  #Save as pdf
  pdf(paste0(file_dir, comparison, "_soft_scale.pdf"), height = 10, width = 10)
  print(sft_r2)
  dev.off()
  
  #Plot mean connectivity as a function of soft threshold
  sft_k <- ggplot(spt$fitIndices, aes(x = Power, y = mean.k., label = Power)) + 
    geom_text(colour = "red") +
    labs(title = paste0(comparison, " Mean connectivity"), x = "Soft Threshold (power)", y = "Mean Connectivity")
  
  #Save as pdf
  pdf(paste0(file_dir, comparison, "_mean_connectivity.pdf"), height = 10, width = 10)
  print(sft_k)
  dev.off()
}


for (i in 1:10) {
  plot_WGCNA_threshold(i)
}
