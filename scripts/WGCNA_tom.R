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

 deseq2_conf_liststr <- "Uptake_33_AbvsCtrl,Uptake_44_AbvsCtrl,Uptake_44vs33_Ab,Uptake_44vs33_Ctrl,Degrade_33_24hrvsCtrl,Degrade_44_24hrvsCtrl,Degrade_33_Saturatedvs24hr,Degrade_44_Saturatedvs24hr,Degrade_33_SaturatedvsCtrl,Degrade_44_SaturatedvsCtrl,Degrade_44vs33_24hr,Degrade_44vs33_Saturated"
 deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

 gmt_dir <- "/projectnb/tcwlab/MSigDB/"

 # Determine the soft threshold based on maximum r2 and first minimum connectivity
 sft <- 10
 
 wgcna_tom <- function() {
   wgcna_dir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/"

   #Read the normalized expression data
   df <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/normalized_expression.csv")
   
   outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/WGCNA/"
   
   TOM <-  blockwiseModules(df[, -1], maxBlockSize=100000, minBlockSize=0, minModuleSize=30,
                            corType="bicor", maxPOutliers=0.10, pearsonFallback="individual",
                            power=sft, networkType="signed", TOMType="signed", reassignThreshold=1E-8,
                            mergeCutHeight=0.3, deepsplit=2, numericLabels=TRUE, verbose=4)
   
   saveRDS(TOM, paste0(outdir, "TOM.rds"))
 }
 
 wgcna_tom()
