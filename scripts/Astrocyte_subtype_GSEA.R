deseq2_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
gsea_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/GSEA/"

if(!file.exists(gsea_root)) {
  dir.create(gsea_root, mode="0755", recursive=TRUE)
}


deseq2_conf_liststr <- "Uptake_33_AbvsCtrl,Uptake_44_AbvsCtrl,Uptake_44vs33_Ab,Uptake_44vs33_Ctrl,Degrade_33_24hrvsCtrl,Degrade_44_24hrvsCtrl,Degrade_33_Saturatedvs24hr,Degrade_44_Saturatedvs24hr,Degrade_33_SaturatedvsCtrl,Degrade_44_SaturatedvsCtrl,Degrade_44vs33_24hr,Degrade_44vs33_Saturated"

deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

uptake_conf_list <- deseq2_conf_list[1:4]

degrade_conf_list <- deseq2_conf_list[5:12]

gmt_dir <- "/projectnb/tcwlab/MSigDB/"


astro_subtype <- read.csv("/projectnb/tcwlab/RefData/Astrocytes_signatures/astrocyte_signatures_trans_human.csv")

astro_subtype_pathway_list <- list()

#Subgroup each GO
for (i in 1:length(unique(astro_subtype$type))) {
  astro_subtype_pathway_list[[i]] <- astro_subtype[astro_subtype$type == unique(astro_subtype$type)[i], ]$Symbol.human
}

names(astro_subtype_pathway_list) <- unique(astro_subtype$type)

#Takes file index in the deseq2_conf_liststr
run_gsea_astrocyte <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  outdir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/"
  
  #Rank data frame by stat
  df_sorted <- df[order(df$stat),]
  ranks <- setNames(df_sorted$stat, df_sorted$gene_name)
  
  pathways <- astro_subtype_pathway_list
  
  #GSEA
  fgseaRes <- fgsea(pathways = astro_subtype_pathway_list, stats = ranks, scoreType='std',nPermSimple = 10000)
  
  #Not including the leading edges
  gsea_stat <- fgseaRes[, -8]
  
  #Leading edges
  gsea_genes <- data.frame(leadingEdge = sapply(fgseaRes$leadingEdge, paste, collapse = ","))
  rownames(gsea_genes) <- fgseaRes$pathway
  
  file_dir <- paste0(outdir, comparison, "/Astrocyte_subtype_signatures/")
  if(!file.exists(file_dir)) {
    dir.create(file_dir, mode="0755", recursive=TRUE)
  }
  
  #Save the data frame
  write.csv(as.data.frame(gsea_stat), paste0(file_dir, "Astrocyte_subtype_signature_pathway_stats.csv"), row.names = FALSE, quote = FALSE)
  write.csv(gsea_genes, paste0(file_dir, "Astrocyte_subtype_signature_pathway_leadingedges.csv"))
}

for (i in 1:length(deseq2_conf_list)) {
  run_gsea_astrocyte(i)
}


#Get a single data frame with the astrocyte subtypes
uptake_astro_subtype_pathway <- data.frame()
degrade_astro_subtype_pathway <- data.frame()

for (i in 1:length(uptake_conf_list)) {
  comparison <- uptake_conf_list[i]
  
  tmp <- data.frame()
  
  tmp <- read.csv(paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/", comparison, "/Astrocyte_subtype_signatures/Astrocyte_subtype_signature_pathway_stats.csv"))
  
  tmp$Comparison <- comparison
  
  uptake_astro_subtype_pathway <- rbind(uptake_astro_subtype_pathway, tmp)
}

ggplot(uptake_astro_subtype_pathway, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_col() +
  scale_fill_gradientn(colors = c("yellow", "red")) +
  facet_grid(. ~ Comparison) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Astrocyte Subtype Signature Pathways Enrichment across Comparisons", x = "Comparison", y = "NES")


for (i in 1:length(degrade_conf_list)) {
  comparison <- degrade_conf_list[i]
  
  tmp <- data.frame()
  
  tmp <- read.csv(paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/", comparison, "/Astrocyte_subtype_signatures/Astrocyte_subtype_signature_pathway_stats.csv"))
  
  tmp$Comparison <- comparison
  
  degrade_astro_subtype_pathway <- rbind(degrade_astro_subtype_pathway, tmp)
}

ggplot(degrade_astro_subtype_pathway, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_col() +
  scale_fill_gradientn(colors = c("yellow", "red")) +
  facet_grid(. ~ Comparison) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Astrocyte Subtype Signature Pathways Enrichment across Comparisons", x = "Comparison", y = "NES")