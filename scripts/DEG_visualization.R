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


library(ggplot2)

#Count DEGs in each comparison
deg_up_df <- data.frame()
deg_down_df <- data.frame()

for (i in 1:4) {
  comparison <- deseq2_conf_list[i]
  dir <- paste0(deseq2_root, comparison, "/")
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #p_value threshold
  df_subset <- df[df$padj <= 0.1, ]
  
  deg_up_df[i, 1] <- comparison
  deg_up_df[i, 2] <- nrow(df_subset)
  deg_up_df[i, 3] <- nrow(df_subset[df_subset$log2FoldChange > 0, ])
  deg_up_df[i, 4] <- nrow(df_subset[df_subset$log2FoldChange < 0, ])
}

for (i in 5:12) {
  comparison <- deseq2_conf_list[i]
  dir <- paste0(deseq2_root, comparison, "/")
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #p_value threshold
  df_subset <- df[df$padj <= 0.1, ]
  
  deg_down_df[i, 1] <- comparison
  deg_down_df[i, 2] <- nrow(df_subset)
  deg_down_df[i, 3] <- nrow(df_subset[df_subset$log2FoldChange > 0, ])
  deg_down_df[i, 4] <- nrow(df_subset[df_subset$log2FoldChange < 0, ])
}

colnames(deg_up_df)[1] <- "Comparison"
colnames(deg_up_df)[2] <- "DEG_Count"
colnames(deg_up_df)[3] <- "Up"
colnames(deg_up_df)[4] <- "Down"

colnames(deg_down_df)[1] <- "Comparison"
colnames(deg_down_df)[2] <- "DEG_Count"
colnames(deg_down_df)[3] <- "Up"
colnames(deg_down_df)[4] <- "Down"

#Melt data
plot_df <- tidyr::pivot_longer(deg_up_df, cols = c(Up, Down), names_to = "DEG", values_to = "value")
plot_df$DEG <- factor(plot_df$DEG, levels = c("Up", "Down"))

# Change the order of the x-axis groups
plot_df$Comparison <- factor(plot_df$Comparison, levels = deseq2_conf_list)

ggplot(plot_df, aes(x = Comparison, y = DEG_Count, fill = DEG)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "", x = "Comparison", y = "Number of DEGs (FDR<0.1)") +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


#Volcano Plot
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#Takes data frame and plot title
volcano <- function(df, title){
  p <- EnhancedVolcano(df,
                       lab = df$gene_name,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = title,
                       pointSize = c(ifelse(df$log2FoldChange>3, 3, 1)),
                       pCutoff = 0.1 ,
                       FCcutoff = 0.5,
                       cutoffLineType = 'twodash',
                       cutoffLineCol = 'black',
                       cutoffLineWidth = 0.8,
                       hline = c(0.1, 1e-10),
                       hlineCol = 'black',
                       hlineType = 'longdash',
                       hlineWidth = 1,
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE,
                       legendLabels=c('Not sig.','Log (base 2) FC','FDR', 'FDR & Log (base 2) FC'),
                       legendPosition = "right",
                       legendLabSize = 12,
                       legendIconSize = 5.0,
                       caption = bquote(~Log[2]~ "fold change cutoff, 0.5; FDR cutoff, 0.1"))
}

for (i in 1:length(deseq2_conf_list)) {
  comparison <- deseq2_conf_list[i]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  filename <- paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/volcano_plot/", comparison, "_Volcano Plot.pdf")
  pdf(filename, height = 10, width = 10)
  print(volcano(df, comparison))
  dev.off()
}


#PCA
library(tidyverse)
#BiocManager::install("ggfortify")
library(ggfortify)

#Takes file index in the deseq2_conf_liststr
make_pca <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "deseq2_normalized_counts.csv"))
  
  #Assign genes as rownames
  rownames(df) <- df[, 1]
  
  matrix_data <- as.matrix(df[, -1])
  
  # Perform the PCA
  pca_result <- prcomp(t(matrix_data))
  
  pc_eigenvalues <- pca_result$sdev^2
  
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), variance = pc_eigenvalues) %>% 
    #Percent variance
    mutate(pct = variance/sum(variance)*100) %>% 
    #Cumulative variance explained
    mutate(pct_cum = cumsum(pct))
  
  #Percentage explained by PC Visualization
  prct_pc <- pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct)) +
    geom_line(aes(y = pct_cum, group = 1)) + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Percent variance explained", title = paste0(comparison, " Top PCs")) +
    theme_minimal()
  
  output_subdir <- paste0("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/PCA/", comparison, "/")
  if(!file.exists(output_subdir)) {
    dir.create(output_subdir, mode="0755", recursive=TRUE)
  }
  
  #PC file name
  plot1 <- paste0(output_subdir, comparison, "_", "PC_percentage.pdf")
  
  #Save as pdf
  ggsave(plot1, plot = prct_pc, height = 10, width = 10)
  
  #Shorten the sample name to fit the label
  substrings_to_replace <- c("TCW", "\\.Ast")
  for (i in substrings_to_replace) {
    row.names(pca_result$x) <- gsub(i, "", row.names(pca_result$x))
  }
  
  #Get PC scores
  pc_scores <- as_tibble(pca_result$x, rownames = "Sample")
  
  #Sample Cluster Visualization
  pca_plot <- autoplot(pca_result, data = pc_scores, colour = 'Sample', size = 5, main = paste0(comparison, "_PCA")) + theme_minimal()
  
  #PCA Cluster file name
  plot2 <- paste0(output_subdir, comparison, "_", "PCA_cluster.pdf")
  
  #Save as pdf
  pdf(plot2, height = 10, width = 10)
  print(pca_plot)
  dev.off()
  
  #Get PC dimensions
  pc_loadings <- as_tibble(pca_result$rotation, rownames = "gene")
  
  #Get top 10 genes by PCs
  top_genes <- pc_loadings %>% 
    # Select the PCs of interest
    dplyr::select(gene, PC1, PC2) %>%
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    dplyr::group_by(PC) %>% 
    dplyr::arrange(desc(abs(loading))) %>%
    # take the 10 top rows
    dplyr::slice(1:10) %>% 
    # pull the gene column as a vector
    dplyr::pull(gene) %>% 
    # ensure only unique genes are retained
    unique()
  
  top_loadings <- pc_loadings %>% 
    filter(gene %in% top_genes)
  
  gene_pca <- ggplot(data = top_loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.1, "in"))) +
    geom_text(aes(x = PC1, y = PC2, label = gene),
              nudge_y = 0.005, size = 3) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    labs(title = paste0(comparison, " PCA"), x = "PC1", y = "PC2") +
    theme_minimal()
  
  pdf(paste0(output_subdir, comparison, "_PCA_Eigengene.pdf"), height = 10, width = 10)
  print(gene_pca)
  dev.off()
}

for (i in 1:length(deseq2_conf_list)) {
  make_pca(i)
}