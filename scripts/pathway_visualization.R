#!/bin/usr/Rscript

deseq2_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
gsea_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/"
if(!file.exists(gsea_root)) {
      dir.create(gsea_root, mode="0755", recursive=TRUE)
    }

deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCtrl,Uptake_33vs44_Ab,Uptake_33vs44_Ctrl,Uptake_44_AbvsCtrl"

deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

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


library(ggplot2)

top_pathway_df_all <- data.frame()

for (i in 1:length(deseq2_conf_list)) {
  comparison <- deseq2_conf_list[i]
  for (j in 1:nrow(pathway_df)) {
    pathway_res <- paste0(gsea_root, comparison, "/", comparison, "_", pathway_df$gmt[j], "_stats.csv")
    df <- read.csv(pathway_res)
    df_sorted <- df[order(df$NES), ]
    top_neg_pathway <- head(df_sorted, 10)
    top_pos_pathway <- tail(df_sorted, 10)
    
    top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
    top_pathways$Pathway <- pathway_df$name[j]
    top_pathways$Comparison <- deseq2_conf_list[i]
    
    #Get a single data frame that contains all top pathways from all comparisons
    top_pathway_df_all <- rbind(top_pathway_df_all, top_pathways)
    
    #Pathway visualization
    outdir <- '/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/'
    
    plot_dir <- paste0(outdir, comparison, "/", pathway_df$gmt[j], "/")
    if(!file.exists(plot_dir)) {
      dir.create(plot_dir, mode="0755", recursive=TRUE)
    }
    
    #Bar plot
    bar <- ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
      geom_col() +
      scale_fill_gradientn(colors = c("yellow", "red")) +
      theme_light() +
      theme(panel.margin=unit(.05, "lines"),  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
      labs(title = comparison, x = "NES", y = "Pathway")
    
    #Save as pdf
    pdf(paste0(plot_dir, comparison, "_", pathway_df$gmt[j], "_bar_plot.pdf"), height = 10, width = 10)
    print(bar)
    dev.off()
  }
}

#Example
df <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/Uptake_44_AbvsCtrl/Uptake_44_AbvsCtrl_c2.cp.v2023.1.Hs.symbols.gmt_stats.csv")

ggplot(df[1:10,], aes(x = NES, y = pathway, fill = -log10(padj))) +
  geom_col() +
  scale_fill_gradientn(colors = c("yellow", "red")) +
  theme_light() + 
  theme(panel.margin=unit(.05, "lines"),  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(title = "test", x = "", y = "Pathway")


# Find c2 pathways that are overlapped between comparisons
library(ggplot2)

#Get all pathways from each gmt for each comparison
pathway_all_df <- data.frame()

for (i in 1:length(degrade_conf_list)) {
  comparison <- deseq2_conf_list[i]
  gsea_dir <- paste0(gsea_root, comparison, "/")

  temp <- readRDS(paste0(gsea_dir, comparison, "_c2.cp.v2023.1.Hs.symbols.gmt_stats.rds"))
  temp_df <- data.frame()
  temp_df <- rbind(temp_df, temp)
  temp_df$Comparison <- comparison
  
  pathway_all_df <- rbind(pathway_all_df, temp_df)
}


library(stringr)
require(igraph)
require(ggraph)

LeadingEdges<-function(res_fgsea){
  if(all(c('term','n.overlap','genes.overlap')%in%colnames(res_fgsea))){
    l_genes<-str_extract_all(res_fgsea$genes.overlap,'[A-Za-z0-9]+')
    l_genes<-lapply(l_genes, function(x)x[x!='c'])
    names(l_genes)<-res_fgsea$pathway
    return(l_genes)
    
  }else{
    l_genes<-str_extract_all(res_fgsea$leadingEdge,'[A-Za-z0-9]+')
    l_genes<-lapply(l_genes, function(x)x[x!='c'])
    names(l_genes)<-res_fgsea$pathway
    return(l_genes)
  }
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

get_similarity_matrix <- function(leading_edge_list) {
  
  n <- length(leading_edge_list)
  ids<-names(leading_edge_list)
  w <- matrix(NA, nrow=n, ncol=n)
  colnames(w) <- rownames(w) <- names(leading_edge_list)
  
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      w[i,j] <- overlap_ratio(leading_edge_list[ids[i]], leading_edge_list[ids[j]])
    }
  }
  return(w)
}


# get graph of sim
get_igraph <- function(res_fgsea, simmat,leading_edge_list,
                       pathway_names, col.var, min_edge) {
  if(any(duplicated(res_fgsea$pathway)))stop('error: duplicated pathways')
  
  wd <- reshape2::melt(simmat[pathway_names,pathway_names])
  wd <- wd[wd[,1] != wd[,2],]
  # remove NA
  wd <- wd[!is.na(wd[,3]),]
  
  g <- graph.data.frame(wd[, -3], directed=FALSE)
  E(g)$width <- sqrt(wd[, 3] * 5) 
  
  
  
  # Use similarity as the weight(length) of an edge
  E(g)$weight <- wd[, 3]
  g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
  
  res_fgseaf<-res_fgsea[V(g)$name,on='pathway']
  #idx <- unlist(sapply(V(g)$name, function(x) which(x == res_fgseaf$pathway)))
  cnt <- sapply(leading_edge_list, length)
  
  V(g)$size <- cnt[V(g)$name]
  
  colVar <- as.numeric(as.vector(res_fgseaf[V(g)$name, on='pathway'][,.SD,.SDcols=col.var][[1]]))
  
  V(g)$colvar <- colVar
  
  return(g)
}


#plot the graphs
add_category_nodes <- function(p,col.var,cols=cols) {
  locol=cols[1]
  midcol=ifelse(length(cols==3),cols[2],NULL)
  hicol=cols[-1]
  
  p<-p + ggnewscale::new_scale_fill() +geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ colvar,
                                                                   size =~ size)) +
    scale_size_continuous(name = "number of genes",
                          range = c(3, 8) )
  
  if(!is.null(midcol))p<-p+scale_fill_gradient2(low = locol,mid=midcol, high = hicol,name=col.var,
                                                guide = guide_colorbar()) 
  else p<-p+scale_fill_continuous(low = locol, high = hicol,name=col.var,
                                  guide = guide_colorbar())
  
  p<-p+theme(legend.title = element_text(size = 10),
             legend.text  = element_text(size = 10)) +
    theme(panel.background = element_blank()) 
  return(p)
}
add_node_label <- function(p,label.size=label.size,max.overlaps=10) {
  
  p <- p + geom_node_text(aes_(label=~name),
                          size = label.size, repel=TRUE,
                          max.overlaps=max.overlaps)
  
  return(p)
}

#Emmaplot (Alexandre's code)
emmaplot<-function(res_fgsea,
                   pathway_names=NULL, 
                   col.var="NES",
                   show_pathway_of=NULL,
                   min_edge=0.2,
                   label.size=2.5,
                   cols=c('blue','white','red'),
                   max.overlaps=10){
  require('ggrepel')
  
  if(all(c('term','n.overlap')%in%colnames(res_fgsea))){
    res_fgsea[,pathway:=term]
    res_fgsea[,NES:=fold.enrichment]
    
  }
  
  if(is.null(pathway_names))pathway_names=res_fgsea[order(res_fgsea$pval),]$pathway
  
  lelist<-LeadingEdges(res_fgsea[res_fgsea$pathway%in%pathway_names, ])
  
  if(!is.null(show_pathway_of)){
    
    lelist<-lelist[sapply(lelist, function(leadingedges)any(show_pathway_of%in%leadingedges))]
    if(length(lelist)>0){
      pathway_names<-names(lelist)
      
    }else{
      stop('This gene are not found in any leading edges of the given pathways')
    }
  }
  
  if(length(lelist)>1){
    
    simat<-get_similarity_matrix(lelist)
    
    g <- get_igraph(res_fgsea = res_fgsea,
                    pathway_names = pathway_names,
                    simmat = simat,
                    leading_edge_list = lelist,
                    min_edge = min_edge,
                    col.var = col.var
    )
    
    
    p <- ggraph(g, layout='nicely')
    
    p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
    ## add dot
    p <- add_category_nodes(p = p,col.var =col.var,cols=cols)
    
    ## add node label
    
    p <- add_node_label(p = p,label.size=label.size,max.overlaps=max.overlaps)
    
  }else{
    # p <- ggplot(res_fgsea[pathway%in%pathway_names][,x:=1][,y:=1],aes_string(x='x',y='x'))+
    #   geom_point(aes_string(size='size',col=col.var))+
    #   geom_text_repel(aes(label=pathway))+
    #   scale_color_gradient2(low = cols[1],high = cols[max(1:length(cols))],midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var]))),
    #                                                                                abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var])))))+
    #   theme_graph()
    
    p <- ggplot(res_fgsea[,x:=1][,y:=1],aes_string(x='x',y='x'))+
      geom_point(aes_string(size='size',col=col.var))+
      geom_text_repel(aes(label=pathway))+
      scale_color_gradient2(low = cols[1],high = cols[max(1:length(cols))],midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[,..col.var]))),
                                                                                                 abs(as.numeric(as.vector(res_fgsea[,..col.var])))))+
      theme_graph()
    
  }
  
  
  if(!is.null(show_pathway_of)){
    if(length(show_pathway_of)>1){
      return(p+ggtitle(paste('Enriched pathways for selected genes')))
      
    }else{
      return(p+ggtitle(paste('Enriched pathways with', show_pathway_of)))
      
    }
    
  }else{
    return(p)
    
  }
}

# Takes a vector of indecies of deseq config files.
# The first indices will be the reference, where all the pathways that are overlapped will be selected based on that.
get_overlap_pathway <- function(x){
  
  df <- pathway_all_df[pathway_all_df$Comparison %in% deseq2_conf_list[x], ]
  
  #Subset overlaped pathways based on p_value threshold 0.1, absolute value of NES threshold 0.5
  df_subset <- df[abs(df$NES) >= 0.5, ]
  
  #Extract the pathways that are overlapped with reference comparison
  df_ref <- df_subset[df_subset$Comparison == deseq2_conf_list[x][1], ]
  
  #Sort the overlapped pathways by abs(NES)
  df_ref_sorted <- df_ref[order(abs(df_ref$NES)), ]
  
  #Save the pathways that are present in all comparisons
  overlapped_pathway <- data.frame()
  
  for (i in unique(df_ref_sorted$pathway)) {
    temp <- df_subset[df_subset$pathway == i, ]
    if (nrow(temp) == length(x)) {
      overlapped_pathway <- rbind(overlapped_pathway, temp)
    }
  }
  
  return(overlapped_pathway)
}

out_dir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/output/Pathway/Overlapped_Pathways/"

#Pathways in 44 that are overlapped in 33 upon treatment
df <- get_overlap_pathway(c(1, 2))

# Change the order of the x-axis groups
df$Comparison <- factor(df$Comparison, levels = c("Uptake_44_AbvsCtrl", "Uptake_33_AbvsCtrl"))

#Plot heatmap for the top 20 overlapped pathways, one asterisk for p values that are less than 0.1 and two asterisks for p values that are less than 0.05
p <- ggplot(df[1:20, ], aes(x = Comparison, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.1 & padj > 0.05, "*", "")), color = "black") +  
  geom_text(aes(label = ifelse(padj < 0.05, "**", "")), color = "black") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradientn(colors = c("blue", "red")) +
  coord_equal()
ggsave(paste0(out_dir, "Uptake_44vs33_overlapped_pathway_bar_plot.pdf"), plot = p, width = 10, height = 10)

#Emma plot for the top 80 overlapped pathways
emma_df <- df[df$Comparison == "Uptake_33_AbvsCtrl", ]
emma_df <- emma_df[emma_df$NES >= 0, ]
emma <- emmaplot(emma_df[1:80, ], label.size = 3)
ggsave(paste0(out_dir, "Uptake_44vs33_overlapped_pathway_emma_plot.pdf"), plot = emma, height = 10, width = 20)


#Pathways in 33_Saturatedvs24hr that are overlapped in 33_Saturatedvs8hr
get_overlap_pathway(c(7, 9)) # No overlaps

#Pathways in 44_Saturatedvs24hr that are overlapped in 44_Saturatedvs8hr
get_overlap_pathway(c(8, 10)) # No overlaps

#Pathways in 44_24hrvs8hr that are overlapped in 33_24hrvs8hr
df <- get_overlap_pathway(c(5, 6))

# Change the order of the x-axis groups
df$Comparison <- factor(df$Comparison, levels = c("Degrade_33_24hrvsCtrl", "Degrade_44_24hrvsCtrl"))

#Plot heatmap for the top 20 overlapped pathways, one asterisk for p values that are less than 0.1 and two asterisks for p values that are less than 0.05
p <- ggplot(df[1:20, ], aes(x = Comparison, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.1 & padj > 0.05, "*", "")), color = "black") +  
  geom_text(aes(label = ifelse(padj < 0.05, "**", "")), color = "black") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradientn(colors = c("blue", "red")) +
  coord_equal()
ggsave(paste0(out_dir, "Degrade_44vs33_24vs8hr_overlapped_pathway_bar_plot.pdf"), plot = p, width = 10, height = 10)

#Emma plot for the top 80 overlapped pathways
emma_df <- df[df$Comparison == "Degrade_33_24hrvsCtrl", ]
emma <- emmaplot(emma_df[1:80, ], label.size = 3)
ggsave(paste0(out_dir, "Degrade_44vs33_24vs8hr_overlapped_pathway_emma_plot.pdf"), plot = emma, height = 10, width = 20)


#Pathways in 44_Saturatedvs8hr that are overlapped in 33_Saturatedvs8hr
get_overlap_pathway(c(9, 10)) # No overlaps


#Pathways in 44_Saturatedvs24hr that are overlapped in 33_Saturatedvs24hr
df <- get_overlap_pathway(c(7, 8))

# Change the order of the x-axis groups
df$Comparison <- factor(df$Comparison, levels = c("Degrade_33_Saturatedvs24hr", "Degrade_44_Saturatedvs24hr"))

#Plot heatmap for the top 20 overlapped pathways, one asterisk for p values that are less than 0.1 and two asterisks for p values that are less than 0.05
p <- ggplot(df[1:20, ], aes(x = Comparison, y = pathway, fill = NES)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(padj < 0.1 & padj > 0.05, "*", "")), color = "black") +  
  geom_text(aes(label = ifelse(padj < 0.05, "**", "")), color = "black") +  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab(NULL) + ylab(NULL) +
  scale_fill_gradientn(colors = c("blue", "red")) +
  coord_equal()
ggsave(paste0(out_dir, "Degrade_44vs33_48vs24hr_overlapped_pathway_bar_plot.pdf"), plot = p, width = 10, height = 10)

#Emma plot for the top 80 overlapped pathways
emma_df <- df[df$Comparison == "Degrade_33_Saturatedvs24hr", ]
emma <- emmaplot(emma_df[1:80, ], label.size = 3)
ggsave(paste0(out_dir, "Degrade_44vs33_48vs24hr_overlapped_pathway_emma_plot.pdf"), plot = emma, height = 10, width = 20)


