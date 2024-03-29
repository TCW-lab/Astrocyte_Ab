#!/bin.usr./Rscript

deseq2_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/"
gsea_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/GSEA/"
if(!file.exists(gsea_root)) {
      dir.create(gsea_root, mode="0755", recursive=TRUE)
    }

deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCtrl,Uptake_33vs44_Ab,Uptake_33vs44_Ctrl,Uptake_44_AbvsCtrl"

deseq2_conf_list <- unlist(strsplit(deseq2_conf_liststr, ","))

gmt_dir <- "/projectnb/tcwlab/MSigDB/"


#####################GSEA#######################


#BiocManager::install('fgsea', lib = library_path)
library(fgsea)

#Takes file index in the deseq2_conf_liststr
run_gsea <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  #Rank data frame by stat
  df_sorted <- df[order(df$stat),]
  ranks <- setNames(df_sorted$stat, df_sorted$gene_name)
  
  for (i in 1:nrow(pathway_df)) {
    pathways <- gmtPathways(paste0(gmt_dir, pathway_df$gmt[i]))
  
    #GSEA
    fgseaRes <- fgsea(pathways, ranks, scoreType='std',nPermSimple = 10000)
    
    #Not including the leading edges
    gsea_stat <- fgseaRes[, -8]
    
    #Leading edges
    gsea_genes <- data.frame(leadingEdge = sapply(fgseaRes$leadingEdge, paste, collapse = ","))
    rownames(gsea_genes) <- fgseaRes$pathway
    
    #Pathways
    topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    
    gsea_subdir <- paste0(gsea_root, comparison, "/")
    if(!file.exists(gsea_subdir)) {
        dir.create(gsea_subdir, mode="0755", recursive=TRUE)
    }
    
    #Name each GSEA result by the comparison name
    gsea_stats_filename <- paste0(gsea_subdir, comparison, "_", pathway_df$gmt[i], "_stats.csv")
    gsea_genes_filename <- paste0(gsea_subdir, comparison, "_", pathway_df$gmt[i], "_leadingedges.csv")
    gsea_pathway_filename <- paste0(gsea_subdir, comparison, "_", pathway_df$gmt[i], "_pathways.csv")
    
    #Save the data frame
    write.csv(as.data.frame(gsea_stat), gsea_stats_filename, row.names = FALSE, quote = FALSE)
    write.csv(gsea_genes, gsea_genes_filename)
    write.csv(topPathways, gsea_pathway_filename)
  }
}

#Define pathways of interest
pathway_df <- data.frame(
      name=c("GO_all", "C2_CP"),
      gmt=c("c5.go.v2023.1.Hs.symbols.gmt","c2.cp.v2023.1.Hs.symbols.gmt"),
      desc=c("GO All" ,"C2 CP"),
      category=c("C5", "C2"))


#Two-way plot for uptake experiment
for (i in 1:10) {
  run_gsea(i)
}


############################Pathway Visualization#########################
library(ggplot2)

pathway_df_all <- data.frame()

for (i in 1:10) {
  comparison <- deseq2_conf_list[i]
  for (j in 1:nrow(pathway_df)) {
    pathway_res <- paste0(gsea_root, comparison, "/", comparison, "_", pathway_df$gmt[j], "_stats.csv")
    df <- read.csv(pathway_res)
    df_sorted <- df[order(df$NES), ]
    top_neg_pathway <- head(df_sorted, 10)
    top_pos_pathway <- tail(df_sorted, 10)
  
    top_pathways <- rbind(top_neg_pathway, top_pos_pathway)
    top_pathways[, 8] <- pathway_df$name[j]
    top_pathways[, 9] <- deseq2_conf_list[i]
    
    #Get a single data frame that contains all top pathways from all comparisons
    pathway_df_all <- rbind(pathway_df_all, top_pathways)
    
    #Pathway visualization
    outdir <- '/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/Pathway/'
    
    plot_dir <- paste0(outdir, comparison, "/", pathway_df$gmt[j], "/")
    if(!file.exists(plot_dir)) {
        dir.create(plot_dir, mode="0755", recursive=TRUE)
    }
    
    #Bar plot
    bar <- ggplot(top_pathways, aes(x = NES, y = pathway, fill = -log10(padj))) +
      geom_col() +
      scale_fill_gradientn(colors = c("blue", "red")) +
      theme(panel.margin=unit(.05, "lines"),  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
      labs(title = comparison, x = "", y = "Pathway")

    #Save as pdf
    pdf(paste0(plot_dir, comparison, "_", pathway_df$gmt[j], "_bar_plot.pdf"), height = 10, width = 10)
    print(bar)
    dev.off()
    
    #Dotplot
    dot <- ggplot(top_pathways, aes(x=NES, y=pathway, colour=padj, size=size)) +
      geom_point() +
      expand_limits(x=0) +
      labs(title = comparison, x="NES", y="GO term", colour="p_adj", size="Count")
  
    #Save as pdf
    pdf(paste0(plot_dir, comparison, "_", pathway_df$gmt[j], "_dot_plot.pdf"), height = 10, width = 10)
    print(dot)
    dev.off()
  }
}

colnames(pathway_df_all)[8] <- "Pathway"
colnames(pathway_df_all)[9] <- "Comparison"

#Example
tmp <- pathway_df_all[pathway_df_all$Comparison == "Degrade_33_8vs24hr", ]


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
  
  if(is.null(pathway_names))pathway_names=res_fgsea[order(pval)]$pathway
  
  lelist<-LeadingEdges(res_fgsea[pathway%in%pathway_names])
  
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
    p <- ggplot(res_fgsea[pathway%in%pathway_names][,x:=1][,y:=1],aes_string(x='x',y='x'))+
      geom_point(aes_string(size='size',col=col.var))+
      geom_text_repel(aes(label=pathway))+
      scale_color_gradient2(low = cols[1],high = cols[max(1:length(cols))],midpoint = 0,limits=c(-abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var]))),
                                                                                                 abs(as.numeric(as.vector(res_fgsea[pathway%in%pathway_names][,..col.var])))))+
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



run_emma <- function(x){
  comparison <- deseq2_conf_list[x]
  dir <- paste0(deseq2_root, comparison)
  subdir <- list.dirs(dir)[2]
  df <- read.csv(paste0(subdir, "/", "Results_GTFAnnotated_NoGeneIDDuplicates.csv"))
  
  outdir <- '/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/output/Pathway/'
  
  #Rank data frame by stat
  df_sorted <- df[order(df$stat),]
  ranks <- setNames(df_sorted$stat, df_sorted$gene_name)
  
  for (i in 1:nrow(pathway_df)) {
    pathways <- gmtPathways(paste0(gmt_dir, pathway_df$gmt[i]))
  
    #GSEA
    fgseaRes <- fgsea(pathways, ranks, scoreType='std',nPermSimple = 10000)
    
    #Top 40 pathways in either directions
    activated <- head(fgseaRes[order(fgseaRes$NES, decreasing = TRUE), ], 40)
    suppressed <- head(fgseaRes[order(fgseaRes$NES, decreasing = FALSE), ], 40)
    
    top_pathways <- data.frame()
    top_pathways <- rbind(activated, suppressed)
    
    emma <- emmaplot(top_pathways, label.size = 1.5)
    
    emma_dir <- paste0(outdir, comparison, "/", pathway_df$gmt[i], "/")
    if(!file.exists(emma_dir)) {
        dir.create(emma_dir, mode="0755", recursive=TRUE)
    }

    #Save as pdf
    ggsave(paste0(emma_dir, comparison, "_", pathway_df$gmt[i], "_emma_plot.pdf"), height = 10, width = 10)
  }
}

for (i in 1:10) {
  run_emma(i)
}


#Example
df <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/DESEQ2/Degrade_33_8vs24hr/33(8vs24hr)/Results_GTFAnnotated_NoGeneIDDuplicates.csv")



