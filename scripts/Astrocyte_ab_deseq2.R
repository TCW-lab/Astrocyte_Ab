#add project rpath to include access to the proj-specific packages:
 proj_rpath <- "/projectnb/tcwlab/software/R/library/4.2.1"
 # BiocManager::install("WGCNA", lib=proj_rpath)
 .libPaths(c(proj_rpath, .libPaths()))

library(RColorBrewer)
library(edgeR)
library(WGCNA) #https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/#cranInstall
library(flashClust)
library(DESeq2)
library(parallel)
ncores <- as.numeric(Sys.getenv("NSLOTS"))
print(paste0("number of cores: ", ncores))
cl <- makeCluster(ncores) # 4)
library(doParallel)
library(org.Hs.eg.db)  #https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
library(AnnotationDbi)  #https://www.rdocumentation.org/packages/AnnotationDbi/versions/1.34.4
#https://github.com/Bioconductor/AnnotationDbi
library(ggplot2) 
library(tidyverse)
library(magrittr)
library(stringr)


run_deseq2 <- function () {
   conf <- read.csv(deseq2_conf, header=T)

   #################### Defines parameters in config file ####################

   for (i in 1:nrow(conf)) {
     refGroup <- str_trim(conf$Reference[i])
     targetGroup <- str_trim(conf$Target[i])
     comparison <- str_trim(conf$Comparison[i])
     namespecifier <- str_trim(conf$NameSpecifier[i])
     cmp_cols <- unlist(strsplit(comparison, '+', fixed=T)) 

     #Redefined to meet the design for these sample IDs
     refSubGroup <- substr(refGroup, 5, nchar(refGroup)) 
     targetSubGroup <- substr(targetGroup, 5, nchar(refGroup))

     # use 'namespecifier in deseq conf to create the subdir:
     deseq2_subdir<-paste0(deseq2_dir, "/", namespecifier, "/")
     if(dir.exists(deseq2_subdir)) {
       unlink(deseq2_subdir, recursive = TRUE)
     }
     dir.create(deseq2_subdir, mode="0755", recursive=TRUE)


     #################### Defines parameters in Metadata and Feature Count ####################

     #cova = metadata
     #expr = feature count

     ref_covaA <- ref_cova %>% 
       filter(grepl(refSubGroup, substr(ref_cova$Sample_Name, 5, nchar(ref_cova$Sample_Name)))) %>%
       select(all_of(select_cols)) #select_cols are the columns specified in the config file (META COL:)

     ref_exprA<- ref_exprs[,colnames(ref_exprs) %in% as.character(ref_covaA$Sample_Name)]

     target_covaA <- target_cova %>% 
       filter(grepl(targetSubGroup, substr(target_cova$Sample_Name, 5, nchar(target_cova$Sample_Name)))) %>%
       select(all_of(select_cols))

     target_exprA<- target_exprs[,colnames(target_exprs) %in% as.character(target_covaA$Sample_Name)]

     # join two data sets:
     exprA <- cbind(target_exprA, ref_exprA)
     covaA <- rbind(target_covaA, ref_covaA)


     ###Remove the low expression genes
     isexpr <- rowSums(cpm(exprA)>1) >= 0.1 * ncol(exprA)
     nolowexprA <- exprA[isexpr,]

     # Defines variates
     # Needs manual adjustment of covariate factor.

     for (j in 1:length(cmp_cols))  {
       switch(cmp_cols[j],
              SampleGroup={
                covaA$Sample_Name <- factor(
                  ifelse(substr(covaA$Sample_Name, 5, nchar(covaA$Sample_Name)) == refSubGroup, refGroup, targetGroup), 
                  levels=c(refGroup, targetGroup))
              },
              RIN={
                covaA$RIN <- as.factor(covaA$RIN)
              }, 
              APOE={
                covaA$APOE <- as.factor(covaA$APOE)
              },
              Time={
                covaA$Time <- as.factor(covaA$Time)
              },

              {
                print(paste0("Warning: Unknown variate, ", cmp_cols[j]))
              }
       )
     }

     print (paste("i=", i, "comparison=", comparison, "nameSpecifier=", namespecifier))


     #################### Defines parameters for Deseq2 analysis ####################

     #  comparison <- "SEX+RIN+disease+A_P_O_E"
     designForm <- formula(paste0("~", comparison))


     # Check variate 
     # varPartResidA <- fitExtractVarPartModel(nolowexprA, designForm, covaA)
     # # here we need to check if the variable returns big enough differences, if so, remove that variable and report it.
     # varPartPDF <- paste0(deseq2_subdir, 'VariancePartion_Covariates.pdf')
     # #    varPartPDF <- "/projectnb/tcwlab/LabMember/yshen16/Project/APOE22/00_fastq_100/DESeq2/4V_Tvs44_M/VariancePartion_Covariates.pdf"
     # pdf(file=varPartPDF)
     # main_label <- paste0("Variance Partitioning on raw counts, Design=", comparison)
     # plotVarPart(varPartResidA ,main=main_label)
     # # main_label <- "Variance Partitioning on raw counts"    
     # # plotVarPart(varPartResidA)
     # dev.off()

     # Main Deseq2 analysis
     tryCatch(
       {
         cov <- unlist(strsplit(comparison, "\\+"))

         xddsMatnolowA <- DESeqDataSetFromMatrix(nolowexprA, covaA[cov], design = designForm)  
         xddsMatA <- DESeqDataSetFromMatrix(exprA, covaA, design=designForm)

         yddsnolowA<-DESeq(xddsMatnolowA)
         resnolowA <- results(yddsnolowA)
         resnolowAordered <- resnolowA[order(resnolowA$padj),]

         ############## Annotate DESeq2 result according to gene id and gene name ##############

         annotation <- t2g[match(rownames(resnolowAordered), t2g$gene_id), ]

         # bind reslult and save a copy: 
         resnolowAorderedAnnotated <- cbind(resnolowAordered,annotation)
         nolowAorderedAnno_file <- paste0(deseq2_subdir, "Results_GTFAnnotated.csv")
         write.csv(as.data.frame(resnolowAorderedAnnotated), file=nolowAorderedAnno_file)

         # order result by gene name and FC Log2foldChange, save the result
         nodup <- resnolowAorderedAnnotated
         nodup$absvalFC <- abs(nodup$log2FoldChange)
         nodup <- nodup[order(nodup$gene_name,-nodup$absvalFC),]
         nodup <- nodup[!duplicated(nodup$gene_name),]
         nodup_file <- paste0(deseq2_subdir, "Results_GTFAnnotated_NoGeneIDDuplicates.csv")
         write.csv(as.data.frame(nodup), file=nodup_file) #final result to deliver

         ########### Normalized Counts from DESeq2 and CPM ##############
         xddsMatA_estimateSizeFactors <- estimateSizeFactors(xddsMatA)
         xddsMatA_deseq2_normalized_counts <- counts(xddsMatA_estimateSizeFactors, normalized=TRUE)
         normalized_count_file <- paste0(deseq2_subdir, "deseq2_normalized_counts.csv")
         write.csv(as.data.frame(xddsMatA_deseq2_normalized_counts), file=normalized_count_file) # final result 2 to deliver


         exprA_count_per_million <- cpm(exprA)
         exprA_cpm_file <- paste0(deseq2_subdir, "count_per_million.csv")
         write.csv(as.data.frame(exprA_count_per_million), file=exprA_cpm_file)

         ########### MA (Mean Average) plot #############
         folderChg_pdf <- paste0(deseq2_subdir, "FoldChange.pdf")
         pdf(folderChg_pdf) # final result 4
         label_str <- paste(targetGroup, " vs ", refGroup, ", design=~", comparison)
         plotMA(resnolowA,ylim=c(-5,5),main=label_str)
         dev.off()
       },
       error = function(e) 
       {
         print(e$message)
       }
     )
   } # end each deseq conf
 }  # end of run_deseq2_analysis



 ################ START of Main Program entry ################


proj_dir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results/"
if(! str_sub(proj_dir, -1, -1) == "/") {
 proj_dir <- paste0(proj_dir, "/")
}

deseq2_conf_liststr <- "APOE33_Treatment_vs_Ctrl,APOE44_Treatment_vs_Ctrl,APOE33_Saturated_vs_Ctrl,APOE44_Saturated_vs_Ctrl,APOE44_vs_APOE33_Treatment,APOE44_vs_APOE33_Saturated,APOE44_vs_APOE33_Ctrl,APOE33_d24_vs_d0,APOE44_d24_vs_d0,APOE44_vs_APOE33_Degradation"
deseq2_conf_list <- unlist(str_split(gsub(" ", "", deseq2_conf_liststr), ","))

conf_root <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/"
if(! str_sub(conf_root, -1, -1) == "/") {
 conf_root <- paste0(conf_root, "/")
}
deseq2_conf_dir <- "deseq2_conf"
if(! str_sub(deseq2_conf_dir, -1, -1) == "/") {
 deseq2_conf_dir <- paste0(deseq2_conf_dir, "/")
}

# now let's read in the config files and redo the analysis: 
# deconf_list <- list.files(paste0(conf_root, "/", deseq2_conf_dir))
# here let's create a new metadata to standardize the config info: 
meta_root <- paste0(proj_dir,  "metadata", "/")
deseq2_root <- paste0(proj_dir, "DESEQ2", "/") 

# load gene anno list: 
anno_dir <- "/projectnb/tcwlab/RawData/Genome2.7.9a/hg38/"
t2g_file <- paste0(anno_dir, "gencode.v26.annotation.gtf")
t2g <- as.data.frame(rtracklayer::import(t2g_file))

#Run genotype covariate
for (c in 1:length(deseq2_conf_list)) {
  
  #Comparison name
  cfg_name <- deseq2_conf_list[c]
  
  #Metadata main dir
  meta_dir <- paste0(meta_root, cfg_name, "/")
  # remove existing metadata/cfg_name folder
  if(dir.exists(meta_dir)) {
   unlink(meta_dir, recursive = TRUE, force=1)
  }
  dir.create(meta_dir, mode="0755", recursive=TRUE)
  
  #Comparison subdir
  orig_cfg_dir <- paste0(conf_root, deseq2_conf_dir, cfg_name)
  
  # Reads config file
  cfg_file <- paste0(orig_cfg_dir, "/", cfg_name, ".cfg")
  cfg_data <- readLines(cfg_file) 
  cfg_data <- str_replace_all(cfg_data, " ", "")
  
  self_ref <- as.numeric(str_split(cfg_data[grep("SELF_REF:", cfg_data)], ":")[[1]][2])
  
  ref_meta_file <- paste0(meta_dir, "ref_metadata.csv")
  file.symlink(paste0(orig_cfg_dir, "/", str_split(cfg_data[grep("META:", cfg_data)], ":")[[1]][2]), ref_meta_file)
  
  target_meta_file <- paste0(meta_dir, "target_metadata.csv")
  file.symlink(paste0(orig_cfg_dir, "/", str_split(cfg_data[grep("META:", cfg_data)], ":")[[2]][2]), target_meta_file)
  
  ref_fc_file <- paste0(meta_dir, "ref_fc.txt")
  file.symlink(paste0(proj_dir, "FeatureCount/featureCounts_clean.txt"), ref_fc_file)
  
  target_fc_file <- paste0(meta_dir, 'target_fc.txt')
  file.symlink(paste0(proj_dir, "FeatureCount/featureCounts_clean.txt"), target_fc_file)
  
  deseq2_dir <- paste0(deseq2_root, cfg_name, "/")  
  if(!dir.exists(deseq2_dir)) {
   dir.create(deseq2_dir, mode="0755", recursive = TRUE)
  }
  
  deseq2_orig <- str_split(cfg_data[grep("DESEQ2_CONF:", cfg_data)], ":")[[1]][2]
  deseq2_conf <- paste0(meta_dir, "deseq2_config.csv")
  file.symlink(paste0(orig_cfg_dir, "/", deseq2_orig), deseq2_conf)
  
  meta_col<- str_split(cfg_data[grep("META_COL:", cfg_data)], ":")[[1]][2]
  select_cols <- unlist(str_split(meta_col,","))
  
  ##################### Reads Metadata and Feature Count #####################
  # load feature counts: 
  ref_exprs <- as.matrix(read.table(ref_fc_file, header=TRUE, sep="\t", row.names=1,as.is=TRUE,check.names = FALSE))
  
  target_exprs <- as.matrix(read.table(target_fc_file, header=TRUE, sep="\t", row.names=1,as.is=TRUE,check.names = FALSE))
  
  # load the entire meta data 
  ref_cova <- read.csv(ref_meta_file, header=TRUE, stringsAsFactors=FALSE, colClasses="character")
  
  target_cova <- read.csv(target_meta_file, header=TRUE, as.is=TRUE, stringsAsFactors=FALSE, colClasses="character")
  
  run_deseq2()
}
