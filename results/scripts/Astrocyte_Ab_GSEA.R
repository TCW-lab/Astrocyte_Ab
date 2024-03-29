#fgsea analysis (on pseudobulk results)

proj_rpath <- "/projectnb/tcwlab/software/R/library/4.2.1"
#BiocManager::install("fgsea", lib=proj_rpath)
.libPaths(c(proj_rpath, .libPaths()))

library(tidyverse)
library(magrittr)
library(stringr)
library(foreach)
library(fgsea)
library(grid)
library(gridExtra)
library(ggplot2)
library(BiocParallel)

# task_id<-2
# total_task <- 48
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
total_task <- as.integer(Sys.getenv("SGE_TASK_LAST"))

# define the main function to do the job:

run_gsea_analysis <- function () {
  
  # first read from config file:
  conf <- read.csv(deseq2_conf, header=T)
  for (i in 1:nrow(conf)) {
    print(paste0("config = ", i))
    refGroup <- str_trim(conf$Reference[i])
    targetGroup <- str_trim(conf$Target[i])
    namespecifier <- str_trim(conf$NameSpecifier[i])
    refSubGroup <- unlist(strsplit(refGroup, '+', fixed=T))
    targetSubGroup <- unlist(strsplit(targetGroup, '+', fixed=T))
    ref_str <- paste(refSubGroup, collapse ="_")
    target_str <- paste(targetSubGroup, collapse ="_")  
    if(is.null(namespecifier)) {
      namespecifier <- paste0(target_str, "vs", ref_str)
    }
    # use 'namespecifier in deseq conf to create the subdir, the dir has been created in deseq stage:
    deseq2_subdir<-paste0(deseq2_dir, "/", namespecifier, "/")
    if(!file.exists(deseq2_subdir)) {
      print("DESeq2 Analysis result is missing. Please run DESeq2 first!")
    }
    # get gtf annotation file from deseq2 analysis:
    gtf_anno_infile <- paste0(deseq2_subdir, "Results_GTFAnnotated_NoGeneIDDuplicates.csv")
    if(!file.exists(gtf_anno_infile) )  { 
      print(paste0("It seems no DESEQ2 result is avaiable for Config: ", namespecifier))
      next # go to next conf
    }
    anno_result <- read.csv(gtf_anno_infile, header=TRUE)
    # order increasingly by stat:
    anno_result <- anno_result[order(anno_result$stat),]
    ranks <- setNames(anno_result$stat, anno_result$gene_name)
    
    ct=str_sub(ref_str, -1, -1)
    switch(ct # get the last character
           ,M={cellType<-"Microglia"}
           ,A={cellType<-"Astrocyte"}
           ,h={cellType<-"h"}
           ,V={cellTyep<-"V"}
           ,{
             print(paste0("Warning: unsupported cell type: ", ct))
           }
    )
    
    
    #C5: GO
    #GO All: "c5.all.v7.5.1.symbols.gmt"
    
    #CP Canonical Pathways
    #C2 CP: "c2.cp.v2023.1.Hs.symbols.gmt"
    
    pw_conf <- data.frame(
      name=c("GO_all", "C2_CP"),
      gmt=c("c5.go.v2023.1.Hs.symbols.gmt"
            ,"c2.cp.v2023.1.Hs.symbols.gmt"
      ),
      desc=c("GO All"
             ,"C2 CP"
      ),
      category=c("C5", "C2")
    )# end of data.frame
    
    
    for (j in (1:nrow(pw_conf))) {
      if(
        (c*i*j-1) %% total_task==(task_id-1)
      ) {
        x <- (c*i*j-1)
        y <- x %% total_task
        z <- y!=(task_id-1)
        print(paste0("task=", task_id, " c=", c, " i=", i, " j=", j, " x=", x, " y=", y, " z=", z))
        print(paste0("pathway=", j))
        # pathways: 
        pathways <- gmtPathways(paste0(gmt_dir, pw_conf$gmt[j]))
        
        gsea_subdir <- paste0(gsea_dir, namespecifier, "/")
        # print(paste0("gsea_subdir=",gsea_subdir))
        if(!file.exists(gsea_subdir)) {
          dir.create(gsea_subdir, mode="0755", recursive=TRUE)
        }
        fgsea_time <- system.time(
          fgseaRes <- fgsea(pathways, ranks) 
        )
        
        print(paste("Scenario#3 FGSEA Time Elapsed: "))
        print(fgsea_time)
        gsea_csv <- paste0(gsea_subdir, pw_conf$name[j], "_", cellType, "_", namespecifier, ".csv")
        gsea.df <- data.frame(lapply(fgseaRes, as.character), stringsAsFactors=FALSE) 
        write.csv(
          as.data.frame(gsea.df),
          file = gsea_csv)
        topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown)) 
        gsea_pdf <- paste0(gsea_subdir,
                           pw_conf$name[j], "_", cellType, "_TopPathways_",namespecifier, ".pdf")
        pdf(gsea_pdf, 
            height = 20, width = 32)
        #note: used "plotGseaTable" function instead of "plotGseaTable2", as "plotGseaTable2" is no longer supported
        plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5, colwidths = c(20, 6, 2, 2, 2))
        dev.off() ####
        
        print ("Done with calc. all three methods of GSEA parameter combinations. ")
      } else {
        print (paste0("skip c=", c, " i=", i, " j=", j))
      } # end task distrib check
    } # end loop pw_conf        
  } # end of loop conf
} # end run_gsea


# start of the main program:
ncores <- as.numeric(Sys.getenv("NSLOTS"))
if(is.na(ncores)) ncores=1
register(BatchtoolsParam(ncores))

# from fgsea vignettes
# Also, fgsea is parallelized using BiocParallel package. By default the first registered backend returned by bpparam() is used. 
# To tweak the parallelization one can either specify BPPARAM parameter used for bplapply of set nproc parameter, which is a 
# shorthand for setting BPPARAM=MulticoreParam(workers = nproc).

# set pathway file directory: 
gmt_dir <- "/projectnb/tcwlab/MSigDB/"

# set the feature count result dir: 
proj_dir <- "/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/results" #"/projectnb/tcwlab/LabMember/yshen16/Project/APOE_Mye/Micro_trim_100"
if(! str_sub(proj_dir, -1, -1) == "/") {
  proj_dir <- paste0(proj_dir, "/")
}


deseq2_conf_num <- 10
deseq2_conf_liststr <- "Degrade_33_8vs24hr,Degrade_33_8vs48hr,Degrade_33vs44_24hr,Degrade_33vs44_48hr,Degrade_44_8vs24hr,Degrade_44_8vs48hr,Uptake_33_AbvsCrtl,Uptake_33vs44_Ab,Uptake_33vs44_Crtl,Uptake_44_AbvsCrtl"
deseq2_conf_list <- unlist(str_split(gsub(" ", "", deseq2_conf_liststr), ","))

subdir <- c("genotype_comparison", "time_comparison", "treatment_comparison")

for (i in subdir) {
  meta_root <- paste0(proj_dir,  "metadata", "/", i, "/")
  
  # loop through the config files:
  for (c in deseq2_conf_liststr) {
    if (file.exists(c)) {
      cfg_name <- deseq2_conf_list[c]
      print(paste0("deseq2_config=",c, " name=", cfg_name))
      meta_dir <- paste0(meta_root, cfg_name, "/")
      ref_meta_file <- paste0(meta_dir, "ref_metadata.csv")
      target_meta_file <- paste0(meta_dir, "target_metadata.csv")
      ref_fc_file <- paste0(meta_dir, "ref_fc.txt")
      target_fc_file <- paste0(meta_dir, 'target_fc.txt')
      deseq2_dir <- paste0(proj_dir, "DESEQ2/", cfg_name, "/")  
      deseq2_conf <- paste0(meta_dir, "/deseq2_config.csv")
      gsea_dir <- paste0(proj_dir, "GSEA/", cfg_name, "/")  
      run_gsea_analysis()
    } # end if-statement
  } # end of config
} # end each sub_dir



