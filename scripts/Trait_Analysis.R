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

trait_data <- read.csv("/projectnb/tcwlab/LabMember/mwu/Project/Astrocyte_Ab/2X100/deseq2_conf/Metadata.csv", header=TRUE)
# Reformat the sample name
trait_data <- trait_data[, c(1, 2, 3,5,6,7,10, 11, 15)]

#Numberical variables
num_vars<-colnames(trait_data)[c(2,6,8,9)]

#Categorical variables
cat_vars<-colnames(trait_data)[c(3,4,5,7)]

for (i in cat_vars) {
  trait_data[, colnames(trait_data) == i] <- as.numeric(factor(trait_data[, colnames(trait_data) == i]))
}

#Change sample names to match
trait_data$Sample_Name <- gsub("-", ".", trait_data$Sample_Name)
traitRows_astro <- match(df[, 1], trait_data$Sample_Name)

trait_data_astro <- trait_data[traitRows_astro, -1]
nAstrocytes <- ncol(df) - 1
rownames(trait_data_astro) <- trait_data[traitRows_astro, 1]

# Correlate modules with traits using Spearman's correlation function:
moduleTraitCor_astro <- cor(module_merged$newMEs, trait_data_astro, use="p", method="spearman") # Obtain correlation value
moduleTraitPvalue_astro <- corPvalueStudent(moduleTraitCor_astro, nAstrocytes) # Obtain Student asymptotic p-values. 
textMatrix_astro_pvalues <- paste(signif(moduleTraitCor_astro, 2), "\n(", signif(moduleTraitPvalue_astro, 2), ")", sep="")
textMatrix_astro <- signif(moduleTraitCor_astro, 2)
dim(textMatrix_astro) <- dim(moduleTraitCor_astro)

pdf(paste0(outdir, "Trait_analysis_with_p_value.pdf"), width = 10, height = 10)
par(mar=c(9, 10, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor_astro,
               xLabels=names(trait_data_astro),
               yLabels=names(module_merged$newMEs),
               ySymbols=names(module_merged$newMEs),
               colorLabels=TRUE,
               colors=blueWhiteRed(30),
               textMatrix=textMatrix_astro_pvalues,
               setStdMargins=FALSE,
               cex.text=0.6,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()

# Create figure with * instead of p-values and no correlation values:
# p<0.001 '***', p<0.01 '**', p<0.05 '*', p<0.1 '.' 
moduleTraitPvalue_astro0 <- moduleTraitPvalue_astro
moduleTraitPvalue_astro0[moduleTraitPvalue_astro0 < 0.001] <- "***"
moduleTraitPvalue_astro0[moduleTraitPvalue_astro0 < 0.01 & moduleTraitPvalue_astro0 >= 0.001] <- "**"
moduleTraitPvalue_astro0[moduleTraitPvalue_astro0 < 0.05 & moduleTraitPvalue_astro0 >= 0.01] <- "*"
moduleTraitPvalue_astro0[moduleTraitPvalue_astro0 < 0.1 & moduleTraitPvalue_astro0 >= 0.05] <- "."
moduleTraitPvalue_astro0[moduleTraitPvalue_astro0 > 0.1] <- ""

pdf(paste0(outdir, "Trait_analysis_with_p_value.pdf"), width = 10, height = 10)
par(mar=c(9, 10, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor_astro,
               xLabels=names(trait_data_astro),
               yLabels=names(module_merged$newMEs),
               ySymbols=names(module_merged$newMEs),
               colorLabels=TRUE,
               colors=blueWhiteRed(30),
               textMatrix=moduleTraitPvalue_astro0,
               setStdMargins=FALSE,
               cex.text=1.4,
               zlim=c(-1,1),
               main=paste("Module-trait relationships with p_value"))
dev.off()