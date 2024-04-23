#-------------------------------------------------------------------------------------------
## ChromVAR anlysis script 
#-------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
## load dependencies
#-------------------------------------------------------------------------------------------
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2016)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gplots)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(attempt)

source("scripts/chromvar.helper.R")
source("scripts/chromvar_analyse_functionsB.R")

#################################################

args <- commandArgs(trailingOnly = TRUE)


if (length(args)>6) {
  metadatafile <- as.character(args[1])
  zscorefile <- as.character(args[2])
  devobjfile <- as.character(args[3])
  outdir <- as.character(args[4])
  outname <- as.character(args[5])
  qval <- as.numeric(args[6])
  groups <- args[7:length(args)]

} else {
  stop("Please provide the all the arguments", call.=FALSE)
}
#################################################

## colors from devtools::install_github("caleblareau/BuenColors")
#library(BuenColors)
brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')
solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D')
Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)

#-------------------------------------------------------------------------------------------
## Read in metadata
## NEEDS TO BE ADJUSTED FOR EACH ANALYSIS
#-------------------------------------------------------------------------------------------
# You can find an example of metadata in the directory containing the scripts
metadata <- read.table(metadatafile, header=T, stringsAsFactors=F)
#metadata$Group[which(metadata$Group=="PCa1" | metadata$Group=="PCaWZ1")] <- "grp1"
#metadata$Group[which(metadata$Group=="PCa2" | metadata$Group=="PCaWZ2")] <- "grp2"

#-------------------------------------------------------------------------------------------
## Here you should have a directory containing 2 subdirectories:
### one named "chromvar" containing the output files of the previous step
### and another called "chromvar_plots" where the output files will be generated
## Read in zscores data for samples
#-------------------------------------------------------------------------------------------
# read Z-scores and deviation object

zscore <- read.table(zscorefile,header=T, sep="\t", stringsAsFactors=F)
load(devobjfile)

annocol <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore) ,2:3])

# this first step take all samples into consideration
## I normally do not use these files but it gives me a general idea that everything is fine
# analyze chromvar
analysechromvar(zscoresdf=zscore,
                dir=outdir,qval=qval,
                annocol1=annocol,
                grouping=as.factor(annocol$group.Group),
                opname=outname,dev1=dev)
#-------------------------------------------------------------------------------------------
## Different comparisons - as many as you wish, by pairs
## This pat will focus the comparison between the groups of samples that you will indicate
# I set the qvalue at 0.01 to be more stringent but you may consider use 0.05 as well
#-------------------------------------------------------------------------------------------
# a) grp1 vs grp2

for (g in groups){
comp <- str_split(g,"_vs_",simplify = TRUE)
metadata_grp1_grp2 <- subset(metadata, metadata$Group == toString(comp[1]) | metadata$Group == toString(comp[2]))
zscore_grp1_grp2 <- zscore[,rownames(metadata_grp1_grp2)]
annocol <- metadata_grp1_grp2[,2:3]

print(paste0("processing ",g))

analysechromvar(zscoresdf=zscore_grp1_grp2,
                dir=outdir,qval=qval,
                annocol1=annocol,
                grouping=as.factor(annocol$Group),
                opname=g,dev1=dev)
}

# The comaprison will generate a heatmap with the most differential TEs between the two groups
# and a txt file containing all the TEs and 5 columns:
# p_value, p_value_adjusted, median_grp1 (Zscore), median_grp2 (Zscore), repname (TE name)
## I have a different script for heatmaps but you may consider using the one coming out of here
