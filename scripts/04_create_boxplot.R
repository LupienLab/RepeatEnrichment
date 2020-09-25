library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gplots)
library(stringr)
library(RColorBrewer)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)>5) {
  metadatafile <- as.character(args[1])
  zscorefile <- as.character(args[2])
  outdir <- as.character(args[3])
  outname <- as.character(args[4])
  groups <- args[5:length(args)]

} else {
  stop("Please provide the all the arguments", call.=FALSE)
}

metadata <- read.table(metadatafile, header=T, stringsAsFactors=F)
metadata <- unique(metadata)

zscore <- read.table(zscorefile, header=T, sep="\t", stringsAsFactors=F)
annocol_file <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore) ,2:3])


zscore_anno <- t(zscore)
zscore_anno <- merge(annocol_file, zscore_anno, by=0, all=TRUE)
rownames(zscore_anno) <- zscore_anno[,1]
zscore_anno <- zscore_anno[,-1]
zscore_anno <- zscore_anno[order(zscore_anno$group.Group),]


zscore_anno$grps <-factor(zscore_anno$group.Group, levels = c("grp1","grp2","grp3","grp4"), ordered=T)
print(head(zscore_anno))
pdf(paste0(outdir,outname,"_boxplot.pdf"), onefile=F)
p <-ggplot((zscore_anno=subset(zscore_anno, !is.na(zscore_anno$grps))), aes(x=zscore_anno$grps, y=zscore_anno$`TE.subfamily.name`)) +
geom_boxplot(aes(fill=zscore_anno$grps)) +
coord_flip() +
scale_fill_manual(values=c(grp1= "#CC0033", grp2 = "#1e71aa", grp3 = "#013220", grp4 = "#68228B")) +
theme(legend.position = "right") +
geom_jitter(color="black", size=0.4, alpha=0.9) +
labs(title="TE.subfamily.name - grp#", x=element_blank() , y="Zscore", fill = "Group")
p
dev.off()

