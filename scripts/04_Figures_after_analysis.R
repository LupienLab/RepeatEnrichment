#-------------------------------------------------------------------------------------------
## ChromVAR anlysis script - on mordor loading R/3.4.1
#-------------------------------------------------------------------------------------------

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

# set working directory
setwd("/mnt/work1/users/lupiengroup/People/ggrillo/chromvar/")

source("scripts/chromvar.helperfunctions.R")
#source("scripts/chromvar_analyse_functionsB.R")
source("scripts/chromvar_analyse_functions.R")
not_all_na <- function(x) any(!is.na(x))


## colors from devtools::install_github("caleblareau/BuenColors")
#library(BuenColors)
brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')
solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D')
Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
scalered <- colorRampPalette(c("white", "red"), space = "rgb")(256)


#-------------------------------------------------------------------------------------------
## Volcano plot - using the txt file generated with script 03 for your groups of interest
#-------------------------------------------------------------------------------------------
grp1_grp2.diffdeviations <- read.table("path.to.diffDeviations.txt.file/chromvar_plots/fileName.DiffDeviations.with.medians.txt", header=T, sep="\t", stringsAsFactors=F)
grp1_grp2.diffdeviations$mediandiff <- grp1_grp2.diffdeviations$Stem - grp1_grp2.diffdeviations$Prostate
grp1_grp2.diffdeviations$neglog10qval <- -log10(grp1_grp2.diffdeviations$p_value_adjusted)
grp1_grp2.diffdeviations$type <- ifelse(grp1_grp2.diffdeviations$p_value_adjusted <= 0.01, "sig", "n.s.")
grp1_grp2.diffdeviations$Enrichment <- ifelse(grp1_grp2.diffdeviations$type=="sig" & grp1_grp2.diffdeviations$mediandiff>0, "grp1",
                                                 ifelse(grp1_grp2.diffdeviations$type=="sig" & grp1_grp2.diffdeviations$mediandiff<0,"grp2", "n.s."))

pdf("path.to.wanted.directory/fileName.pdf", useDingbats=F)
ggplot(data=grp1_grp2.diffdeviations, aes(x=mediandiff, y=neglog10qval)) +
  geom_point(aes(color=Enrichment)) +  
  theme_bw() +
  labs(x = "Median Zscore difference (grp1 - grp2)" , y="-log10 (q-value)", colour = "Enrichment") + 
  scale_colour_manual(values=c(grp1="#CC0033",n.s.="#000000", grp2="#1e71aa")) +
  geom_hline(yintercept=2,linetype="dashed", size=0.5, col="grey") +
  xlim(-15, 10) +
  theme(legend.position = c(1, 1), legend.background=element_rect(size=0.2, linetype="solid", colour ="#000000")) +
  theme_bw()
dev.off()


#-------------------------------------------------------------------------------------------
## Plot TE families for grp1 and grp2
#-------------------------------------------------------------------------------------------
grp1.grp1_grp2 <- fread("path.to.diffDeviations.txt.file/chromvar_plots/fileName.DiffDeviations.with.medians.txt", header=T, sep="\t", stringsAsFactors=F)

grp1.grp1_grp2.reps <- subset(grp1.grp1_grp2$repname,grp1.grp1_grp2$p_value_adjusted<=0.01 & 
                                            grp1.grp1_grp2$median_grp1 > grp1.grp1_grp2$median_grp2)

repeat_metadata <- read.table("data/repeat.subfam/repeat_metadata.onlyTEs.txt", header=T, sep="\t", stringsAsFactors=F)
repeat_metadata$grp1.grp1_grp2.reps <- ifelse(repeat_metadata$repname %in% grp1.grp1_grp2.reps,1,0)
#repeat_metadata$lsc.lsc_csc_cscneg_diff.reps <- ifelse(repeat_metadata$repname %in% lsc.lsc_csc_cscneg_diff.reps,1,0)

#repeat_metadata <- repeat_metadata[rowSums(repeat_metadata[,3:ncol(repeat_metadata)])>0,]

colnames(repeat_metadata)[3:ncol(repeat_metadata)] <- gsub(".reps","", colnames(repeat_metadata)[3:ncol(repeat_metadata)] )
colnames(repeat_metadata) <- c("repname", "repname2","Stem_Prostate")
rownames(repeat_metadata) <- repeat_metadata$repname2
repfam <- read.table("data/repeat.subfam/repeat_fam_mapping.onlyTEs.txt", header=F , sep="\t", stringsAsFactors=F)
repeat_metadata <- merge(repeat_metadata, repfam[,1:2], by.x="repname2",by.y="V2", all.x=T)
rownames(repeat_metadata) <- repeat_metadata$repname2
#repeat_metadata[is.na(repeat_metadata)] <- "UNKNOWN"

# merge repeat data with z-scores and stats'
repfam_stats <- merge(grp1_grp2.diffdeviations, repeat_metadata, by="repname", all=TRUE)
repfam_stats <- na.omit(repfam_stats)
repfam_sig <- repfam_stats
repfam_sig$sig <- ifelse(repfam_sig$type == "sig", 1,0)
repfam_sig$siggrp1 <- ifelse(repfam_sig$Enrichment == "grp1", 1,0)
repfam_sig$siggrp2 <- ifelse(repfam_sig$Enrichment == "grp2", 1,0)
summary_families <- plyr::ddply(repfam_sig, "V1", summarise, count=length(type), sig=sum(sig), grp1 =sum(siggrp1), grp2=sum(siggrp2), .drop=FALSE)
summary_families[summary_families$sig == 0,][,1]

#sigstem and signonstem repeats
grp1_repeats <- subset(repfam_sig, repfam_sig$siggrp1==1)[,1]
grp2_repeats <- subset(repfam_sig, repfam_sig$siggrp2==1)[,1]

#"LINE3" "LINE5" "rRNA"  "SINE"  "SINE3"
#repfam_sig <- repfam_sig[!(repfam_sig$V1 %in% c("LINE3", "LINE5", "rRNA", "SINE", "SINE3")),]
repfam_sig$famSort <- factor(repfam_sig$V1, )

# plot this wonderful overview plot 
pdf(paste0("path.to.wanted.directory/fileName.pdf" ), onefile=F, useDingbats=FALSE)
ggplot(repfam_sig,aes(x=famSort,y=mediandiff))+
geom_point(aes(colour=Enrichment, size=neglog10qval))+coord_flip()+theme_bw() + 
theme(legend.position = c(0.8, 0.2)) +
labs(x = "Repeat Family" , y="Median Zscore difference (grp1 - grp2)", colour = "Enrichment q-value <= 0.01") +
scale_colour_manual(values=c(grp1="#CC0033",n.s.="#000000", grp2="#1e71aa")) +
scale_x_discrete(limits = rev(levels(repfam_sig$famSort))) 
dev.off()

# plot this version as well with size proportional to -log10(qvalue) 
pdf(paste0("path.to.wanted.directory/fileName.pdf" ), onefile=F, useDingbats=FALSE)
ggplot(repfam_sig,aes(x=famSort,y=mediandiff))+
geom_point(aes(colour=Enrichment))+coord_flip()+theme_bw() + 
theme(legend.position = c(0.8, 0.1), legend.background=element_rect(size=0.2, linetype="solid", colour ="#000000")) +
labs(x = "Repeat Family" , y="Median Zscore difference (grp1 - grp2)", colour = "Enrichment q-value <= 0.01") +
scale_colour_manual(values=c(grp1="#CC0033",n.s.="#000000", grp2="#1e71aa")) +
scale_x_discrete(limits = rev(levels(repfam_sig$famSort)))
dev.off()


#-------------------------------------------------------------------------------------------
## Plot box plot for TE of interest
## this can be used to generate boxplots for two or more groups
#-------------------------------------------------------------------------------------------

metadata <- read.table("path.to.same.metadata.used.in.previous.step/fileName.txt", header=T, stringsAsFactors=F)
metadata <- unique(metadata)

zscore_file <- read.table("path.to.Zscore.file/chromvar/fileName.Zscore.txt",header=T, check.names=F,sep="\t", stringsAsFactors=F)
annocol_file <- data.frame(stringsAsFactors=F, group=metadata[colnames(zscore_file) ,2:3])


zscore_anno <- t(zscore_file)
zscore_anno <- merge(annocol_file, zscore_anno, by=0, all=TRUE)
rownames(zscore_anno) <- zscore_anno[,1]
zscore_anno <- zscore_anno[,-1]
zscore_anno <- zscore_anno[order(zscore_anno$group.Group),]
#zscore_anno <- unique(zscore_anno)

zscore_anno$grps <-factor(zscore_anno$group.Group, levels = c("grp1","grp2","grp3","grp4"), ordered=T)

pdf(paste0("path.to.wanted.directory/fileName.pdf" ), onefile=F)
p <-ggplot((zscore_anno=subset(zscore_anno, !is.na(zscore_anno$grps))), aes(x=zscore_anno$grps, y=zscore_anno$`TE.subfamily.name`)) + 
geom_boxplot(aes(fill=zscore_anno$grps)) + 
coord_flip() + 
scale_fill_manual(values=c(grp1= "#CC0033", grp2 = "#1e71aa", grp3 = "#013220", grp4 = "#68228B")) +
#theme_bw() + 
theme(legend.position = "right") +
geom_jitter(color="black", size=0.4, alpha=0.9) +
labs(title="TE.subfamily.name - grp#", x=element_blank() , y="Zscore", fill = "Group")
p
dev.off()

