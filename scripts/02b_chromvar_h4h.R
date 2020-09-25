#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# run on h4h - CHANGE ONLY THE BOTTOM PART!
# Objective : Generic chromVAR script to run on any set of binary peakcount and annotation matrix and compare X conditions
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------
## load dependencies
#----------------------------------------------------------------
library(chromVAR)
#library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
#register(SerialParam())
#library(JASPAR2016)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(data.table)
library(gplots)
library(stringr)
source("scripts/chromvar.helper.R")

###############################################################




args <- commandArgs(trailingOnly = TRUE)

if (length(args)==4) {
  peaksbinarymat <- as.character(args[1])
  peaksbinarymatrepeats <- as.character(args[2])
  opname <- as.character(args[3])
  outdir <- as.character(args[4])

} else if (length(args)==0) {
  stop("Please provide the following arguments - Merged peak file, Directory containing all peak files, pattern to identify peak files,
       Output dir and Binary matrix file name", call.=FALSE)
}

###############################################################



## colors from devtools::install_github("caleblareau/BuenColors")
#library(BuenColors)
brewer_yes = c('#053061', '#2971B1', '#6AACD0', '#C1DDEB', '#F7F7F7', '#FACDB5', '#E58267', '#BB2933', '#67001F')
solar_basic = c('#214B85', '#1873CC', '#1E90FF', '#00BFFF', '#ACD8E5', '#D2D2D2', '#FFD700', '#ED2C2C', '#A31D1D')
Zissou = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")


run_chromvar=function(peakrds,
                       gsubpattern,
                       repeatrds,
                       opname,
                       outdir,
                       run_diff,
                       grouping,...){

    #----------------------------------------------------------------
    ## Load Countdata
    #----------------------------------------------------------------
    set.seed(2017)

    print("read in count data")
    my_counts_matrix <- readRDS(peakrds)
    rownames(my_counts_matrix) <- paste(my_counts_matrix[,1],my_counts_matrix[,2],my_counts_matrix[,3], sep="_")
    colnames(my_counts_matrix) <- gsub(gsubpattern, "",  colnames(my_counts_matrix))


    ## Remove random chromosomes

    toMatch <- c("random", "alt", "Un", "chrM", "EBV")
    my_counts_matrix <- subset(my_counts_matrix, !(grepl(paste(toMatch, collapse="|"), my_counts_matrix$seqnames)))

    fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix)
	#fragment_counts <- sort(fragment_counts, decreasing=FALSE)
    assayNames(fragment_counts) <- "counts"

    #----------------------------------------------------------------
    ## add gc content
    #----------------------------------------------------------------
    print("correcting for gc content")


    fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
    counts_filtered <- filterPeaks(fragment_counts,min_fragments_per_peak = 1, non_overlapping = TRUE)
    save(counts_filtered,file=paste0(outdir, "/chromvar/", opname, ".counts_filtered.Rdata"))
    rm(fragment_counts)

    # #----------------------------------------------------------------
    # ## Get motifs and what peaks contain motifs --
    # #----------------------------------------------------------------
    my_annotation_df <- readRDS(repeatrds)
    rownames(my_annotation_df) <- paste(my_annotation_df[,1], my_annotation_df[,2], my_annotation_df[,3], sep="_")
    my_annotation_df <- my_annotation_df[rownames(counts_filtered),]
    rm(my_counts_matrix)
    print("generate annotation database")

    anno_ix <- getAnnotations(as.matrix(my_annotation_df[,4:ncol(my_annotation_df)]), rowRanges = rowRanges(counts_filtered))
    save(anno_ix,file=paste0(outdir, "/chromvar/", opname, ".anno_ix.Rdata"))

    #----------------------------------------------------------------
    ## compute deviation
    #----------------------------------------------------------------
    set.seed(2017)
    print("Computing Deviation")
    dev <- computeDeviations(object = counts_filtered, annotations = anno_ix)
    save(dev,file=paste0(outdir, "/chromvar/", opname, ".devobj.Rdata"))

    z.scores = deviationScores(dev) ## deviation Z-score
    dev.scores = deviations(dev) ## bias corrected deviations

    write.table(z.scores, file=paste0(outdir, "/chromvar/", opname, ".Zscore.txt"), col.names=T, row.names=T, sep="\t", quote=F)
    write.table(dev.scores, file=paste0(outdir, "/chromvar/", opname, ".Deviations.txt"), col.names=T, row.names=T, sep="\t", quote=F)

    #----------------------------------------------------------------
    ## compute variablity
    #----------------------------------------------------------------
    variability <- computeVariability(dev)
    write.table(variability, file=paste0(outdir, "/chromvar/", opname, ".Variability.txt"), col.names=T, row.names=F, sep="\t", quote=F)

    #----------------------------------------------------------------
    ## Heatmap of repeats diff in CSC and TCGA
    #----------------------------------------------------------------

    if(run_diff==TRUE){

        diff_var <- differentialDeviations3(z.scores, grouping, parametric = FALSE)
        diff_var <- diff_var[order(diff_var$p_value_adjusted),]
        diff_var$repname <- rownames(diff_var)
        colnames(diff_var)[3:4] <- (levels(as.factor(grouping)))
        write.table(diff_var, file=paste0(outdir, "/chromvar/", opname, ".DiffDeviations.txt"), col.names=T, row.names=F, sep="\t", quote=F)

        top_motifs = subset(diff_var$repname,  diff_var$p_value_adjusted < 0.001)

        top_devs = z.scores[which(rownames(z.scores) %in% (top_motifs)), ]
        top_devs[!is.finite(top_devs)] <- NA
        top_devs[top_devs>10] <- 10
        top_devs[top_devs < -10] <- -10

        annocol1 <- data.frame(stringsAsFactors=F, group=grouping)
        rownames(annocol1) <- colnames(top_devs)

        pdf(paste0(outdir, "/chromvar/Heatmap", opname, ".DiffDeviations.qval1e-3.txt"), onefile=F)
        pheatmap(na.omit(top_devs),
                clustering_method="ward.D2",
                clustering_distance_rows="euclidean",
                clustering_distance_cols="euclidean",
                col = brewer_yes,
                annotation_col=annocol1,
                cluster_rows=T, cluster_cols=F,
                fontsize_col=5,fontsize_row=3,
                scale="none")
        dev.off()
    }
}


## In this part you need to enter the path to the two binary matrices and the directory where you binary matrices are located
## The directory containing the two binary matrices must contain a directory called "chromvar" where the output files will be generated

run_chromvar(peakrds=peaksbinarymat,
             gsubpattern="",
             repeatrds=peaksbinarymatrepeats,
             opname=opname,
             outdir=outdir,
             run_diff=F,
             grouping="")
