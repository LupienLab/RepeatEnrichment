#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
## Objective : function to analyse chromvar in different conditions
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------



analysechromvar=function( zscoresdf,dir,qval,annocol1,grouping,opname,dev1,...){
    diff_var <- differentialDeviations3(zscoresdf, grouping, parametric = FALSE)
    #colnames(diff_var)[3:4] <- levels(as.factor(grouping))
    diff_var <- diff_var[order(diff_var$p_value_adjusted),]
    diff_var$repname <- rownames(diff_var)
    write.table(diff_var, file=paste0(dir, "/plots_chromvar/", opname, ".DiffDeviations.with.medians.txt" ), col.names=T, row.names=F, sep="\t", quote=F)
    #diff_var <- diff_var[,c(1:2,5)]
    #write.table(diff_var, file=paste0(dir, "/plots_chromvar/", opname, ".DiffDeviations.txt" ), col.names=T, row.names=F, sep="\t", quote=F)

    top_motifs <- subset(diff_var$repname,  diff_var$p_value_adjusted < qval)
    top_devs <- as.matrix(zscoresdf[which(rownames(zscoresdf) %in% (top_motifs)), ])

    top_devs[!is.finite(top_devs)] <- NA
    #my.breaks <- c(seq(-20, 20, length.out=10))
    top_devs[top_devs>15] <- 15
    top_devs[top_devs < -15] <- -15
    cols <- colorRampPalette(brewer.pal(12, "Paired"))
    mycolors <- cols(length(unique(annocol1$Group)))
    names(mycolors) <- unique(annocol1$Group)
    mycolors2 <- cols(length(unique(annocol1$SampleName)))
    names(mycolors2) <- unique(annocol1$SampleName)
    mycolors <- list(Group=mycolors, SampleName=mycolors2)
    
    pdf(paste0(dir,"/plots_chromvar/Heatmap.",opname, ".DiffDeviations.qval", qval,".pdf" ), onefile=F)
    if(length(top_devs)>0){
    pheatmap(na.omit(top_devs),
    clustering_method="ward.D2",
    clustering_distance_rows="euclidean",
    clustering_distance_cols="euclidean",
    col = brewer_yes,
    annotation_col=annocol1,
    annotation_colors=mycolors,
    cluster_rows=T, cluster_cols=T,
    fontsize_col=5,fontsize_row=3,
    scale="none")
    }
    else{
     print(paste0(opname," does not pass threshold ",qval))  
    }
    dev.off()

    #breaks=my.breaks,
    # tsne_results <- deviationsTsne(dev1, threshold = 1.5, perplexity = 30)
    # colData(dev1)$group2 <- as.factor(annocol1$group.V2)
    # colData(dev1)$group3 <- as.factor(annocol1$group.V3)

    # tsne_plots <- plotDeviationsTsne(dev1, tsne_results,   annotation_name = NULL,
    #           sample_column = c("group2","group3"), shiny = FALSE)
    # pdf(paste0(dir, "/plots_chromvar/Tsne.",opname, ".pdf" ));
    # print(tsne_plots) ; dev.off()
}



repeatsig=function(opname, repnamelist,opname2,... ){

    load(paste0(dir,"/chromvar/",opname,".devobj.Rdata"))
    dev1 <- dev
    dev1_zscore <-  deviationScores(dev1) ## deviation Z-score
    dev1_metadata <- metadata[colnames(dev1_zscore),]

    dev1_annocol <- data.frame(stringsAsFactors=F,grp=dev1_metadata$V2,grp2=dev1_metadata$V3 )
    rownames(dev1_annocol) <- colnames(dev1_zscore)

    #rownames(dev1_zscore) <- str_split_fixed(rownames(dev1_zscore), "_",2)[,2]
    dev1_zscore_use <- dev1_zscore[repnamelist[which(repnamelist %in% rownames(dev1_zscore) )],]

    #dev1_zscore_use[dev1_zscore_use>10] <- 10
    #dev1_zscore_use[dev1_zscore_use<-10] <- -10
    my.breaks <- c(seq(-10, 10, length.out=10))

    cols <- colorRampPalette(brewer.pal(12, "Paired"))
    mycolors <- cols(length(unique(dev1_annocol$grp)))
    names(mycolors) <- unique(dev1_annocol$grp)
    mycolors2 <- cols(length(unique(dev1_annocol$grp2)))
    names(mycolors2) <- unique(dev1_annocol$grp2)
    mycolors <- list(grp=mycolors,grp2= mycolors2)

    pdf(paste0(dir,"/plots_chromvar/",opname2,".pdf"), onefile=F)
    pheatmap(na.omit(dev1_zscore_use),
    breaks=my.breaks,
    clustering_method="ward.D2",
    clustering_distance_rows="euclidean",
    clustering_distance_cols="euclidean",
    col = brewer_yes,
    annotation_col=dev1_annocol,
    annotation_colors=mycolors,
    cluster_rows=T, cluster_cols=T,
    fontsize_col=5,fontsize_row=3,
    scale="none")
    dev.off()

    pcaPRComp <- prcomp(t(as.matrix(na.omit(dev1_zscore_use))))
    df_out <- as.data.frame(pcaPRComp$x)
    df_out$group <- dev1_annocol$grp

    p <- ggplot(df_out,aes(x=PC1,y=PC2, color=as.factor(df_out$group))) +
    geom_text_repel(size =1,segment.size  = 0.2,aes(label = df_out$group)) +
    #theme(legend.position="none") +
    geom_point()

    pdf(paste0(dir,"/plots_chromvar/PCA.", opname2,".pdf"), onefile=F)
    print(p)
    dev.off()


    umap1 <- uwot::umap(t(na.omit(dev1_zscore_use)))
    umap1df <- as.data.frame(umap1)
    umap1df$group <- dev1_annocol$grp
    colnames(umap1df) <- c("PC1","PC2", "group")

    p1 <- ggplot(umap1df,aes(x=PC1,y=PC2, color=as.factor(umap1df$group))) +
    geom_text_repel(size =1,segment.size  = 0.2,aes(label = umap1df$group)) +
    #theme(legend.position="none") +
    geom_point()

    pdf(paste0(dir,"/plots_chromvar/Umap.", opname2,".pdf"), onefile=F)
    print(p1)
    dev.off()

}

