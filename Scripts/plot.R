rankHeatmap2 <- function(res, rankBreaks=c(1,10,30,75,150,300,600), ..., cap=NULL,
                         column_title="Datasets", doDraw=TRUE, rowann=NULL, colann=NULL,
                         datasetInfo=data.frame(
                          type=sapply(getDatasets(), FUN=function(x) x$type)),
                         cellLabelFontsize=8, row_order=NULL, column_order=NULL){
  library(ComplexHeatmap)
  rankm <- reshape2::dcast(res, method~dataset, value.var="rank")
  row.names(rankm) <- rankm[,1]
  rankm <- sqrt(as.matrix(rankm[,-1]))
  rankmImputed <- .getImputed(rankm)
  if(!is.null(cap)) rankm[which(rankm>cap)] <- cap
  if(is.null(row_order)) 
    row_order <- row.names(rankm)[order(rowMeans(rankmImputed, na.rm=TRUE))]
  if(is.null(column_order))
    column_order <- names(sort(-colMeans(rankm, na.rm=TRUE)))
  ancols <- list(type=c(CRISPRi="darkorange3", dTag="black", ligand="darkslateblue"))
  if(is.null(colann)){
    colann <- HeatmapAnnotation(df = datasetInfo[colnames(rankm),,drop=FALSE],
                                col=ancols, show_annotation_name=FALSE)
  }else if(all(is.na(colann))){
    colann <- NULL
  }
  bdist <- log10(rankBreaks[-1]-rankBreaks[-length(rankBreaks)])
  rann <- getMethodAnno(row.names(rankm))
  if(is.null(rowann)){
    rowann <- rowAnnotation(df=rann$df, col=rann$col,
                            show_legend=colnames(rann$df)=="family")
  }else if(all(is.na(rowann))){
    rowann <- NULL
  }
  h <- ComplexHeatmap::Heatmap(rankm, cluster_rows=FALSE, cluster_columns=FALSE,
          cell_fun=function(j, i, x, y, width, height, fill) {
            l <- sprintf("%1.0f", rankm[i, j]^2)
            if(!is.null(cap) && rankm[i,j]==cap) l <- ""
            grid.text(l, x, y,
                      gp = gpar(fontsize=ifelse(is.na(rankm[i, j]),6,
                                                cellLabelFontsize)))
          },
     col=viridisLite::viridis(100, direction=-1), name="Rank of\ntrue TF", ...,
     column_title=column_title, row_order=row_order, column_order=column_order, 
     left_annotation=rowann, row_names_side="left", top_annotation=colann,
     heatmap_legend_param=list(at=rev(sqrt(rankBreaks)),
                               labels=rev(c("top",rankBreaks[-1])),
                               break_dist=rev(bdist), legend_height=unit(3.5, "cm"),
                               labels_gp=gpar(fontsize=9)))
  if(!doDraw) return(h)
  draw(h, merge=TRUE)
}


sensFDRplot <- function(res, fade=NULL, PR=TRUE, hull=TRUE, label.size=3.5, longTitles=TRUE, useArch=FALSE){
  res$Sensitivity <- res$trueQ<0.05
  res$FDR[which(is.na(res$FDR))] <- 0
  if(!is.null(res$archFDR)) res$archFDR[which(is.na(res$archFDR))] <- 0
  if(is.null(res$type)){
    res$type <- ifelse(grepl("GSEA|ulm|msViper|decoupleR|MonaLisa|diffTF|-lm|Lasso|meirlop",
                             res$method, ignore.case=TRUE), "LFC-based", "sample-wise")
    res$type[grep("limma",res$method)] <- "sample-wise"
  }
  cols <- setNames(c("#CC6677", "#4477AA"), unique(res$type))
  res2 <- aggregate(res[,intersect(c("Sensitivity", "FDR","archFDR"),colnames(res))],
                    by=res[,c("method","type")], na.rm=TRUE, FUN=mean)
  if(useArch){
    res2$precision <- 1-res2$archFDR
  }else{
    res2$precision <- 1-res2$FDR
  }
  res2$fade <- FALSE
  if(!is.null(fade)) res2$fade <- grepl(fade, res2$method)
  res3 <- res2[chull(res2$Sensitivity, res2$precision),]
  res3 <- res3[which(res3$Sensitivity==max(res3$Sensitivity) | 
                       res3$precision==max(res3$precision) |
                       (res3$Sensitivity+res3$precision)>
                         mean(res2$Sensitivity+res2$precision)),]
  if(PR){
    p <- ggplot(res2, aes(Sensitivity, precision, label=method, colour=type, alpha=fade))
    if(longTitles){
      yl <- ifelse(useArch,
                   "Precision\n(Proportion in network archetypes)",
                   "Precision\n(Proportion of network motifs among significant)"
                   )
      p <- p +
      labs(x="Recall (i.e. sensitivity)\n(Proportion of datasets in which the true motif is significant)",
           y=yl)
    }else{
      p <- p +
        labs(x="Recall (i.e. sensitivity)",
             y="Precision")
    }
  }else{
    p <- ggplot(res2, aes(FDR, Sensitivity, label=method, colour=type, alpha=fade))
  }
  if(hull) p <- p + geom_line(data=res3, colour="lightgrey", linetype="dashed", alpha=1)
  p +
    geom_point(show.legend=FALSE) +  theme_bw() + 
    scale_x_continuous(breaks=scales::pretty_breaks()) + 
    scale_y_continuous(breaks=scales::pretty_breaks()) + 
    ggrepel::geom_text_repel(min.segment.length=0, show.legend=FALSE, size=label.size) +
    scale_color_manual(values=cols) + scale_alpha_manual(values=c("TRUE"=0.5, "FALSE"=1)) +
    theme(legend.position = "bottom")
}

relAUCplot <- function(res, squeeze=FALSE){
  res$auc2 <- gsub("^0\\.",".",round(res$relAUC,2))
  p <- ggplot(res, aes(reorder(dataset,relAUC), reorder(method,relAUC), 
                  fill=relAUC, label=auc2)) + 
    geom_tile() + theme_minimal() + geom_text() +
    scale_fill_viridis_c(option="A") +
    labs(fill="relative\nAUC", y="Methods", x="Datasets") + theme_minimal()
  if(squeeze) p <- p + 
    theme(legend.position="bottom", legend.key.width = unit(1,"cm"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p
}

.getImputed <- function(x){
  imput <- pmin(colMeans(x,na.rm=TRUE),
                matrixStats::colMedians(x,na.rm=TRUE))
  rankmImputed <- x
  for(i in seq_len(ncol(x)))
    rankmImputed[which(is.na(x[,i])),i] <- imput[i]
  rankmImputed
}

getOrder <- function(res, values=c("transRank","relAUC","archRelAUC"),
                     columns=FALSE){
  x <- -sqrt(res$rank)
  res$transRank <- 2*sqrt(exp(x)/(1+exp(x)))
  y <- lapply(values, FUN=function(x){
    x <- as.data.frame(reshape2::dcast(res, method~dataset, value.var=x))
    row.names(x) <- x[,1]
    .getImputed(as.matrix(x[,-1]))
  })
  if(columns){
    y <- do.call(rbind, y)
    return(colnames(y)[order(colMeans(y))])
  }
  y <- do.call(cbind, y)
  row.names(y)[order(rowMeans(y))]
}

relAUCplot2 <- function(res, doDraw=TRUE, column_title="Datasets", name=NULL,
                        column_order=NULL, row_order=NULL, ..., val="relAUC", col="AUC",
                        type=c("propOfMax","MAD","relAUC"), cellLabelFontsize=8,
                        hmcolors=viridisLite::magma(100),
                        datasetInfo=data.frame(
                           type=sapply(getDatasets(), FUN=function(x) x$type))){
  type <- match.arg(type)
  m1 <- reshape2::dcast(res, method~dataset, value.var=val)
  if(type=="propOfMax"){
    dmax <- aggregate(res[[col]], by=list(ds=res$dataset), na.rm=TRUE, FUN=max)
    dmax <- setNames(dmax$x, dmax$ds)
    res$AUCprop <- res[[col]]/dmax[as.character(res$dataset)]
    m <- reshape2::dcast(res, method~dataset, value.var="AUCprop")
  }else if(type=="MAD"){
    dmed <- aggregate(res[[val]], by=list(ds=res$dataset), na.rm=TRUE, FUN=median)
    dmed <- setNames(dmed$x, dmed$ds)
    res$AUCdiff <- res[[val]]-dmed[as.character(res$dataset)]
    res$AUCdiff <- res$AUCdiff/median(abs(res$AUCdiff),na.rm=TRUE)
    m <- reshape2::dcast(res, method~dataset, value.var="AUCdiff")
  }else{
    m <- m1
  }
  row.names(m1) <- row.names(m) <- m[,1]
  m <- as.matrix(m[,-1])
  m1 <- as.matrix(m1[,-1])
  if(is.null(column_order)) column_order <- names(sort(colMeans(m, na.rm=TRUE)))
  if(is.null(row_order)) row_order <- row.names(m)[order(-rowMeans(m, na.rm=TRUE))]

  ancols <- list(type=c(CRISPRi="darkorange3", dTag="black", ligand="darkslateblue"),
                 LFCbased=c("FALSE"="white", "TRUE"="brown4"))
  colan <- HeatmapAnnotation(df = datasetInfo[colnames(m),,drop=FALSE],
                             col=ancols, show_annotation_name=FALSE)
  if(is.null(name))
    name <- ifelse(type=="relAUC","Relative\nnetwork AUC","Relative\nnetwork\nscore")
  h <- ComplexHeatmap::Heatmap(
      m, cluster_rows=FALSE, cluster_columns=FALSE, ..., column_order=column_order,
      cell_fun=function(j, i, x, y, width, height, fill) {
        v <- gsub("^0","",round(m1[i,j],2))
        if(is.na(v)) v <- "NA"
        if(v=="") v <- "0"
        if(v=="NaN") v <- "NA"
        grid.text(v, x, y, gp=gpar(fontsize=ifelse(is.na(m1[i,j]),6,
                                                   cellLabelFontsize),
                                   col=ifelse(m[i,j]<0.5 & !is.na(m[i,j]),
                                              "lightgrey","black")))
      },
      col=hmcolors, name=name,
      column_title=column_title, row_names_side="left", top_annotation=colan)
  if(!doDraw) return(h)
  draw(h, merge=TRUE)
}

runtimePlot <- function(res, removeNAs=FALSE, cpu=TRUE){
  if(removeNAs) res <- res[!is.na(res$elapsed) & !is.na(res$cpu),]
  if(cpu){
    res$runtime <- pmax(res$elapsed, res$cpu, na.rm=TRUE)
  }else{
    res$runtime <- res$elapsed
  }
  ggplot(res, aes(runtime/60, reorder(method,runtime))) + stat_summary() +
    #geom_boxplot(fill="lightgrey", outlier.size=1) + 
    theme_bw() + labs(x=ifelse(cpu,"CPU time (min)","Total elapsed time (min)"),y="")
}

runtimePlot2 <- function(res, removeNAs=FALSE, removeCPU="Lasso|decoupleR|VIPER"){
  if(removeNAs) res <- res[!is.na(res$elapsed) & !is.na(res$cpu),]
  if(!is.null(removeCPU)){
    w <- which(grepl(removeCPU,res$method))
    res$elapsed[w] <- pmax(res$elapsed[w], res$cpu[w], na.rm=TRUE)
    res$cpu[w] <- NA
  }
  res3 <- reshape2::melt(res, id.vars=c("method","dataset"),
                         measure.vars=c("elapsed","cpu"))
  ggplot(res3, aes(value/60, reorder(method,value,na.rm=TRUE), colour=variable)) + 
    stat_summary(position=position_dodge(width=0.25)) +
    theme_bw() + labs(x="Time per dataset (min)", y="") +
    theme(axis.text.y=element_text(color="black"))
  
}

plotAllMetrics <- function(res){
  cowplot::plot_grid(
    rankHeatmap(res, squeeze = TRUE),
    relAUCplot(res, squeeze = TRUE),
    sensFDRplot(res), runtimePlot(res),
    nrow=2, labels="AUTO", scale=0.95, rel_heights = c(5,4)
  )
}

renameMethods <- function(x, renaming=NULL){
  if(is.data.frame(x)){
    stopifnot("method" %in% colnames(x))
    x$method <- renameMethods(x$method, renaming=renaming)
    return(x)
  }
  if(is.null(renaming)) renaming <- c(
    "CVoriginal"="chromVAR::differentialDeviations", "CV"="chromVAR(z)>limma",
    "CVqt"="chromVAR(z)>Qt>limma", "CVcentered"="chromVAR(z)>center>limma",
    "CVnorm"="chromVAR(z)>scale>limma", MLzero="monaLisa.vsZero",
    MLsp="monaLisa.vsOthers+spearman", MLfewerBins="monaLisa.vsOthers(fewer bins)",
    MLStabSel="monaLisa.StabSel", MLlower="monaLisa.vsOthers(smaller zeroBin)",
    ML="monaLisa.vsOthers", "regreg"="Lasso-lm", "regregR"="Ridge-lm",
    VIPER="VIPER(scores)>limma", VIPERb="VIPER(binary)>limma",
    msVIPER="msVIPER(scores)", msVIPERb="msVIPER(binary)", ulmGC="ulm+GC", 
    ulmB="ulm(binary)", ulm="ulm(scores)", MBA="insertionModel", GSEA="fGSEA",
    CVdev="chromVAR(deviations)>limma", CVdevCentered="chromVAR(deviations)>center>limma",
    CVdevNorm="chromVAR(deviations)>scale>limma", CVdevqt="chromVAR(deviations)>Qt>limma",
    diffTF="diffTF(permutations)", diffTF_noPerm="diffTF(analytic)",
    fastMLM="GC>fastMLM>limma", meirlop="MEIRLOP"
  )
  for(i in names(renaming)) x <- replace(x, x==i, renaming[[i]])
  x <- gsub("decoupleR","decoupleR(",x,fixed=TRUE)
  w <- grepl("decoupleR",x,fixed=TRUE)
  x[w] <- paste0(x[w],")")
  w <- grepl("decoupleR(limma",x,fixed=TRUE)
  x <- gsub("decoupleR(limma","decoupleR(",x,fixed=TRUE)
  x[w] <- paste0(x[w],">limma")
  x <- gsub("::",":",x,fixed=TRUE)
  x <- gsub("((","(",x,fixed=TRUE)
  x <- gsub("))",")",x,fixed=TRUE)
  x <- gsub(">limma)",">limma",x,fixed=TRUE)
  x
}

filterMethods <- function(x, meth=getMainMethods(), rename=FALSE){
  w <- which(x$method=="CV" | x$method=="chromVAR>scale>limma")
  CVcpu <- setNames(x$cpu[w], x$dataset[w])
  CVelapsed <- setNames(x$elapsed[w], x$dataset[w])
  w <- which(grepl("^CV|^chromVAR",x$method) & is.na(x$cpu))
  x$cpu[w] <- CVcpu[x$dataset[w]]
  x$elapsed[w] <- CVelapsed[x$dataset[w]]
  x$method <- as.character(x$method)
  renames <- meth
  if(rename) x <- renameMethods(x)
  x <- x[which(x$method %in% names(renames)),]
  renameMethods(x, renames)
}

getTopMethods <- function(){
  c(`chromVAR::differentialDeviations` = "chromVAR", 
    `chromVAR(z)>Qt>limma` = "chromVAR-adjusted",
    monaLisa.vsOthers = "monaLisa", monaLisa.StabSel = "monaLisa.StabSel", 
    `msVIPER(scores)` = "msVIPER", `VIPER(binary)>limma` = "VIPER", 
    `decoupleR:consensus` = "decoupleR", fGSEA = "fGSEA")
}

getMainMethods <- function(){
  c(`chromVAR:differentialDeviations` = "chromVAR", 
    `chromVAR(z)>Qt>limma` = "chromVAR-adjusted", diffTF="diffTF",
    monaLisa.vsOthers = "monaLisa", monaLisa.StabSel = "monaLisa.StabSel", 
    `msVIPER(scores)` = "msVIPER", `VIPER(binary)>limma` = "VIPER", 
    `decoupleR(consensus)` = "decoupleR(consensus)", `decoupleR(mlm)>limma` = "decoupleR(mlm)>limma", 
    fGSEA = "fGSEA", insertionModel="insertionModel", ulmB="ulm",
    "Lasso-lm"="Lasso", "BaGFootLike"="BaGFootLike",
    "diffTF(analytic)"="diffTF(analytic)", MEIRLOP="MEIRLOP",
    "diffTF(permutations)"="diffTF(permutations)")
}


getMethodAnno <- function(methods){
  methods <- unique(methods)
  LFCbased <- grepl("GSEA|ulm|msViper|decoupleR:|MonaLisa|diffTF|meirlop|-lm",
                    methods,ignore.case=TRUE) &
    !grepl("decoupleR.+limma",methods,ignore.case=TRUE)
  family <- grepl("chromVAR",methods) + 10*grepl("monaLisa",methods) +
    100*grepl("decoupleR",methods) + 1000*grepl("VIPER",methods)
  family <- factor(family, c(1,10,100,1000,0),
                   c("chromVAR", "monaLisa", "decoupleR", "VIPER", "other"))
  ancols <- list(type=c(CRISPRi="darkorange3", dTag="black", ligand="darkslateblue"),
                 family=setNames(c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
                                 levels(family)),
                 LFCbased=c("FALSE"="white", "TRUE"="brown4"))
  df <- data.frame(row.names=methods, LFCbased=LFCbased, family=family)
  if(any(grepl("^>GC", methods))){
    df$GCsmoothQ <- grepl("^GC>",methods)
    ancols$GCsmoothQ <- ancols$LFCbased
  }
  list(df=df, col=ancols)
}

resHeatmap <- function(res, row_names_gp=gpar(fontsize=10), ...){
  res$relAUC[is.na(res$relAUC)] <- ifelse(is.na(res$rank[which(is.na(res$relAUC))]), NA, 0)
  res$archRelAUC[is.na(res$archRelAUC)] <- ifelse(is.na(res$rank[which(is.na(res$archRelAUC))]), NA, 0)
  ro <- rev(getOrder(res))
  co <- getOrder(res, columns=TRUE)
  h1 <- rankHeatmap2(res, doDraw=FALSE, column_title="Rank of true TF", 
                     row_names_gp=row_names_gp, row_order=ro, column_order = co, ...)
  h2 <- suppressWarnings(relAUCplot2(res, row_order=ro, doDraw=FALSE,
                                     column_title="Network score", column_order=co, ...))
  h3 <- suppressWarnings(relAUCplot2(res, row_order=ro, doDraw=FALSE, name="Archetype\nscore",
                                     column_title="Archetype score", val="archRelAUC", ...,
                                     col="archRelAUC", column_order=co, hmcolors=viridisLite::inferno(100)))
  h1 + h2 + h3
}
