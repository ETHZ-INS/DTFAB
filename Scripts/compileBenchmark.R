#' Get benchmark metrics for all dataset
#'
#' @param datasets A list of dataset objects (see getDatasets() )
#' @param path The path in which the datasets' results folders are
#' @param resin The name of the results subfolder
#' @param interactors A named list of interactors per TF
#'
#' @return A data.frame of compiled metrics
compileBenchmark <- function(datasets, rootPath=".", resin="runATAC_results",
                             interactors, archetypes){
  res <- list()
  for(dn in names(datasets)){
    ds <- datasets[[dn]]
    if(!is.null(ds$folder)) dn <- ds$folder
    archs <- lapply(archetypes, FUN=function(archs){
      setNames(strsplit(names(archs),"/"), names(archs))
    })
    print(dn)
    tryCatch(
      res[[dn]] <- getBenchmarkMetrics(ds, file.path(rootPath, dn),
                                       interactors=interactors, archetypes=archs),
      error=function(e){
        message(e);
        #traceback()
      })
  }
  dplyr::bind_rows(res, .id="dataset")
}


#' Get benchmark metrics for one dataset
#'
#' @param dataset A dataset object (see getDatasets() )
#' @param path The path in which the dataset's results are
#' @param resin The name of the results subfolder
#' @param interactors A named list of interactors per TF
#' @param archetypes An optional list of archetypes for each species. Each list
#'  should itself be a list of archetypes, with each archetype being
#'   a character vector of the motifs in it.
#'
#' @return A data.frame of compiled metrics for one dataset
getBenchmarkMetrics <- function(dataset, path=head(dataset$truth,1),
                                resin="runATAC_results", interactors,
                                archetypes=NULL){
  stopifnot(is.null(archetypes) || is.list(archetypes))
  stopifnot(is.list(interactors))
  # grab runtimes
  fl <- list.files(file.path(path, resin, "raw"), pattern="\\.rds$", full=TRUE)
  names(fl) <- gsub("_raw\\.rds$","",basename(fl))
  rt <- sapply(fl, FUN=function(x){
    x <- readRDS(x)
    if(!is.null(names(x)) && "runtime" %in% names(x)){
      x <- x$runtime
    }else{
      x <- x[[2]]
    }
    if(isTRUE(try(is(x,"proc_time"),silent=TRUE))) x <- summary(x)
    x
  })
  rt <- rbind(elapsed=rt["elapsed",], cpu=colSums(rt[which(row.names(rt)!="elapsed"),]))
  rt <- as.data.frame(t(rt))
  if(any(row.names(rt)=="decoupleR")){
    # grab times for individual methods
    dt <- readRDS(file.path(path, resin, "raw","decoupleR_raw.rds"))$runtime2
    rt <- rbind(rt, data.frame(row.names=paste0("decoupleR",names(dt)),
                               elapsed=dt, cpu=dt))
  }
  if(any(row.names(rt)=="decoupleRlimma")){
    # grab times for individual methods
    dt <- readRDS(file.path(path, resin, "raw","decoupleRlimma_raw.rds"))$runtime2
    rt <- rbind(rt, data.frame(row.names=paste0("decoupleRlimma",names(dt)),
                               elapsed=dt, cpu=dt))
  }

  # get actual results
  fl <- list.files(file.path(path, resin, "with_pvalues"), full=TRUE)
  names(fl) <- gsub("\\.rds$","",basename(fl))
  res <- lapply(fl, readRDS)
  res <- res[sapply(res, FUN=function(x) isTRUE(nrow(x)>1))]
  
  res <- dplyr::bind_rows(lapply(res, FUN=function(x){
    m <- .getBenchmarkMetrics(x, truth=dataset$truth, interactors=interactors)
    if(!is.null(archetypes)){
      arch <- archetypes[[dataset$species]]
      arch <- unlist(arch[which(lengths(lapply(arch, y=dataset$truth,
                                               FUN=intersect))>0)])
      m2 <- .getBenchmarkMetrics(x, truth=dataset$truth, interactors=arch)
      m$archRelAUC <- m2$relAUC
      m3 <- .getArchMetrics(x, dataset$truth, cofactors=interactors,
                            archetypes=archetypes[[dataset$species]])
      cbind(m, m3)
    }
  }), .id="method")
  res$elapsed <- rt[res$method, "elapsed"]
  res$cpu <- rt[res$method, "cpu"]
  res
}

.getArchMetrics <- function(x, truth, archetypes, cofactors){
  if(is.list(cofactors))
    cofactors <- unique(unlist(cofactors[intersect(truth,names(cofactors))]))
  tf2arch <- setNames(unlist(archetypes), names(archetypes))
  x$arch <- tf2arch[row.names(x)]
  w <- which(is.na(x$arch))
  x$arch[w] <- row.names(x)[w]
  trueArch <- unique(x$arch[which(row.names(x) %in% c(truth))])
  notFalseArch <- unique(x$arch[which(row.names(x) %in% c(truth,cofactors))])
  x$trueArch <- x$arch %in% trueArch
  x$notFalseArch <- x$arch %in% notFalseArch
  w <- head(which(x$trueArch),1)
  xsig <- x[which(x$padj<0.05),,drop=FALSE]
  if(length(w)==0){
    w <- nrow(x)+1L
    return(data.frame(archRank=w, archSens=0, archFDR=1))
  }
  data.frame(archRank=w, archSens=as.numeric(x[w,"padj"]<0.05),
             archFDR=sum(!xsig$notFalseArch)/nrow(xsig))
}

.getBenchmarkMetrics <- function(x, truth, interactors){
  if(is.list(interactors)){
    cofactors <- unique(c(truth, unlist(interactors[truth])))
  }else{
    cofactors <- union(truth, interactors)
  }
  allTFs <- row.names(x)
  cofactors <- intersect(cofactors, allTFs)
  cappedNCof <- min(length(cofactors),100)
  optimalAUC <- sum(cumsum(c(rep(1,cappedNCof),
                             rep(0,100-cappedNCof)))/seq_len(100))
  w <- head(which(row.names(x) %in% truth),1)
  if(length(w)==0){
    w <- nrow(x)+1L
    trueQ <- 1
  }else{
    trueQ <- x[w,"padj"]
  }
  x$isCofactor <- row.names(x) %in% cofactors
  tmp <- head(x$isCofactor,100)
  auc <- sum(cumsum(tmp)/seq_along(tmp))
  if(is.null(x$padj)) x$padj <- x$adj.P.Val
  xsig <- x[which(x$padj<0.05),,drop=FALSE]
  data.frame(rank=w, trueQ=trueQ,
             FDR=sum(!xsig$isCofactor)/nrow(xsig),
             AUC=auc, relAUC=auc/optimalAUC)
}

# fetch a list of interactors for all requested TFs
getAllInteractors <- function(datasets, extra=c("CEBPB", "MAZ", "ZNF143", "NR3C1")){
  names(ds2) <- ds2 <- extra
  ds2 <- lapply(ds2, FUN=function(x) list(truth=x, species="h"))
  ints <- lapply(c(datasets,ds2), FUN=function(x){
    tf <- x$truth
    sp <- "9606"
    if(x$species=="m"){
      tf <- tools::toTitleCase(tolower(tf))
      sp <- "10090"
    }
    unique(toupper(unlist(lapply(tf, species=sp, FUN=getInteractors))))
  })
  return(ints)
}



transferDiffTFres <- function(infolder=".", outfolder, noPerm=FALSE){
  lf <- list.files(infolder, pattern="\\.summary\\.tsv\\.gz", full=TRUE, recursive=TRUE)
  if(noPerm){
    lf <- lf[grepl("noPerm",lf)]
  }else{
    lf <- lf[!grepl("noPerm",lf)]
  }
  m <- sapply(strsplit(dirname(lf),"/"),FUN=identity)
  names(lf) <- m[which.max(apply(m,1,FUN=function(x) length(unique(x)))),]
  
  if(file.exists(rtf <- paste0(infolder,"/diffTF_runtimes"))){
    tt <- read.delim(rtf, sep=" ", row.names=1, header = FALSE)
    colnames(tt)[1] <- "time"
    tt$dataset <- gsub("/.+$","",row.names(tt))
    tt$method <- sapply(strsplit(row.names(tt),"/"), FUN=function(x) x[2])
    tt <- tt[tt$time>0,]
    if(noPerm){
      tt <- tt[grepl("noPerm",tt$method),]
    }else{
      tt <- tt[!grepl("noPerm",tt$method),]
    }
    row.names(tt) <- tt$dataset
  }else{
    tt <- NULL
  }
  
  for(x in names(lf)){
    a <- read.delim(lf[[x]])
    b <- data.frame(row.names=a$TF, wmdiff=a$weighted_meanDifference, p=a$pvalue, padj=a$pvalueAdj)
    b <- b[order(b$p, -abs(b$wmdiff)),]
    b$rank <- seq_len(nrow(b))
    saveRDS(b, paste0(outfolder,"/",x,"/runATAC_results/with_pvalues/diffTF",
                      ifelse(noPerm,"_noPerm",""),".rds"))
    if(!is.null(tt) && x %in% row.names(tt)){
      rt <- c(user=NA_integer_, system=NA_integer_, elapsed=tt[x,1])
      rt <- list(res=NA, runtime=rt, runtime2=rt)
      saveRDS(rt, paste0(outfolder,"/",x,"/runATAC_results/raw/diffTF",
                         ifelse(noPerm,"_noPerm",""),"_raw.rds"))
    }
  }
}