## Normal correlation (default pearson)
Adjcor <- function(x,method='pearson',...){

  ## Check for use argument
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }
  Adj <- abs(cor(x,method=method, use=use))
  diag(Adj) <- 0

  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
      Adj[is.na(Adj)] <- 0
  
  return(Adj)
}

## WGCNA method with cor^P
AdjWGCNA <- function(x,P,...){
  ## Default P = 6 send a warning if P!=6
  if (P!=6)
    warning(paste("WGCNA computed with P =",P,"!!!"))

  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## Compute the adjacency matrix (no auto-loop)
  Adj <- unsignedAdjacency(x,
                           power = P,
                           corFnc = 'cor',
                           corOptions=paste("use=\'", use,"\'", sep=""))
  diag(Adj) <- 0
  
  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
    Adj[is.na(Adj)] <- 0

  ## return Adj
  return(Adj)
}


## WGCNA FDR
AdjWGCNAFDR <- function(x,FDR,P,...){
  Adj <- AdjWGCNA(x,P=P, ...)
  idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
  
  ## check for multiple cores
  if (!is.na(match("n.cores", names(list(...))))){
    n.cores <- list(...)[["n.cores"]]
    if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
      if (n.cores > detectCores())
        n.cores <- detectCores() - 1
      cl <- makeCluster(n.cores)
    } else {
      cl <- NULL
    }
  } else {
    cl <- NULL
  }

  ## Start the FDR computation using fdrrun function
  for (i in seq(1/FDR)){
    if (nrow(idx) > 2){
      cormat <- fdrrun(x,idx,FUN=cor,cl,...)
      idx <- which(Adj>cormat,arr.ind=TRUE)
    }
  }

  ## create the adjacency matrix after the fdr filter
  adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                     dimnames=list(rownames(Adj),colnames(Adj)))
  if(dim(idx)[1]>0){
    for (i in seq(dim(idx)[1])){
      adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
    }
  }
  return(adjfinal)
}

## Bicor
Adjbicor <- function(x,...){
  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## compute the matrix
  Adj <- abs(bicor(x, use=use))
  diag(Adj) <- 0

  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
    Adj[is.na(Adj)] <- 0

  ## Return the Adj matrix
  return(Adj)
}

## Bicor FDR
AdjbicorFDR <- function(x,FDR,P,...){
    Adj <- Adjbicor(x, ...)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))

    ## Check for multicore
    if (!is.na(match("n.cores", names(list(...))))){
      n.cores <- list(...)[["n.cores"]]
      if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
        if (n.cores > detectCores())
          n.cores <- detectCores() - 1
        cl <- makeCluster(n.cores)
      } else {
        cl <- NULL
      }
    } else {
      cl <- NULL
    }

    ## Start the FDR computation using fdrrun function
    for (i in seq(1/FDR)){
      if (nrow(idx) > 2){
        cormat <- fdrrun(x,idx,FUN='bicor', cl, ...)
        idx <- which(Adj>cormat,arr.ind=TRUE)
      }
    }

    ## create the adjacency matrix after the fdr filter
    adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                       dimnames=list(rownames(Adj),colnames(Adj)))
    if(dim(idx)[1]>0){
        for (i in seq(dim(idx)[1])){
            adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
        }
    }

    ## Return the Adj matrix
    return(adjfinal)
}

## TOM
AdjTOM <- function(x,P,...){
  if (P!=6)
    warning(paste("WGCNA computed with P =",P,"!!!"))
  Adj <- TOMsimilarityFromExpr(datExpr=x,
                               power=P,
                               corType="pearson",
                               TOMType="unsigned",
                               verbose=0)
  diag(Adj) <- 0
  return(Adj)
}

## Aracne
AdjARACNE <- function(x,...){
  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## mutual information estimation does not accept use paramter!
  Adj <- minet(x,method="aracne",estimator="spearman")
  diag(Adj) <- 0

  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
    Adj[is.na(Adj)] <- 0
  
  return(Adj)
}

## CLR normalized
AdjCLR <- function(x,...){
  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## mutual information estimation does not accept use paramter!
  Adj <- minet(x,method="clr",estimator="spearman")
  diag(Adj) <- 0

  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
    Adj[is.na(Adj)] <- 0

  ## return Adj
  return(Adj)
}

## MINE
AdjMINE <- function(x,measure,alpha,C, var.thr=1e-10,...){
  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## Infer the adjacency matrix
  Adj <- mine(x,alpha=alpha, C=C,var.thr=var.thr, ...)[[measure]]
  
  ## Check whether the measure for MINE is available 
  if (is.null(Adj)){
    stop("Invalid measure for method MINE")
  } else {
    ## Set diagonal to 0, no auto-loop
    diag(Adj) <- 0

    ## use argument, remove NA from Adj and set to 0
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                               "everything"))
    if (na.method %in% c(2L, 3L))
      Adj[is.na(Adj)] <- 0
    return(Adj)
  }
}

## MINEFDR
AdjMINEFDR <- function(x,measure,alpha,C,FDR,...){
  ## Compute the step 0 for FDR filter
  Adj <- AdjMINE(x,measure,alpha,C,...)
  idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
  
  ## check for multicore
  if (!is.na(match("n.cores", names(list(...))))){
    n.cores <- list(...)[["n.cores"]]
    if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
      if (n.cores > detectCores())
        n.cores <- detectCores() - 1
      cl <- makeCluster(n.cores)
    } else {
      cl <- NULL
    }
  } else {
    cl <- NULL
  }
  
  for (i in seq(1/FDR)){
    if (nrow(idx) > 2){
      cormat <- fdrrun(x,idx,FUN='mine',measure=measure, cl,...)
      idx <- which(Adj>cormat,arr.ind=TRUE)
    }
  }
  adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                     dimnames=list(rownames(Adj),colnames(Adj)))
  if(nrow(idx)>2){
    for (i in seq(dim(idx)[1])){
      adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
    }
  }
  return(adjfinal)
}


## DTWMIC
AdjDTWMIC <- function(x,DP,...){
  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  
  ## Infer the adjacency matrix
  Adj <- matrix(0,nrow=ncol(x),ncol=ncol(x))
  for(i in 1:(ncol(x)-1)){
    for(j in (i+1):ncol(x)){
      d <- 1/(1+dtw(x[,i],x[,j],distance.only=TRUE)$normalizedDistance)
      m <- mine(x[,i],x[,j])$MIC
      Adj[i,j] <- Adj[j,i]<- sqrt(d**2/2+m**2/2)**DP
    }
  }
  rownames(Adj) <- colnames(Adj) <- colnames(x)

  ## use argument, remove NA from Adj and set to 0
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (na.method %in% c(2L, 3L))
    Adj[is.na(Adj)] <- 0

  ## Return Adj
  return(Adj)
}

## Function for computing the FDR value on the null hypothesis
fdrrun <- function(x,idx,FUN='cor', cl=NULL, ...){
  FUN <- match.fun(FUN)

  ## Set the measure for MINE statistics
  measure <- 1
  if (!is.na(match("measure",names(list(...))))){
    measure <- list(...)[["measure"]]
    if (length(measure)==0)
      measure <- 1
  }

  ## Check for use argument in ...
  param <- list(...)
  nn <- names(param)
  if (!is.na(match("use", nn))){
    use <- param$use
  } else {
    use <- "all.obs"
  }

  ## Allocate space for the matrix
  cormat <- matrix(Inf,ncol=dim(x)[2],nrow=dim(x)[2])
  myfun <- function(x, y, idx, tmpfun, measure, use){
    idi <- idx[x,]
    if (idi[1]>idi[2]){
      return(abs(tmpfun(sample(y[,idi[2]]),y[,idi[1]], use=use)[[measure]]))
    } else {
      return(NA)
    }
  }
  if(nrow(idx)>2){
    if (is.null(cl)){
      aa <- sapply(seq(dim(idx)[1]),myfun, y=x, idx=idx, tmpfun=FUN, measure=measure, use=use)
    } else {
      ## NB need to load package on all the cores?!??!?!
      aa <- parSapply(cl, seq(dim(idx)[1]),myfun, y=x, idx=idx, tmpfun=FUN, measure=measure, use=use)
    }
    idx <- cbind(idx,aa)
    tmpid <- subset(idx,!is.na(idx[,3]))
    cormat[tmpid[,c(1,2)]] <- cormat[tmpid[,c(2,1)]] <- tmpid[,3]
  } else {
    stop("Not conformable array")
  }
  invisible(cormat)  
}

## Function for check the variance by features
checkvar <- function(x, tol=1e-5, ...){
  
  ## Compute the variance by columns
  feat.var <- apply(x,2,var)
  
  ## Get the indexes of the feature with low variance
  idx <- which(feat.var<tol)
  if (length(idx) == 0L){
    return (NULL)
  } else {
    return (idx)
  }
}
