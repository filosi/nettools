## Normal correlation (default pearson)
Adjcor <- function(x,method='pearson',...){
    Adj <- abs(cor(x,method=method))
    diag(Adj) <- 0
    return(Adj)
}

## WGCNA method with cor^P
AdjWGCNA <- function(x,P,...){
    ## Default P = 6
    if (P!=6)
        warning(paste("WGCNA computed with P =",P,"!!!"))
    ## Adj <- WGCNA::unsignedAdjacency(x,
    ##                                 power = P,
    ##                                 corFnc = 'cor')
    Adj <- unsignedAdjacency(x,
                             power = P,
                             corFnc = 'cor')
    diag(Adj) <- 0
    return(Adj)
}

## Function for computing the FDR value on the null hypothesis
fdrrun <- function(Data,idx,FUN='cor',...){
  FUN <- match.fun(FUN)
  measure <- 1
  if (!is.na(match('measure',names(list(...)))))
      measure <- list(...)[['measure']]
  cormat <- matrix(Inf,ncol=dim(Data)[2],nrow=dim(Data)[2])
  if(dim(idx)[1]>0){
      for (i in seq(dim(idx)[1])){
          idi <- idx[i,]
          if (idi[1]>idi[2])
              cormat[idi[1],idi[2]] <- cormat[idi[2],idi[1]] <- abs(FUN(sample(Data[,idi[2]]),Data[,idi[1]])[[measure]])
      }

  }
  invisible(cormat)
}

## WGCNA FDR
AdjWGCNAFDR <- function(x,FDR,P,...){
    Adj <- AdjWGCNA(x,P=P)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
    for (i in seq(1/FDR)){
        cormat <- fdrrun(x,idx,FUN=cor)
        idx <- which(Adj>cormat,arr.ind=TRUE)
    }
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
    ## Adj <- WGCNA::bicor(x)
    Adj <- abs(bicor(x))
    diag(Adj) <- 0
    return(Adj)
}

## Bicor FDR
AdjbicorFDR <- function(x,FDR,P,...){
    Adj <- Adjbicor(x)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
    for (i in seq(1/FDR)){
        cormat <- fdrrun(x,idx,FUN='bicor')
        idx <- which(Adj>cormat,arr.ind=TRUE)
    }
    adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                       dimnames=list(rownames(Adj),colnames(Adj)))
    if(dim(idx)[1]>0){
        for (i in seq(dim(idx)[1])){
            adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
        }
    }
    return(adjfinal)
}

## TOM
AdjTOM <- function(x,P,...){
  if (P!=6)
    warning(paste("WGCNA computed with P =",P,"!!!"))
  ## Adj <- WGCNA::TOMsimilarityFromExpr(datExpr=x,
  ##                                     power=P,
  ##                                     corType="pearson",
  ##                                     TOMType="unsigned",
  ##                                     verbose=0)
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
    ##anet <- minet::minet(x,method="aracne",estimator="spearman")
    Adj <- minet(x,method="aracne",estimator="spearman")
    diag(Adj) <- 0
    return(Adj)
}

## CLR normalized
AdjCLR <- function(x,...){
    ## cnet <- minet::minet(x,method="clr",estimator="spearman")
    Adj <- minet(x,method="clr",estimator="spearman")
    diag(Adj) <- 0
    return(Adj)
}

## MINE
AdjMINE <- function(x,measure,alpha,C,...){
    Adj <- mine(x,alpha=alpha, C=C)[[measure]]
    if (!is.null(Adj)){
      diag(Adj) <- 0
      return(Adj)
    } else {
      stop("Invalid measure for method MINE")
    }
    
}

## MINEFDR
AdjMINEFDR <- function(x,measure,alpha,C,FDR,...){
    Adj <- AdjMINE(x,measure,alpha,C,...)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
    for (i in seq(1/FDR)){
        cormat <- fdrrun(x,idx,FUN='mine',measure=measure)
        idx <- which(Adj>cormat,arr.ind=TRUE)
    }
    adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                       dimnames=list(rownames(Adj),colnames(Adj)))
    if(dim(idx)[1]>0){
        for (i in seq(dim(idx)[1])){
            adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
        }
    }
    return(adjfinal)
}


## DTWMIC
AdjDTWMIC <- function(x,DP,...){
    Adj <- matrix(0,nrow=ncol(x),ncol=ncol(x))
    for(i in 1:(ncol(x)-1)){
        for(j in (i+1):ncol(x)){
            d <- 1/(1+dtw(x[,i],x[,j],distance.only=TRUE)$normalizedDistance)
            m <- mine(x[,i],x[,j])$MIC
            Adj[i,j] <- Adj[j,i]<- sqrt(d**2/2+m**2/2)**DP
        }
    }
    rownames(Adj) <- colnames(Adj) <- colnames(x)
    return(Adj)
}

