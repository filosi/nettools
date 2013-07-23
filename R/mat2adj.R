## Compute the adjacency matrix with multiple methods on a N x P matrix (x)
## x =

mat2adj <- function(x,...) UseMethod("mat2adj")

mat2adj.data.frame <- function(x,...){
  x <- apply(x,2,as.numeric)
  x <- as.matrix(x)
  UseMethod("mat2adj", x)
}
setMethod("mat2adj","data.frame",mat2adj.data.frame)

mat2adj.matrix <- function(x,method='cor',FDR=1e-3,P=6,measure='MIC',alpha=0.6,C=15,DP=1,...){

    METHODS <- c('cor','WGCNA','WGCNAFDR','bicor',
                 'bicorFDR','TOM','ARACNE','CLR',
                 'MINE','MINEFDR','DTWMIC')
    MEASURE <- c('MIC','MAS','MCN','MEV','MIC-R2')
    method <- pmatch(method, METHODS)
    if(is.na(method))
        stop("invalid distance method")
    if(method == -1)
        stop("ambiguous distance method")
    if(method == 3L)
        if (P > 1){
            P <- 1
            warning("Using WGCNAFDR method with P > 1, not yet implemented,\n P will be ignored")
        }
    ## If mine suite is require, choose the metric
    if(method >= 9L && method <= 10L)
        measure <- pmatch(measure,MEASURE)
    myfun <- paste('Adj',METHODS[method],sep='')
    tmp <- do.call(myfun,list(x=x,FDR=FDR,P=P,measure=MEASURE[measure],alpha=alpha,C=C,DP=DP,...))
    return(tmp)
}
setMethod("mat2adj","matrix",mat2adj.matrix)
