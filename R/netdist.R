netdist <- function(x, ...) UseMethod("netdist")

## Method netdist for matrix
netdist.matrix <- function(x, h=NULL, d="HIM", ga=NULL, components=TRUE, ...){

  if (is.null(h))
    stop("Need to provide a second matrix to compute the distance")
  
  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L
  
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)
  
  Call <- match.call()
  
  ##add a check so that an unexisting parameter cannot be passed
  id.Call <- match( names(Call),c("x", "h", "d", "ga","components","n.cores","verbose", "rho"), nomatch=0)
  if(sum(id.Call[-1]==0)==1){
    warning("The parameter '",names(Call)[which(id.Call==0)[2]],"' will be ignored",call.=FALSE)
  }
  if(sum(id.Call[-1]==0)>1){
    msg <- "The following parameters will be ignored:\n"
    for(i in which(id.Call==0)[-1]){
      msg <- paste(msg,"'",names(Call)[i],"'\n",sep="")
    }
    warning(msg,call.=FALSE)
  }
  
  if(is.na(d))
    stop("invalid distance")
  if(d == -1)
    stop("ambiguous distance")
  
  ##check on components argument
  if(is.null(components)){
    if(d==1){
      comp <- TRUE
    }else{
      comp <- FALSE
    }
  }else{
    comp <- components
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    }else{
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }
  
  ##check on ga passing through ipsen function
  if(is.null(ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    if(!is.numeric(ga) && !is.null(ga))
      stop("ga must be numeric",call.=FALSE)
  }
  
  ## Create the class
  g1 <- g2adj(x)
  g2 <- g2adj(h)
  
  ## Check for dimension and directionality of g and h
  ## Return two different error messages
  if ((g1$N == g2$N)){
    if((g1$tag == g2$tag)){
      ## Create the list of adjacency matrix
      myadj <- list(d=DISTANCE[d],G=list(g1$adj,g2$adj),N=g1$N,tag=g1$tag)
      ## myadj <- list(d=DISTANCE[d],G1=g1$adj,G2=g2$adj,N=g1$N,tag=g1$tag)
    } else {
      stop("Not conformable graph g and h: one is directed while the other is undirected", call.=FALSE)
    }
  } else {
    stop("Not conformable graph g and h: they have different dimensions", call.=FALSE)
  }
  
  ## Check the distance chosen
  if (myadj$d=="HIM"){
    mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- him(list(ADJ=myadj,LAP=mylap),  ga=ga,  components=comp, ltag=FALSE, ...)
  }
  if (myadj$d=="IM"){
    mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- ipsen(mylap,ga=ga, ...)
  }
  if (myadj$d=="H"){
    dd <- hamming(myadj)
  }
  
  return(dd)
}
netdist.Matrix <- netdist.matrix

setMethod("netdist","matrix", netdist.matrix)
setMethod("netdist","Matrix", netdist.matrix)
setMethod("netdist","data.frame", netdist.matrix)

## Method netdist for list of adjacency matrices
##--------------------------------------------------
## Implementation for kernel distance
## Output a distance matrix
netdist.list <- function(x, d="HIM", ga=NULL, components=TRUE, ...){
  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L
  
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)
  
  Call <- match.call()
  
  ## add a check so that an unexisting parameter cannot be passed
  id.Call <- match( names(Call),c("x", "d", "ga","components", "n.cores", "verbose", "rho"), nomatch=0)
  if(sum(id.Call[-1]==0)==1){
    warning("The parameter '",names(Call)[which(id.Call==0)[2]],"' will be ignored",call.=FALSE)
  }
  if(sum(id.Call[-1]==0)>1){
    msg <- "The following parameters will be ignored:\n"
    for(i in which(id.Call==0)[-1]){
      msg <- paste(msg,"'",names(Call)[i],"'\n",sep="")
    }
    warning(msg,call.=FALSE)
  }
  
  if(is.na(d))
    stop("invalid distance")
  if(d == -1)
    stop("ambiguous distance")
  
  ##check on components argument
  if(is.null(components)){
    if(d==1){
      comp <- TRUE
    }else{
      comp <- FALSE
    }
  }else{
    comp <- components
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    }else{
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }
  
  ##check on ga passing through ipsen function
  if(is.null(ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    if(!is.numeric(ga) && !is.null(ga))
      stop("ga must be numeric",call.=FALSE)
  }
  
  ## Write the distance computation here!
  ## Create the class
  tmp <- lapply(x,g2adj)
  if (d == 1L | d == 2L)
    laplist <- list()
  
  adjlist <- list()
  
  N <- tmp[[1]]$N
  tag <- tmp[[1]]$tag
  for (i in 1:length(tmp)){
    g <- tmp[[i]]
    if (g$N == N && g$tag == tag){
        if (d == 1L | d == 2L)
          laplist[[i]] <- Lap(g$adj)
        adjlist[[i]] <- g$adj
    } else {
      stop(paste("Element",i,"not of the same length!"), call.=FALSE)
    }
    myadj <- list(d=DISTANCE[d],G=adjlist,N=N,tag=tag)
  }
  
  ## Check the distance chosen
  if (myadj$d=="HIM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- him(list(ADJ=myadj,LAP=mylap), ga=ga,  components=comp, ltag=TRUE, ...)
  }
  if (myadj$d=="IM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- ipsen(mylap,ga=ga, ...)
    ## if (!is.null(names(x)))
    ##   colnames(dd) <- rownames(dd) <- names(x)
  }
  if (myadj$d=="H"){
    dd <- hamming(myadj)
    ## if (!is.null(names(x)))
    ##   colnames(dd) <- rownames(dd) <- names(x)
  }
  
  ## Give names to the matrices
  if (!is.null(names(x))){
    if (is.list(dd)){
      dd <- lapply(dd,function(y,x){colnames(y) <- rownames(y) <- names(x)
                                    return(y)}, x=x)
    } else {
      colnames(dd) <- rownames(dd) <- names(x)
    }
  }
  
  return(dd)
}
setMethod("netdist", "list", netdist.list)


## Prepare the matrix for computing distance if the graph is directed
transfmat <- function(x){
  ## Check if the x matrix is symmetric (undirected graph)
  ## Otherwise it returns a list with a matrix like:
  ## |zeros    t(A)|
  ## |  A     zeros|
  ##
  Adj <- x
  n <- ncol(Adj)
  tag <- "undir"
  
  ## If the graph is directed create a new matrix 
  if (!isSymmetric(x,check.attributes=FALSE, check.names=FALSE)){
    zero <- matrix(0, nrow=n, ncol=n)
    tmp <- Matrix::rBind(Matrix::cBind(zero,t(Adj)),Matrix::cBind(Adj,zero))
    Adj <- tmp
    tag <- "dir"
  }

  ## Normalize the weights
  if (any(Adj > 1) || any(Adj < 0)){
    warning("Edge weight should be >= 0 and <= 1, scaling has been automatically applied!", call.=FALSE)
    Adj <- (Adj - min(Adj)) / (max(Adj) - min(Adj))
  }
  return(list(adj = Adj, tag = tag, N=n))
}

## Create the adjacency matrix structure
##----------------------------------------
g2adj <- function(x,...) UseMethod("g2adj")

g2adj.igraph <- function(x,...,type="both"){
  if (!is.null(get.edge.attribute(x,"weight"))){
    WW <- "weight"
  } else {
    warning("No weight attribute to the graph object\nCompute binary adjacency matrix", call.=FALSE)
    WW <- NULL
  }
  Adj <- get.adjacency(x,type=type,attr=WW,sparse=TRUE)
  diag(Adj) <- 0
  ll <- transfmat(Adj)

  return(ll)
}
setMethod("g2adj","igraph",g2adj.igraph)

g2adj.matrix <- function(x,...){
  ll <- transfmat(x)
  
  return(ll)
}
setMethod("g2adj","matrix",g2adj.matrix)
setMethod("g2adj","Matrix",g2adj.matrix)

g2adj.data.frame <- function(x, ...){
  x <- apply(x,2,as.numeric)
  ll <- transfmat(x)
  return(ll)
}
setMethod("g2adj","data.frame",g2adj.data.frame)



## Generical Laplacian
##----------------------------------------
Lap <- function(x,...){
  D <- apply(x,2,sum)
  L <- -x
  diag(L) <- D
  return(L)
}


## Him distance
##----------------------------------------
him <- function(object,ga=NULL, components=TRUE, ltag=FALSE, rho=1, ...){
  ipd <- ipsen(object$LAP, ga, ...)
  had <- hamming(object$ADJ)
  gloc <- sqrt(1/(1+rho)) * sqrt(had**2+ rho*(ipd**2))
  if (length(gloc)==1)
    names(gloc) <- "HIM"
  if(components==TRUE){
    if (ltag){
      dist <- list(H=had,IM=ipd,HIM=gloc)
    } else {
      dist <- c(had, ipd, gloc)
      names(dist) <- c("H","IM","HIM")
    }
  } else {
    dist <- gloc
    if (!ltag)
      names(dist) <- "HIM"
  }
  return(dist)
}
