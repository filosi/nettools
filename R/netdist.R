netdist <- function(x, ...) UseMethod("netdist")

## Method netdist for matrix
netdist.matrix <- function(x, h, d="HIM", ga=NULL, components=TRUE, ...){
  
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
  
  #add a check so that an unexisting parameter cannot be passed
  id.Call <- match( names(Call),c("x", "h", "d", "ga","components","n.nodes"), nomatch=0)
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
  
  ## n.nodes (1000 by default) number of nodes where to start parallel computation for
  ## ipsen and him distance
  ## if number of nodes < n.nodes it won't start,
  ## otherwise a 2 process will be launched
  if (is.null(Call$n.nodes))
    n.nodes <- 1000 else n.nodes <- eval(Call$n.nodes)
  
  if(is.na(d))
    stop("invalid distance")
  if(d == -1)
    stop("ambiguous distance")
  
  ##check on components argument
  if(is.null(Call$components)){
    if(d==1){
      comp <- TRUE
    }else{
      comp <- FALSE
    }
  }else{
    comp <- eval(Call$components)
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    }else{
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }
  
  ##check on ga passing through ipsen function
  if(is.null(Call$ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    ga <- eval(Call$ga)
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
    dd <- him(list(ADJ=myadj,LAP=mylap), n.nodes, ga=ga,  components=comp, ltag=FALSE, ...)
  }
  if (myadj$d=="IM"){
    mylap <- list(L=list(Lap(g1$adj), Lap(g2$adj)),N=g1$N, tag=g1$tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- ipsen(mylap,ga=ga, n.nodes, ...)
  }
  if (myadj$d=="H"){
    dd <- hamming(myadj)
  }
  
  return(dd)
}
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
  id.Call <- match( names(Call),c("x", "d", "ga","components","n.nodes"), nomatch=0)
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
  
  ## n.nodes (1000 by default) number of nodes where to start parallel computation for
  ## ipsen and him distance
  ## if number of nodes < n.nodes it won't start,
  ## otherwise a 2 process will be launched
  if (is.null(Call$n.nodes))
    n.nodes <- 1000 else n.nodes <- eval(Call$n.nodes)
  
  if(is.na(d))
    stop("invalid distance")
  if(d == -1)
    stop("ambiguous distance")
  
  ##check on components argument
  if(is.null(Call$components)){
    if(d==1){
      comp <- TRUE
    }else{
      comp <- FALSE
    }
  }else{
    comp <- eval(Call$components)
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    }else{
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }
  
  ##check on ga passing through ipsen function
  if(is.null(Call$ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    ga <- eval(Call$ga)
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
      stop("Not the same length!")
    }
    myadj <- list(d=DISTANCE[d],G=adjlist,N=N,tag=tag)
  }
  
  ## Check for dimension and directionality of g and h
  ## Return two different error messages
  ## if ((g1$N == g2$N)){
  ##   if((g1$tag == g2$tag)){
  ##     ## Create the list of adjacency matrix
  ##     myadj <- list(d=DISTANCE[d],G=list(g1$adj,g2$adj),N=g1$N,tag=g1$tag)
  ##     ## myadj <- list(d=DISTANCE[d],G1=g1$adj,G2=g2$adj,N=g1$N,tag=g1$tag)
  ##   } else {
  ##     stop("Not conformable graph g and h: one is directed while the other is undirected", call.=FALSE)
  ##   }
  ## } else {
  ##   stop("Not conformable graph g and h: they have different dimensions", call.=FALSE)
  ## }
  
  ## Check the distance chosen
  if (myadj$d=="HIM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- him(list(ADJ=myadj,LAP=mylap), n.nodes, ga=ga,  components=comp, ltag=TRUE, ...)
  }
  if (myadj$d=="IM"){
    mylap <- list(L=laplist,N=N, tag=tag)
    ## mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
    dd <- ipsen(mylap,ga=ga, n.nodes, ...)
  }
  if (myadj$d=="H"){
    dd <- hamming(myadj)
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


## Generical Laplacian
##----------------------------------------
Lap <- function(x,...){
  D <- apply(x,2,sum)
  L <- -x
  diag(L) <- D
  return(L)
}

## Ipsen distance
##----------------------------------------
## ipsen <- function(object,...) UseMethod("ipsen")
ipsen <- function(object, ga=NULL, n.nodes=1000, ...){
  if (is.null(ga)){
    if (object$tag == "undir"){
      optgamma <- optimal_gamma(object$N)
    } else {
      optgamma <- optimal_gamma_dir(object$N)
    }
  } else {
    optgamma <- ga
  }

  laplist <- object$L

  ## Check for parallelization
  n.cores <- NULL
  if (!is.na(match("n.cores",names(list(...)))))
    n.cores <- list(...)[["n.cores"]]

  ## Should I use multiple cores or not?
  if(object$N>n.nodes && detectCores() >= 2){
    if (is.null(n.cores) || n.cores >= detectCores()){
      if (length(laplist) < detectCores()){
        n.cores <- length(laplist)
      } else {
        n.cores <- detectCores() - 1
      }
    }
    cl <- makeCluster(n.cores)
    ## Eval needed function on nodes
    clusterEvalQ(cl,{K <- nettools:::K
                     rho <- nettools:::rho
                     lorentz <- nettools:::lorentz
                   })

    cat("Start computing eigenvalues with ",n.cores, "cores\n")
    ## Actual computation of eigen-values/vectors
    ll <- clusterApply(cl,laplist,function(x,mygamma=optgamma,...){
      myomega <- sqrt(abs(round(spec(x),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    })
    stopCluster(cl)
  } else {
    ## Computation on 1 CPU
    ll <- lapply(1:length(laplist),function(x,mygamma,laplist, ...){
      print(paste("Computing eigen for ",x))
      aa <- laplist[[x]]
      myomega <- sqrt(abs(round(spec(aa),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    }, mygamma=optgamma, laplist=laplist, ...)
  }
  mydistfun <- function(a,b, optgamma){
    integrand <- function(omega, mygamma, given_omega_G, given_omega_H){
      (rho(omega, optgamma,a)-rho(omega,optgamma,b))**2
    }
    tmp <- sqrt(integrate(integrand,lower=0,upper=Inf,mygamma=optgamma,given_omega_G=a[[1]],given_omega_H=b[[1]], stop.on.error=FALSE,rel.tol=.Machine$double.eps,subdivisions=1e4)$value)
    return(tmp)
  }
  cat("Start computing mutual distances")
  if (length(laplist) == 2){
    ## Compute distance between 2 adjacency matrices
    dist <- mydistfun(ll[[1]], ll[[2]], optgamma=optgamma)
    names(dist) <- "IM"
  } else {
    ## Compute mutual distances between all the matrices in the list
    idx <- combn(length(ll),2)
    tmpdist <- sapply(1:dim(idx)[2], function(x,ll,optgamma, idx){
      print(paste("Distance",idx[1,x],"vs", idx[2,x]))
      mydistfun(ll[[idx[1,x]]], ll[[idx[2,x]]], optgamma)
    }, ll=ll, optgamma=optgamma, idx=idx)
    dist <- matrix(NA,ncol=length(ll), nrow=length(ll))
    dist[t(idx)] <- dist[t(idx)[,c(2,1)]] <- tmpdist
    diag(dist) <- 0
  }
  return(dist)
}

## Hamming distance
##----------------------------------------
hamming <- function(object,...){
  
  adjlist <- object$G
  
  ## for weighted networks, weights must be in [0,1]
  if (object$tag == "undir"){
    dist <- ham.undir(adjlist, object$N, ...)
  } else{
    dist <- ham.dir(adjlist, object$N, ...)
  }
  return(dist)
}
## setMethod("hamming","list",hamming.list)

## Him distance
##----------------------------------------
him <- function(object,ga=NULL, components=TRUE, n.nodes=1000, ltag=FALSE, ...){
  ipd <- ipsen(object$LAP, ga, n.nodes, ...)
  had <- hamming(object$ADJ)
  gloc <- sqrt(had**2/2+ipd**2/2)
  if(components==TRUE){
    if (ltag){
      dist <- list(H=had,I=ipd,HIM=gloc)
    } else {
      dist <- c(had, ipd, gloc)
      names(dist) <- c("H","IM","HIM")
    }
  } else {
    dist <- gloc
    names(dist) <- "HIM"
  }
  return(dist)
}
