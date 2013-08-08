netdist <- function(g, h, method="HIM", gamma=NULL){

    METHODS <- c('HIM','ipsen','hamming')
    method <- pmatch(method, METHODS)

    if(is.na(method))
      stop("invalid distance method")
    if(method == -1)
      stop("ambiguous distance method")

    ## Create the class
    g1 <- g2adj(g)
    g2 <- g2adj(h)
    
    ## Check for dimension and directionality of g and h
    if ((g1$N == g2$N) & (g1$tag == g2$tag)){
      ## Create the list of adjacency matrix
      myadj <- list(method=METHODS[method],G1=g1$adj,G2=g2$adj,N=g1$N,tag=g1$tag)
    } else {
      stop("Not conformable graph g and h", call.=FALSE)
    }

    ## Check method for distances
    if (myadj$method=="HIM"){
      mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
      dd <- him(list(ADJ=myadj,LAP=mylap))
    }
    if (myadj$method=="ipsen"){
      mylap <- list(L1=Lap(g1$adj),L2=Lap(g2$adj),N=g1$N, tag=g1$tag)
      dd <- ipsen(mylap,gamma=gamma)
    }
    if (myadj$method=="hamming"){
      dd <- hamming(myadj)
    }
    return(dd)
  }



## Check if a matrix is symmetric
## cksymm <- function(x,...){
##   if (all((x - t(x)) == 0))
##     return(TRUE)
##   else
##     return(FALSE)
## }


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
  
  ## If the graph is directed create a new matrix (function undir)
  if (!isSymmetric(x,check.attributes=FALSE, check.names=FALSE)){
    zero <- matrix(0, nrow=n, ncol=n)
    tmp <- Matrix::rBind(Matrix::cBind(zero,t(Adj)),Matrix::cBind(Adj,zero))
    Adj <- tmp
    tag <- "dir"
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


## Generical Laplacian
##----------------------------------------
Lap <- function(x,...) UseMethod("Lap")
## Lap.default <- function(x,...){
##   D <- apply(x,2,sum)
##   L <- -x
##   diag(L) <- D
##   return(L)
##   #return((D * diag(dim(x)[1])) - x)
## }
setMethod("Lap","matrix",Lap.default)
setMethod("Lap","Matrix",Lap.default)

## Ipsen distance
##----------------------------------------
ipsen <- function(object,...) UseMethod("ipsen")
ipsen.list <- function(object,...,gamma=NULL){
  if (is.null(gamma)){
    if (object$tag == "undir"){
      optgamma <- optimal_gamma(object$N)
    } else {
      optgamma <- optimal_gamma_dir(object$N)
    }
  } else {
    optgamma <- gamma
  }
  
  ## Check if network is directed or not
  if(object$N>1000 && detectCores() >= 2){
    cl <- makeCluster(getOption("cl.cores",2))
    clusterEvalQ(cl,{K <- nettools:::K
                     rho <- nettools:::rho
                     lorentz <- nettools:::lorentz
                   })
    ll <- clusterApply(cl,1:2,function(x,mygamma=optgamma,mylist=object,...){
      myomega <- sqrt(abs(round(spec(object[[x]]),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    })
    stopCluster(cl)
  } else {
    ll <- lapply(list(object$L1,object$L2),function(x,mygamma=optgamma,...){
                myomega <- sqrt(abs(round(spec(x),5)))
                myk <- K(mygamma,myomega)
                return(list(myomega,myk))
                })
  }
  integrand <- function(omega, mygamma, given_omega_G, given_omega_H){
    (rho(omega, optgamma,ll[[1]])-rho(omega,optgamma,ll[[2]]))**2
  }
  dist <- sqrt(integrate(integrand,lower=0,upper=Inf,mygamma=optgamma,given_omega_G=ll[[1]][[1]],given_omega_H=ll[[2]][[1]], stop.on.error=FALSE,rel.tol=.Machine$double.eps,subdivisions=1e4)$value)
  return(dist)
}
setMethod("ipsen","list",ipsen.list)

## Hamming distance
##----------------------------------------
hamming <- function(object,...) UseMethod("hamming")
hamming.list <- function(object,...){
  ## for weighted networks, weights must be in [0,1]
  if (object$tag == "undir"){
    return(ham.undir(object, ...))
  } else{
    return(ham.dir(object, ...))
  }
}
setMethod("hamming","list",hamming.list)

## Him distance
##----------------------------------------
him <- function(object,...) UseMethod("him")
him.list <- function(object,...){
  ipd <- ipsen(object$LAP,gamma=NULL)
  had <- hamming(object$ADJ)
  gloc <- sqrt(had**2/2+ipd**2/2)
  return(c("H"=had,"IM"=ipd,"HIM"=gloc))
}
setMethod("him","list",him.list)
