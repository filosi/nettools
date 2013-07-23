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

    ## Create the list of adjacency matrix
    myadj <- list(method=METHODS[method],G1=g1,G2=g2,N=ncol(g1))

    ## Check method for distances
    if (myadj$method=="HIM"){
      mylap <- list(L1=Lap(g1),L2=Lap(g2),N=ncol(g1))
      dd <- him(list(myadj,mylap))
    }
    if (myadj$method=="ipsen"){
      mylap <- list(L1=Lap(g1),L2=Lap(g2),N=ncol(g1))
      dd <- ipsen(mylap,gamma=gamma)
    }
    if (myadj$method=="hamming"){
      dd <- hamming(myadj)
    }
    return(dd)
  }

## Create the adjacency matrix structure
g2adj <- function(x,...) UseMethod("g2adj")

g2adj.igraph <- function(x,...,type="both"){
  if (!is.null(get.edge.attribute(x,"weight"))){
    WW <- "weight"
  } else {
    warning("No weight attribute to the graph object\nCompute binary adjacency matrix")
    WW <- NULL
  }
  Adj <- get.adjacency(x,type=type,attr=WW,sparse=TRUE)
  diag(Adj) <- 0
  return(Adj)
}
setMethod("g2adj","igraph",g2adj.igraph)

## Generical Laplacian
Lap <- function(x,...) UseMethod("Lap")
Lap.default <- function(x,...){
  D <- apply(x,2,sum)
  L <- -x
  diag(L) <- D
  return(L)
  #return((D * diag(dim(x)[1])) - x)
}
setMethod("Lap","matrix",Lap.default)
setMethod("Lap","Matrix",Lap.default)

## Ipsen distance
ipsen <- function(object,...) UseMethod("ipsen")
ipsen.list <- function(object,...,gamma=NULL){
  if (is.null(gamma))
    optgamma <- optimal_gamma(object$N)
  else
    optgamma <- gamma
  if(object$N>1000){
    cl <- makeCluster(getOption("cl.cores",2))
    clusterEvalQ(cl,require("nettools",quietly=TRUE))
    ll <- clusterApply(cl,1:2,function(x,mygamma=optgamma,mylist=object,...){
      myomega <- sqrt(abs(round(spec(object[[x]]),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    })
    stopCluster(cl)
  } else {
    ll <- lapply(list(object[[1]],object[[2]]),function(x,mygamma=optgamma,...){
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
hamming <- function(object,...) UseMethod("hamming")
hamming.list <- function(object,...){
  ## for weighted networks, weights must be in [0,1]
  n <- object$N
  return(sum(abs(object$G1-object$G2))/(n*(n-1))  )
}
setMethod("hamming","list",hamming.list)


## Him distance
him <- function(object,...) UseMethod("him")
him.list <- function(object,...){
  ipd <- ipsen(object[[2]],gamma=NA)
  had <- hamming(object[[1]])
  gloc <- sqrt(had**2/2+ipd**2/2)
  return(c("H"=had,"IM"=ipd,"HIM"=gloc))
}
setMethod("him","list",him.list)

