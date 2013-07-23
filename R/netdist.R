## netdist <- function(g, h, method="HIM", gamma=0.08){

##     METHODS <- c('HIM','ipsen','hamming')
##     method <- pmatch(method, METHODS)
    
##     if(is.na(method))
##       stop("invalid distance method")
##     if(method == -1)
##       stop("ambiguous distance method")
        
##     g1 <- g2adj(g)
##     g2 <- g2adj(h)
    
##     ## Create the class
##     myadj <- new("adj2",G1=g1,G2=g2,
##                  N1=ncol(g1),N2=ncol(g2),
##                  L1=Lap(g1),#matrix(0,ncol=ncol(g1),nrow=ncol(g1)),
##                  L2=Lap(g2),#matrix(0,ncol=ncol(g2),nrow=ncol(g2)),
##                  method=METHODS[method])
    
##     if (myadj@method=="HIM")
##         dd <- him(myadj)
##     if (myadj@method=="ipsen")
##         dd <- ipsen(myadj,gamma=gamma)
##     if (myadj@method=="hamming")
##         dd <- hamming(myadj)      
    
##     return(dd)
## }

## ## Create the class structure
## g2adj <- function(x,...) UseMethod("g2adj")
## g2adj.igraph <- function(x,...,type="both"){
##     Adj <- get.adjacency(x,type=type,attr="weight")
##     Adj <- as.matrix(Adj)
##     diag(Adj) <- 0
##     return(Adj)
## }
## setMethod("g2adj","igraph",g2adj.igraph)
## setMethod("g2adj","data.frame",g2adj.default)

## ## Generical Laplacian
## Lap <- function(x,...) UseMethod("Lap")
## Lap.default <- function(x,...){
##   D <- apply(x,2,sum)
##   return((D * diag(dim(x)[1])) - x)
## }
## setMethod("Lap","matrix",Lap.default)

## ## Ipsen distance
## ipsen <- function(object,...) UseMethod("ipsen")
## ipsen.adj2 <- function(object,...,gamma=0.08){
##   if (is.na(gamma))
##     optgamma <- optimal_gamma(object@N1)
##   else
##     optgamma <- gamma
##   if(object@N1>1000){
##     cl <- makeCluster(getOption("cl.cores",2))
##     clusterEvalQ(cl,require("nettools",quietly=TRUE))
##     ll <- clusterApply(cl,list(object@L1,object@L2),function(x,mygamma=optgamma,...){
##       myomega <- sqrt(abs(round(spec(x),5)))
##       myk <- K(mygamma,myomega)
##       return(list(myomega,myk))
##     })
##     stopCluster(cl)
##   } else {
##     ll <- lapply(list(object@L1,object@L2),function(x,mygamma=optgamma,...){
##                 myomega <- sqrt(abs(round(spec(x),5)))
##                 myk <- K(mygamma,myomega)
##                 return(list(myomega,myk))
##                 })
##   }
##   integrand <- function(omega, mygamma, given_omega_G, given_omega_H){
##     (rho(omega, optgamma,ll[[1]])-rho(omega,optgamma,ll[[2]]))**2
##   }
##   dist <- sqrt(integrate(integrand,lower=0,upper=Inf,mygamma=optgamma,given_omega_G=ll[[1]][[1]],given_omega_H=ll[[2]][[1]], stop.on.error=FALSE,rel.tol=.Machine$double.eps,subdivisions=1e4)$value)
##   return(dist)
## }
## setMethod("ipsen","adj2",ipsen.adj2)


## ## Hamming distance
## hamming <- function(object,...) UseMethod("hamming")
## hamming.adj2 <- function(object,...){
##   ## for weighted networks, weights must be in [0,1]
##   n <- object@N1
##   return(sum(abs(object@G1-object@G2))/(n*(n-1)))
## }
## setMethod("hamming","adj2",hamming.adj2)


## ## Him distance
## him <- function(object,...) UseMethod("him")
## him.adj2 <- function(object,...){
##   ipd <- ipsen(object,gamma=NA)
##   had <- hamming(object)
##   gloc <- sqrt(had**2/2+ipd**2/2)
##   return(c("H"=had,"IM"=ipd,"HIM"=gloc))
## }
## setMethod("him","adj2",him.adj2)

