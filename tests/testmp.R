library(nettools)

library(parallel)
library(minerva)

source("../R/adj_function.R")
source("../R/mat2adj.R")

## Create toy dataset
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)
cc <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a, method="MINE", n.cores=1)
badj <- mat2adj(b, method="MINE", n.cores=1)
cadj <- mat2adj(cc, method="MINE", n.cores=1)

myarr <- array(data=NA, c(100,100,3))
myarr[,,1] <- aadj
myarr[,,2] <- badj
myarr[,,3] <- cadj

iadj <- matrix(c(0,1,0,1,0,1,0,1,0), ncol=3, byrow=TRUE)

net1 <- list(I=iadj, L=myarr)
net2 <- net1

S <- diag(c(2,3,3,2))
W <- matrix(c(0,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0), ncol=4)

LI <- S - W

dd <- mpnetdist(net1, net2, d="HIM",n.cores=1, components=FALSE)

ll <- IntraAdj(myarr)
ll2 <- IntraAdj(myarr)



## Ll <- intraLaplacian(myarr)
## Li <- directProd(LI, n=7)

## mylap <- function(myarray){
##   nn <- dim(myarray)
##   mymatrix <- matrix(0,nrow=nn[1]*nn[3], ncol=nn[1]*nn[3])
  
##   for (i in 1: dim(myarr)[3]){
##     tmp <- Lap(myarray[,,i])
##     idx <- ((i-1) * nn[1]) + 1
##     mymatrix[seq(idx,nn[1]*i), seq(idx,nn[1]*i)] <- tmp
##   }
##   return(mymatrix)
## }



## library(microbenchmark)

## microbenchmark(
##   mylap(myarr), intraLaplacian(myarr)
##   )
 
