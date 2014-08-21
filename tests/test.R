##--------------------------------------------------
## Testing nettools usage
##--------------------------------------------------

library(nettools)

## Create toy dataset
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a, method="MINE")
badj <- mat2adj(b, method="MINE")

## Compute the him distance between the 2 matrix

## Him distance
netdist(aadj,badj,method="HIM",n.cores=1,components=FALSE)
## HIM distance using list of matrices as input
netdist(list(aadj,badj),method="HIM", components=FALSE, n.cores=1)


## Hamming distance
dd <- netdist(aadj,badj,method="ham")
## [1] 0.2134906
##

## Ipsen distance
netdist(aadj,badj,method="ips",gamma=0.09)
## [1] 0.1358396
##

##-------------------------------------------------
## Using igraph
##-------------------------------------------------
require(igraph)

g1 <- barabasi.game(200,0.8)
g2 <- barabasi.game(200,0.2)

a1 <- get.adjacency(g1,type="both",sparse=TRUE)
a2 <- get.adjacency(g2,type="both",sparse=TRUE)

## HIM distance
netdist(g1,g2,method="HIM")


##-------------------------------------------------
## Computing stability indicators
##-------------------------------------------------

ssind <- netSI(a)


