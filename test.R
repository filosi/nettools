##--------------------------------------------------
## Testing nettools usage
##--------------------------------------------------

library(nettools)

## Create toy dataset
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a,method="MINE",measure="MICR2")
badj <- mat2adj(b)

## Compute the him distance between the 2 matrix

## Him distance
netdist(aadj,badj,method="HIM")
##


## Hamming distance
netdist(aadj,badj,method="ham")
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
require(nettools)

aadj[aadj<0.1] <- 0
aadj[aadj!=0] <- 1
badj[badj<0.1] <- 0
badj[badj!=0] <- 1

g1 <- graph.adjacency(aadj)#,weighted=TRUE)
g2 <- graph.adjacency(badj)#,weighted=TRUE)


g1 <- barabasi.game(100,0.8)
g2 <- barabasi.game(100,0.2)

a1 <- get.adjacency(g1,type="both",sparse=TRUE)
a2 <- get.adjacency(g2,type="both",sparse=TRUE)

## HIM distance
netdist(g1,g2,method="HIM")


