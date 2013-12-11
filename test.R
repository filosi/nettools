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
##
netdist(aadj,badj,method="HIM", components=FALSE)


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


g1 <- barabasi.game(200,0.8)
g2 <- barabasi.game(200,0.2)

a1 <- get.adjacency(g1,type="both",sparse=TRUE)
a2 <- get.adjacency(g2,type="both",sparse=TRUE)

## HIM distance
netdist(g1,g2,method="HIM")


##-------------------------------------------------
## Computing stability indicators
##-------------------------------------------------

netSI(a)

##is the same as

netSI(a,indicator="all", dist='HIM', adj.method='cor', 
      method="montecarlo", k=3, h=20)

netSI(a, components=TRUE)
