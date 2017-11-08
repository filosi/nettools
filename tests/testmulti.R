library(nettools)

## Create toy dataset
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)
cc <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a, infer.method="MINE", n.cores=1)
badj <- mat2adj(b, infer.method="MINE", n.cores=1)
cadj <- mat2adj(cc, infer.method="MINE", n.cores=1)

myarr <- array(data=NA, c(100,100,3))
myarr[,,1] <- aadj
myarr[,,2] <- badj
myarr[,,3] <- cadj

iadj <- matrix(c(0,1,0,1,0,1,0,1,0), ncol=3, byrow=TRUE)
net1 <- list(I=iadj, L=myarr)

## Second network
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)
cc <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a, infer.method="MINE", n.cores=1)
badj <- mat2adj(b, infer.method="MINE", n.cores=1)
cadj <- mat2adj(cc, infer.method="MINE", n.cores=1)

myarr <- array(data=NA, c(100,100,3))
myarr[,,1] <- aadj
myarr[,,2] <- badj
myarr[,,3] <- cadj

net2 <- list(I=iadj, L=myarr)

S <- diag(c(2,3,3,2))
W <- matrix(c(0,1,1,0,1,0,1,1,1,1,0,1,0,1,1,0), ncol=4)

LI <- S - W

dd <- mpnetdist(net1, net2, d="HIM",n.cores=1, components=TRUE)
