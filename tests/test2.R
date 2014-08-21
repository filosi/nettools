library(nettools)

a <- matrix(rnorm(2000),ncol=200)
b <- matrix(rnorm(2000),ncol=200)
d <- matrix(rnorm(2000),ncol=200)
e <- matrix(rnorm(2000),ncol=200)
f <- matrix(rnorm(2000),ncol=200)
g <- matrix(rnorm(2000),ncol=200)
h <- matrix(rnorm(2000),ncol=200)
i <- matrix(rnorm(2000),ncol=200)
mydlist <- list(a,b,d,e,f,g,h,i)


## system.time(aa <- mat2adj(a,method="WGCNAFDR", FDR=1e-2, n.cores=6))
## aa <- mat2adj(a,method="bicorFDR", FDR=1e-2)
## bb <- mat2adj(b)

aa <- mat2adj(a, method="cor")
bb <- mat2adj(b, method="cor")
dd <- mat2adj(d, method="cor")
ee <- mat2adj(e, method="cor")
ff <- mat2adj(f, method="cor")
gg <- mat2adj(g, method="cor")
hh <- mat2adj(h, method="cor")
ii <- mat2adj(i, method="cor")


## diag(aa) <- 0
## diag(bb) <- 0
## write.table(aa, "net1.csv", row.names=FALSE,col.names=FALSE,sep="\t")
## write.table(bb, "net2.csv", row.names=FALSE,col.names=FALSE,sep="\t")
mylist <- lapply(mydlist, mat2adj, method="cor")

am <- "WGCNA"

sstab <- netSI(a, d="HIM", n.cores=1, save=FALSE)

# lapply(1:length(mylist),function(x,ga,h,adj.method, P){netdist(mylist[[x]],h=h, d="HIM", n.cores=1, adj.method=adj.method, P=P)}, ga=myga, h=aa, adj.method=am, P=4)

## mylist <- list(aa,bb,dd, ee, ff, gg, hh, ii)
names(mylist) <- paste(letters[c(1:2,4:9)],letters[c(1:2,4:9)],sep="")
names(mylist) <- NULL

## for (i in 1:10)
DDD <- netdist(mylist,d="HIM",n.cores=1, components=TRUE, verbose=TRUE)
DDDr <- netdist(mylist,d="HIM",n.cores=6, components=TRUE, verbose=FALSE, rho=100)

ssa <- netSI(a,d="HIM",n.cores=6,adj.method="MINE")
ssa <- lapply(mylist,netSI,adj.method="cor",method="kCV",k=4,h=10,n.cores=1)

## mydf <- data.frame(RHO=c(0,10**seq(0,10)),DIST=0)
## mydf[r+2,2] <- netdist(aa,bb,d="HIM", n.cores=1, components=FALSE, verbose=TRUE,rho=0)
## for (r in seq(0,10))
##   mydf[r+2,2] <- netdist(aa,bb,d="HIM", n.cores=1, components=FALSE, verbose=TRUE,rho=10**r)

## library(lattice)
## pp <- xyplot(DIST~RHO,data=mydf, scales=list(x=list(log=TRUE)),pch=19,
##              panel=function(...){
##                y1 <- netdist(aa,bb,d="hamming")
##                y2 <- netdist(aa,bb,d="ipsen")
##                panel.abline(h=y1, col="red")
##                panel.abline(h=y2, col="green")
##                panel.text(x=c(1,1),y=c(y1,y2),labels=c("hamming","ipsen"),adj=c(0,1),offset=0.8)
##                panel.xyplot(...)
##              },ylim=c(0.01,0.25),
##              xlab=expression(rho),ylab="HIM"
##              )

## pdf("test_rho.pdf")
## print(pp)
## dev.off()

system.time(netdist(mylist,d="ipsen"))
system.time(netdist(mylist,d="ipsen", n.cores=1))
system.time(netdist(mylist,d="ipsen", n.cores=4))


netdist(aa,bb,d="HIM", n.cores=1, components=FALSE)

aas <- aa + 2
bbs <- bb + 2

netdist(aa,bb, ga=0.1, d="ipsen", components=FALSE)
netdist(aa,bb)
## netdist(aas,bbs)
a[,10] <- 0
aa <- netSI(a,  indicator="all", sseed=10, adj.method="MINE",d="ipsen",save=TRUE, n.cores=4)
bb <- netSI(a,  indicator="all", sseed=10, save=FALSE, n.cores=4, ga=0.4)
save.image(file="computation.RData")

