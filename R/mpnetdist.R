mpnetdist <- function(x, y, ga=NULL, ...){

  ## x, y should be lists
  ## x, y <- list(I = matrix(#layers x #layers), L = array(dim= (#nodes, #nodes, #layers)))

  ## Check x and y dimension.
  ## x,y should have compatible dimension
  n <- ifelse((dim(x[[2]])[1] == dim(x[[2]])[2]) && (dim(y[[2]])[1] == dim(y[[2]])[2]), dim(x$L)[1], 0)
  l <- ifelse((dim(x[[1]])[1] == dim(x[[1]])[2]) && (dim(y[[1]])[1] == dim(y[[1]])[2]), dim(x[[1]])[1], 0)
  if (!(n & l)){
    stop("Not conformable input arrays!")
  }
  
  ## Compute laplacian for inter-layer matrix  
  lx <- Lap(x[[1]])
  ly <- Lap(y[[1]])

  Lix <- directProd(lx, n)
  Liy <- directProd(ly, n)
  
  ## Compute intra-layer Laplacian
  ## Compute laplacian for each layer
  Llx <- intraLaplacian(x[[2]])
  Lly <- intraLaplacian(y[[2]])
  
  ## Compute the supra-laplacian
  Lapx <- Llx + Lix
  Lapy <- Lly + Liy
  
  ## Create the correct object for ipsen computation
  mylap <- list(L=list(Lapx, Lapy), N=n*l, tag="undir")

  ## Compute the ipsen distance
  mydist <- ipsen(mylap, ga=ga, ...)

  return(mydist)
}
