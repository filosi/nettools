## Hamming distance for undirected graph
ham.undir <- function(object, ...){
  n <- object$N
  return(sum(abs(object$G1-object$G2))/(n*(n-1)))
}

## Hamming distance for directed graph
ham.dir <- function(object, ...){
  n <- object$N
  return(sum(abs(object$G1-object$G2))/(2*n*(n-1)))
}
