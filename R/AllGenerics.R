## Generic methods
mat2adj.default <- function(x,...){
  stop("Wrong input, data.frame or matrix is required")
}
setGeneric("mat2adj",mat2adj.default)

g2adj.default <- function(x,...){
  Adj <- apply(x,2,as.numeric)
  diag(Adj) <- 0
  if (!all(Adj - t(Adj)) == 0)
    stop("Not a symmetric matrix!")
  return(Adj)
}
setGeneric("g2adj",g2adj.default)

Lap.default <- function(x,...){
  D <- apply(x,2,sum)
  return((D * diag(dim(x)[1])) - x)
}
setGeneric("Lap",Lap.default)

ipsen.default <- function(object,...){
  warning("Not the correct input")
}
setGeneric("ipsen",ipsen.default)

hamming.default <- function(object,...){
  warning("Not the correctinput!")
}
setGeneric("hamming",hamming.default)

him.default <- function(object,...){
  warning("Not a correct input object!")
}
setGeneric("him",him.default)

