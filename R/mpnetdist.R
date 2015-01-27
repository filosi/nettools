mpnetdist <- function(x, y, d="HIM", ga=NULL, components=TRUE, ...){

  ## x, y should be lists
  ## x, y <- list(I = matrix(#layers x #layers), L = array(dim= (#nodes, #nodes, #layers)))

  ## Check x and y dimension.
  ## x,y should have compatible dimension
  n <- ifelse((dim(x[[2]])[1] == dim(x[[2]])[2]) && (dim(y[[2]])[1] == dim(y[[2]])[2]), dim(x$L)[1], 0)
  l <- ifelse((dim(x[[1]])[1] == dim(x[[1]])[2]) && (dim(y[[1]])[1] == dim(y[[1]])[2]), dim(x[[1]])[1], 0)
  if (!(n & l)){
    stop("Not conformable input arrays!")
  }


  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L
  
  ## Check distance type
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)
  
  ##check if need to return all components
  if(is.null(components)){
    if(d==1){
      comp <- TRUE
    } else {
      comp <- FALSE
    }
  } else {
    comp <- components
    if(d==1){
      if(!is.logical(comp))
        stop("components must be TRUE or FALSE")
    } else {
      comp <- FALSE
      warning("components parameter will be ignored", call. = FALSE)
    }
  }

  ##check on ga pass through ipsen function
  if(is.null(ga)){
    if(d==2){
      warning("The ga parameter will be automatically defined.", call.=FALSE)
    }
  }else{
    if(!is.numeric(ga) && !is.null(ga))
      stop("ga must be numeric",call.=FALSE)
  }


  ## Distance computation
  ##------------------------------
  
  
  ## Compute the supra-Laplacian
  Lapx <- mpLap(x, n)
  Lapy <- mpLap(y, n)
  mylap <- list(L=list(Lapx, Lapy), N=n*l, tag="undir")
  
  ## Compute the supra-Adjacency
  Adjx <- mpAdj(x, n)
  Adjy <- mpAdj(y, n)
  myadj <- list(G=list(Adjx, Adjy), N=n*l, tag="undir")

  ## Check the distance chosen and compute it
  ## HIM
  if (DISTANCE[d]=="HIM"){
    dd <- him(list(ADJ=myadj,LAP=mylap),  ga=ga,  components=comp, ltag=FALSE, ...)
  }
  ## IM
  if (DISTANCE[d]=="IM"){
    dd <- ipsen(mylap,ga=ga, ...)
  }
  ## H
  if (DISTANCE[d]=="H"){
    dd <- hamming(myadj)
  }

  ## Return distance list
  return(dd)
}



## This function computes the Laplacian for a multiplex network
mpLap <- function(x, nodes, ...){
  lx <- Lap(x[[1]])
  Lix <- directProd(lx, nodes)
  
  ## Compute intra-layer Laplacian
  ## Compute laplacian for each layer
  Llx <- intraLaplacian(x[[2]])
  
  ## Compute the supra-laplacian
  Lapx <- Llx + Lix

  return(Lapx)
}

## This function computes the Supra-Adj for a multiplex network
mpAdj <- function(x, nodes, ...){
  ## Compute directProduct for inter-layer connections
  Aix <- directProd(x[[1]], nodes)
  
  ## Supra adjacency for the Layers
  Alx <- IntraAdj(x[[2]])

  ## Supra-adjacency for the whole network
  Adjx <- Alx + Aix

  return(Adjx)
}
