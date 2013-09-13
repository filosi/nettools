netSI <- function(d,indicator="all", dist='HIM', adj.method='cor', 
                  method="montecarlo", k=3, h=20, n.cores,save=TRUE,
                  verbose=TRUE, ...){

  ## Get the function call parameters
  ## NB the parameters should be evaluated with eval()
  Call <- match.call()
  id.Call <- match(c("d", "indicator", "dist", "adj.method","method","k","h","n.cores"), 
                   names(Call), nomatch=0)
  if(id.Call[1]==0){
    stop("A dataset must be provided",call.=FALSE)
  }else{
    d <- eval(Call$d)
    if(!(is.matrix(d) | is.data.frame(d))){
      stop("d must be a matrix or a data frame", call.=FALSE)
    }
  }
    
  ## Choose the indicators
  if(id.Call[2]!=0){
    indicator <- eval(Call$indicator)
    if(!is.character(indicator))
      stop("indicator must be one of 'S','SI','Sw','Sd','all'.", call.=FALSE)
  }
  
  INDICATORS <- c('S','SI','Sw','Sd',"all")
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator", call. =FALSE)
  if(indicator == -1)
    stop("ambiguous indicator", call. =FALSE)
  
  ##check on save and verbose
  if(is.null(Call$save)){
    save <- TRUE
  }else{
    save <- eval(Call$save)
    if(!is.logical(save))
      stop("save must be TRUE or FALSE",call.=FALSE)
  }
  
  if(is.null(Call$verbose)){
    verbose <- TRUE
  }else{
    verbose <- eval(Call$verbose)
    if(!is.logical(verbose))
      stop("verbose must be TRUE or FALSE", call.=FALSE)
  }

  ## It can be passed through ...
  if (is.null(Call$sseed)) sseed <- 0
  else sseed <- eval(Call$sseed)
  set.seed(sseed)

  ## Pass parameter gamma to netdist functions
  if(!is.null(ga)){
    ga <- eval(Call$ga)
  }
  
  ## Pass parameter components to netdist functions
  if(!is.null(components)){
    components <- eval(Call$components)
    if(dist=="HIM" & components==TRUE){
      warning(paste("component parameter will be ignored. \n
            The stability indicators will be computed just for",
                    dist, "distance.\n
            For computing them for Hamming or Ipsen-Mikhailov 
            distance use dist=hamming or dist=ipsen"), call.=FALSE)
      components <- FALSE
    }
    
  ## Get the dimension of the input matrix
  ddim <- nrow(d)

  ## Get the resampling indexes
  if(verbose==TRUE) cat("computing resampling...\n")
  idxs <-  resamplingIDX(ddim,method=method, k=k, h=h)
  
  ## length of the list for optimization purposes
  ADJcv <- vector("list",length=length(idxs))

  ## evaluate the number of cores to be used (Call$n.cores)
  n.cores <- eval(Call$n.cores)
  
  ## Check availability of cores, otherwise set a default
  if(is.null(n.cores)){
    if(detectCores()>=2){
      n.cores <- detectCores()-1
      cl <- makeCluster(n.cores)
      warning("The computation has been automatically parallelized", call.=FALSE)
    } else {
      cl <- NULL
    }
  } else {
    if (n.cores==1){
      cl <- NULL
    } else {
      if (n.cores<detectCores()){
        cl <- makeCluster(n.cores)
      } else {
        if(detectCores()>=2){
          n.cores <- detectCores()-1
          cl <- makeCluster(n.cores)
          warning("The computation has been automatically parallelized", call.=FALSE)
        } 
      }
    }
  }
  
  if(verbose==TRUE) cat("computing adjacency matrices...\n")

  ## Compute the adjacency matrices for each resampling index
  if(!is.null(cl)){
    ## Parallel computation
    ADJcv <- parLapply(cl=cl,X=idxs,fun=function(x,d,method,...){
      ss <- d[x,]
      tmp <- nettools:::mat2adj(ss, method=method, ...)
      return(tmp)
    },d=d,method=adj.method,...)
  } else {
    ## One core computation
    ADJcv <- lapply(X=idxs,FUN=function(x,d,method,...){
      ss <- d[x,]
      tmp <- mat2adj(ss, method=method, ...)
      return(tmp)
    },d=d,method=adj.method,...)
  }
  
  ##computing the adjacency matrix on the whole dataset
  ADJall <- mat2adj(x=d,method=adj.method,...)
    
  ## Compute the stability indicators
  netsi <- list()
  if(indicator==1L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator S...\n")
    netsi[["S"]] <- netsiS(ADJall, ADJcv, dist=dist, cl=cl, ...)
  }
  if(indicator==2L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator SI...\n")
    netsi[["SI"]] <- netsiSI(ADJcv, dist=dist, cl=cl, ...)
  }
  if(indicator==3L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sw...\n")
    netsi[["Sw"]] <- netsiSw(ADJcv, cl=cl)
  }
  if(indicator==4L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sd...\n")
    netsi[["Sd"]] <- netsiSd(ADJcv, cl=cl)
  }
  
  if(save==TRUE){
    results <- list("call"=Call,"ADJlist"=ADJcv,
                    "S"=netsi[["S"]],
                    "SI"=netsi[["SI"]],
                    "Sw"=netsi[["Sw"]],
                    "Sd"=netsi[["Sd"]])
  } else {
    results <- list("S"=mean(netsi[["S"]]),
                    "SI"=mean(netsi[["SI"]]),
                    "Sw"=apply(netsi[["Sw"]], 1, mean, na.rm=TRUE),
                    "Sd"=apply(netsi[["Sd"]], 2, mean, na.rm=TRUE)
                    )
  }
  
  if (!is.null(cl))
    stopCluster(cl)
  
  return(results)
}


## Stability indicator S
netsiS <- function(g,H,dist,cl, ...){
  
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 1
  
  if(!is.null(cl)){
    s <- parLapply(cl=cl,X=H,fun=function(x,g,dist,type,  ...){
      res <- nettools:::netdist(g,x,dist, ...)[[type]]
      return(res)
    },g=g,dist=dist,type=type, ...)
  }else{
    s <- lapply(X=H,FUN=function(x,g,dist,type,...){
      res <- netdist(g,x,dist, ...)[[type]]
      return(res)
    },g=g,dist=dist,type=type,...)
  }
  return(unlist(s))
}

## Stability indicator SI
netsiSI <- function(H,dist,cl, ...){
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 1

  com <- combn(1:length(H), 2)
  
  if(!is.null(cl)){
    ## Parallel computation
    s <- parLapply(cl=cl,X=1:ncol(com),fun=function(x,com,H,dist,type, ...){
      res <- nettools:::netdist(H[[com[1,x]]],H[[com[2,x]]],dist, ...)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type, ...)
  }else{ ## One core computation
    s <- lapply(X=1:ncol(com),FUN=function(x,com,H,dist,type, ...){
      res <- netdist(H[[com[1,x]]],H[[com[2,x]]],dist, ...)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type, ...)
  }
  return(unlist(s))
}

## Degree stability
netsiSd <- function(H,cl){
  if (length(H))
    n <- ncol(H[[1]]) else stop("No adjacency matrix computed",call.=FALSE)
  
  if (!is.null(cl)){
    ## Parallel computation
    dd <- parLapply(cl=cl, X=H, rowSums)
  } else {
    ## One core computation
    dd <- lapply(H, rowSums)
  }
  dd <- matrix(unlist(dd),ncol=n)
  return(dd)
}

## Edges stability
netsiSw <- function(H,cl){
  if (length(H))
    n <- nrow(H[[1]]) else stop("List of adjacency matrices do not exist")
  
  com <- combn(1:n, 2)

  if (!is.null(cl)){
    ## Parallel computation
    tmp <- parLapply(cl, H,
                  function(x,com){
                    sapply(1:ncol(com),
                           function(y, x, allcom){
                             x[allcom[1,y],allcom[2,y]]},
                           x=x, allcom=com)
                  },
                  com=com)
  } else {
    ## One core computation
    tmp <- lapply(H,
                  function(x,com){
                    sapply(1:ncol(com),
                           function(y, x, allcom){
                             x[allcom[1,y],allcom[2,y]]},
                           x=x, allcom=com)
                  },
                  com=com)
  }

  ## Set up the results in a matrix
  m <- matrix(unlist(tmp), ncol=length(H))
  rownames(m) <- paste(com[1,],com[2,],sep="-")
  colnames(m) <- names(H)
  
  return(m)
}

## Function for the computation of resampling indexes
resamplingIDX <- function(N,method="montecarlo", k=3, h=20){

  ## Check of resampling methods
  METHODS <- c('montecarlo','LOO','kCV')
  method <- pmatch(method, METHODS)

  if(is.na(method))
    stop("invalid distance method")
  if(method == -1)
    stop("ambiguous distance method")
  
  ## Montecarlo
  if (method==1L)
    take <- lapply(1:h,function(x){tmp <- sample(seq(N),floor(N*(1-(1/k)))) 
                                   return(tmp)})
  
  ## Leave One Out
  if (method==2L){
    if (k!=1)
      warning("h is different than 1 but method is set to LOO (Leave One Out cross-validation schema).\nh will be ignored.")
    h <- N
    take <- lapply(1:N,function(x,allid){return(allid[which(allid!=x)])},allid=1:N)
  }
  
  ## K-fold cross-validation
  if (method==3L){
    if (k>=N)
      stop("Number of fold bigger than samples in the dataset!!!")
    else{
      take <- list()
      for (H in seq(h)){
        tmpcv <- sample(seq(N),N)
        s <- split(tmpcv,cut(seq(N),k))
        names(s) <- NULL
        take <- c(take,lapply(s,function(x,idx){tmp <- idx[!is.element(idx,x)]
                                                return(tmp)},idx=tmpcv))
      }
    }
  }
  
  ## return a list with indexes
  return(take)
}
