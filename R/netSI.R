# debug(netSI)
# a <- matrix(rnorm(1000),ncol=100)
# b <- matrix(rnorm(1000),ncol=100)
# netSI(a,adj.method="MINE",measure="MIC",alpha=1)
# netSI(a,adj.method="cor",method="LOO")

netSI <- function(d,indicator="all", dist='HIM', adj.method='cor', 
                  method="montecarlo", k=3, h=20, n.cores,save=TRUE, verbose=TRUE,...){
  
  Call <- match.call()
  id.Call <- match(c("d", "indicator", "dist", "adj.method","method","k","h","n.cores"), 
                   names(Call), nomatch=0)
  if(id.Call[1]==0){
    stop("A dataset must be provided",call.=FALSE)
  }
  
  
  ##retrieving the arguments passed through ...
  ##extraArgs <- list(...)
  ##are not needed in the lapply calls below...
  
  INDICATORS <- c('S','SI','Sw','Sd',"all")
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator")
  if(indicator == -1)
    stop("ambiguous indicator")
  
  ##check on save and verbose are missing!!!!

  ##to be set by the user???
  sseed <- 0
  set.seed(sseed)
  
  ddim <- nrow(d)
  
  if(verbose==TRUE) cat("computing resampling...\n")
  idxs <-  resamplingIDX(ddim,method=method, k=k, h=h)
  
  ##hardcoded the length of the list for optimization purposes
  ADJcv <- vector("list",length=length(idxs))

  ## Call$n.cores should be evaluated first
  n.cores <- eval(Call$n.cores)
  ##check if it is user-friendly this way
  if(is.null(n.cores)){
    if(detectCores()>=2){
      n.cores <- detectCores()-1
      cl <- makeCluster(n.cores)
      warning("The computation has been automatically parallelized", call.=FALSE)
    }
    else{
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
  
  if(!is.null(cl)){
    ##from our checks the argument passed through ... are really used for 
    ##building the adjacency matrices during the parallel computation
    ADJcv <- parLapply(cl=cl,X=1:length(idxs),fun=function(x,d,method,...){
      ss <- d[idxs[[x]],]
      tmp <- nettools:::mat2adj(ss,...)
      return(tmp)
    },d=d,method=adj.method,...)
  } 
  else{
    ADJcv <- lapply(X=idxs,FUN=function(x,d,method,...){
      ss <- d[x,]
      tmp <- mat2adj(ss,...)
      return(tmp)
    },d=d,method=adj.method,...)
  }

  ##computing the adjacency matrix on the whole dataset
  ADJall <- mat2adj(x=d,method=adj.method,...)
  
  
  #here the computation of the stability indicators
  netsi <- list()
  if(indicator==1L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator S...\n")
    netsi[["S"]] <- netsiS(ADJall,ADJcv,dist=dist,cl=cl)
  }
  if(indicator==2L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator SI...\n")
    netsi[["SI"]] <- netsiSI(ADJcv,dist=dist,cl=cl)
  }
  if(indicator==3L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sw...\n")
    netsi[["Sw"]] <- netsiSw(ADJcv,cl=cl)
  }
  if(indicator==4L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sd...\n")
    netsi[["Sd"]] <- netsiSd(ADJcv,cl=cl)
  }
  
  if(save==TRUE)
    results <- list("call"=Call,"ADJlist"=ADJcv,
                    "S"=netsi[["S"]],"SI"=netsi[["SI"]],"Sw"=netsi[["Sw"]],"Sd"=netsi[["Sd"]])
  else
    results <- list("S"=netsi[["S"]],"SI"=netsi[["SI"]],"Sw"=netsi[["Sw"]],"Sd"=netsi[["Sd"]])

  if (!is.null(cl))
    stopCluster(cl)
  
  return(results)
}

##need to do some checkings in order to pass also gamma through...
netsiS <- function(g,H,dist,cl){
  
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 2
  
  if(!is.null(cl)){
    s <- parLapply(cl=cl,X=H,fun=function(x,g,dist,type){
      res <- nettools:::netdist(g,x,dist)[[type]]
      return(res)
    },g=g,dist=dist,type=type)
  }else{
    s <- lapply(X=H,FUN=function(x,g,dist,type){
      res <- netdist(g,x,dist)[[type]]
      return(res)
    },g=g,dist=dist,type=type)
  }
  return(unlist(s))
}

netsiSI <- function(H,dist,cl){
  
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 2
  
  com <- combn(1:length(H), 2)
  if(!is.null(cl)){
    s <- parLapply(cl=cl,X=1:ncol(com),fun=function(x,com,H,dist,type){
      res <- nettools:::netdist(H[[com[1,x]]],H[[com[2,x]]],dist)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type)
  }else{
    s <- lapply(X=1:ncol(com),FUN=function(x,com,H,dist,type){
      res <- netdist(H[[com[1,x]]],H[[com[2,x]]],dist)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type)
  }
  return(unlist(s))
}

## NB dist parameter not needed
netsiSd <- function(H,cl){
  if (length(H))
    n <- ncol(H[[1]])
  else
    stop("No adjacency matrix computed",call.=FALSE)
  
  if (!is.null(cl)){
    dd <- parLapply(cl=cl, X=H, rowSums)
  } else {
    dd <- lapply(H, rowSums)
  }
  dd <- matrix(unlist(dd),ncol=n)
  
  return(dd)
}

netsiSw <- function(H,cl){
  if (length(H))
    n <- nrow(H[[1]])
  else
    stop("List of adjacency matrices do not exist")
  
  com <- combn(1:n, 2)
  ## m <- matrix(NA,ncol=length(H),nrow=ncol(com))
  ## rownames(m) <- paste(com[1,],com[2,],sep="-")
  ## colnames(m) <- names(H)
  if (!is.null(cl)){
    tmp <- parLapply(cl, H,
                  function(x,com){
                    sapply(1:ncol(com),
                           function(y, x, allcom){
                             x[allcom[1,y],allcom[2,y]]},
                           x=x, allcom=com)
                  },
                  com=com)
  } else {
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
  
  ## for(i in 1:length(H)){
  ##   for(j in 1:ncol(com)){
  ##     m[j,i] <- H[[i]][com[1,j],com[2,j]]
  ##   }
  ## }
  return(m)
}
#   s <- lapply(X=1:length(H),FUN=function(x,com,H){
#     res <- lapply(X=1:ncol(com),FUN=function(y,com,G){
#       m[x,y] <- G[com[1,y],com[2,y]]
#       return(m)
#     },com=com,G=H[[x]])
#     return(res)
#   },com=com,H=H)
#   
#   return(s)
# }

resamplingIDX <- function(N,method="montecarlo", k=3, h=20){
  
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
  
  return(take)
  
}
