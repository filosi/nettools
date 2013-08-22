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
  
  ##check if it is user-friendly this way
  if(is.null(Call$n.cores)){
    if(detectCores()>=2){
      n.cores <- detectCores()-1
      warning("The computation of the adjacency matrices has been automatically parallelized")
    }
    else{
      n.cores <- 1
    }
  }
  
  if(verbose==TRUE) cat("computing adjacency matrices...\n")
  
  if(n.cores>1){
    
    cl <- makeCluster(n.cores)
    
    ##from our checks the argument passed through ... are really used for 
    ##building the adjacency matrices during the parallel computation
    ADJcv <- parLapply(cl=cl,X=1:length(idxs),fun=function(x,d,method,...){
      ss <- d[idxs[[x]],]
      tmp <- nettools:::mat2adj(ss,...)
      return(tmp)
    },d=d,method=adj.method,...)
    stopCluster(cl)
  } 
  else{
    ADJcv <- lapply(X=1:length(idxs),FUN=function(x,d,method,...){
      ss <- d[idxs[[x]],]
      tmp <- mat2adj(ss,...)
      return(tmp)
    },d=d,method=adj.method,...)
  }
  
  ##computing the adjacency matrix on the whole dataset
  ADJall <- mat2adj(x=d,method=adj.method,...)
  
  
  #here the computation of the stability indicators is still missing...
  netsi <- list()
  if(indicator==1L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator S...\n")   
    #fare il conto di tutte le distanze tra ADJcv[["all"]] e i restanti ADJcv
      netsi[["S"]] <- netsiS(ADJall,ADJcv,dist=dist,n.cores=n.cores)
  }
  if(indicator==2L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator SI...\n")
    #fare il conto di tutte le distanze tra ADJcv[["all"]] e i restanti ADJcv
    netsi[["SI"]] <- netsiSI(ADJcv,dist=dist,n.cores=n.cores)
  }
  if(indicator==3L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sw...\n")
    #fare il conto di tutte le distanze tra ADJcv[["all"]] e i restanti ADJcv
    netsi[["Sw"]] <- netsiSw(ADJcv,n.cores=n.cores)
  }
  if(indicator==4L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sd...\n")
    #fare il conto di tutte le distanze tra ADJcv[["all"]] e i restanti ADJcv
    netsi[["Sd"]] <- netsiSd(ADJcv,dist=dist,n.cores=n.cores)
  }
  
  if(save==TRUE)
    results <- list("call"=Call,"ADJlist"=ADJcv,
                    "S"=netsi[["S"]],"SI"=netsi[["SI"]],"Sw"=netsi[["Sw"]],"Sd"=netsi[["Sd"]])
  else
    results <- list("S"=netsi[["S"]],"SI"=netsi[["SI"]],"Sw"=netsi[["Sw"]],"Sd"=netsi[["Sd"]])
  
  return(results)
}

##need to do some checkings in order to pass also gamma through...
netsiS <- function(g,H,dist,n.cores){
  
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 2
  
  if(n.cores>1){
    cl <- makeCluster(n.cores)
    s <- parLapply(cl=cl,X=1:length(H),fun=function(x,g,H,dist,type){
      res <- nettools:::netdist(g,H[[x]],dist)[[type]]
      return(res)
    },g=g,H=H,dist=dist,type=type)
    stopCluster(cl)
  }else{
    s <- lapply(X=1:length(H),FUN=function(x,g,H,dist,type){
      print(x)
      res <- netdist(g,H[[x]],dist)[[type]]
      return(res)
    },g=g,H=H,dist=dist,type=type)
  }
  return(unlist(s))
}

netsiSI <- function(H,dist,n.cores){
  
  type <- pmatch(dist,c("H","IM","HIM","hamming","ipsen"))
  if(type==4L) type <- 1
  if(type==5L) type <- 2
  
  com <- combn(1:length(H), 2)
  if(n.cores>1){
    cl <- makeCluster(n.cores)
    s <- parLapply(cl=cl,X=1:ncol(com),fun=function(x,com,H,dist,type){
      res <- nettools:::netdist(H[[com[1,x]]],H[[com[2,x]]],dist)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type)
    stopCluster(cl)
  }else{
    s <- lapply(X=1:ncol(com),FUN=function(x,com,H,dist,type){
      res <- netdist(H[[com[1,x]]],H[[com[2,x]]],dist)[[type]]
      return(res)
    },com=com,H=H,dist=dist,type=type)
  }
  return(unlist(s))
}

netsiSw <- function(H,n.cores){
  com <- combn(1:nrow(H[[1]]), 2)
  m <- matrix(NA,nrow=length(H),ncol=ncol(com))
  colnames(m) <- seq(ncol(m))
  rownames(m) <- seq(nrow(m))

  for(j in 1:ncol(com)){
    for(i in 1:length(H)){
      m[i,j] <- H[[i]][com[1,j],com[2,j]]
      
    }
    colnames(m)[j] <- paste(com[1,j],com[2,j],sep="-")
  }
  rownames(m) <- names(H)
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

netsiSd <- function(H,n.cores){}


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
    if (h!=1)
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
