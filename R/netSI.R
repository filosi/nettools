netSI <- function(d,indicator="S", dist='HIM', adj.method='cor', adj.measure=NULL,
                  method="montecarlo", k=3, h=20, n.cores=NULL, ...){
  

  INDICATORS <- c('S','SI','Sw','Sd')
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator")
  if(indicator == -1)
    stop("ambiguous indicator")
  
  ##to be set by the user???
  sseed <- 0
  set.seed(sseed)
  
  ddim <- nrow(d)
  
  ##for all the indicators but the first one we have to perform resampling
  if(indicator!=1L){
    idxs <-  resamplingIDX(ddim,method=method, k=k, h=h)
    ADJcv <- vector(list,length=length(idxs))
    
    ##this cycle must be parallelized
    for (H in 1:length(idxs)){
      i <- idxs[[H]]
      ADJcv[[paste("res",H,sep="")]] <- mat2adj(x=d[i,],method=adj.method,measure=adj.measure,...)
      ##here we have to insert a check on NULL ADJs due to variance check
    }
    if(is.null(n.cores)){
      if(detectCores()>=2){
        n.cores <- detectCores()-1
      warning("The computation has been automatically parallelized")
      }
      else{
        n.cores <- 1
      }
    }
    if(n.cores>1){
      cl <- makeCluster(getOption("cl.cores",n.cores))
      clusterEvalQ(cl,{mat2adj <- nettools:::mat2adj})
      ll <- clusterApply(cl,1:length(idxs),FUN=mat2adj,)
      stopCluster(cl)
    } else {
      ll <- lapply(list(object$L1,object$L2),function(x,mygamma=optgamma,...){
        myomega <- sqrt(abs(round(spec(x),5)))
        myk <- K(mygamma,myomega)
        return(list(myomega,myk))
      })
    }
  
  }else{
    ##computing the adjacency matrix on the whole dataset
    ADJcv[['all']] <- mat2adj(x=d,method=adj.method,measure=adj.measure,...)
  }


}


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
