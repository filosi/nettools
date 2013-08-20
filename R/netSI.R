netSI <- function(d,indicator="S", dist='HIM', adj.method='cor', adj.measure=NULL,
                  method="montecarlo", k=3, h=20, n.cores=1, ...){
  

  INDICATORS <- c('S','SI','Sw','Sd')
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator")
  if(indicator == -1)
    stop("ambiguous indicator")
  
  method <- pmatch(method, METHODS)
  
                                        #to be set by the user???
  sseed <- 0
  set.seed(sseed)
  
  ADJcv <- vector(list,length=)
  ddim <- dim(d)[1]
  
  
                                        #for all the indicators but the first one we have to perform resampling
  if(indicator!=1L){
    indexes <-  resampling.index 
  }else{
                                        #computing the adjacency matrix on the whole dataset
    ADJcv[['all']] <- mat2adj(x=d,method=adj.method,measure=adj.measure,...)
    ## if(ADJcv[['all']])
  }
  
  
  ##verificare la definizione di k

  ## Computation!!
  for (H in 1:length(take)){
    ti <- take[[H]]
    ADJcv[[paste("res",H,sep="")]] <- mat2adj(x=d[ti,],method=adj.method,measure=adj.measure,...)
  }
}


resampling.index <- function(N,method="montecarlo", k=3, h=20){
  
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
