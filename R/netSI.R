netSI <- function(d,indicator="all", dist='HIM', adj.method='cor', 
                  method="montecarlo", k=3, h=20, n.cores=NULL,save=TRUE, ...){
  
  ##to handle arguments passed through the dots we need to change a little bit here... 
  ##see coxph function of package survival
  Call <- match.call()
  id.Call <- match(c("d", "indicator", "dist", "adj.method","method","k","h","n.cores"), 
              names(Call), nomatch=0)
  if(id.Call[1]==0){
    stop("A dataset must be provided",call.=FALSE)
  }
  
  ##more time needed here
  ##----------------------------------------------
  
  INDICATORS <- c('S','SI','Sw','Sd',"all")
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator")
  if(indicator == -1)
    stop("ambiguous indicator")
  
  ##to be set by the user???
  sseed <- 0
  set.seed(sseed)
  
  ddim <- nrow(d)
  
  idxs <-  resamplingIDX(ddim,method=method, k=k, h=h)
  
  ##hardcoded the length of the list for optimization purposes
  ADJcv <- vector("list",length=length(idxs)+1)

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
    ##not working at the moment... possibly there are problems with the objects evaluated on the slaves
    cl <- makeCluster(getOption("cl.cores",n.cores))
    clusterEvalQ(cl,{
      d <- nettools:::d
      idxs <- nettools:::idxs
      mat2adj <- nettools:::mat2adj
      adj.method <-  nettools:::adj.method
      adj.measure <- nettools:::adj.measure
    })
    ##here we probably have to eval other parameters passed to the mat2adj function in case they are used... 
    ##probably we have to do some checkings on the arguments passed
    ADJcv <- clusterApply(cl,1:length(idxs),FUN=function(x,...){
      ss <- d[idxs[[x]],]
      tmp <- mat2adj(ss,method=adj.method,measure=adj.measure,...)
      return(tmp)
    })
    stopCluster(cl)
  } else{  
    ##we need to check that all the arguments of mat2adj are passed through the lapply
  ADJcv <- lapply(1:length(idxs),FUN=function(x){
    ss <- d[idxs[[x]],]
    tmp <- mat2adj(ss,method=adj.method,measure=adj.measure,...)
    return(tmp)})
  }
  #-----------------

  ##computing the adjacency matrix on the whole dataset
  ADJcv[['all']] <- mat2adj(x=d,method=adj.method,measure=adj.measure,...)

  
  #here the computation of the stability indicators is still missing...
  
  if(save==TRUE)
    results <- list("ADJlist"=ADJcv)
  else
    results <- list()
  
  return(results)
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
