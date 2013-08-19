netSI <- function(d,indicator="S", dist='HIM', adj='cor', adj.method=NULL,
                    method="montecarlo", k=3, h=20, FDR=1e-3, n.cores=1, ...){

    METHODS <- c('montecarlo','LOO','kCV')
    method <- pmatch(method, METHODS)

    if(is.na(method))
      stop("invalid distance method")
    if(method == -1)
      stop("ambiguous distance method")
    
    INDICATORS <- c('S','SI','Sw','Sd')
    indicator <- pmatch(indicator,INDICATORS)
    
    if(is.na(indicator))
      stop("invalid indicator")
    if(indicator == -1)
      stop("ambiguous indicator")
    
    #to be set by the user???
    sseed <- 0
    set.seed(sseed)
    
    ADJcv <- list()
    ddim <- dim(d)[1]
    
    ##verificare la definizione di k
    ## Montecarlo
    if (method==1)
      take <- lapply(1:h,function(x){tmp <- sample(seq(ddim),floor(ddim*(1-(1/k)))) 
                                     return(tmp)})
    
    ## Leave One Out
    if (method==2){
      if (h!=1)
        warning("h is different than 1 but method is set to LOO (Leave One Out cross-validation schema).\nh will be ignored.")
      h <- ddim
      take <- lapply(1:ddim,function(x,allid){return(allid[which(allid!=x)])},allid=1:ddim)
    }
    
    ## K-fold cross-validation
    if (method==3){
      if (k>=ddim)
        stop("Number of fold bigger than samples in the dataset!!!")
      else{
        take <- list()
        for (H in seq(h)){
          tmpcv <- sample(seq(ddim),ddim)
          s <- split(tmpcv,cut(seq(ddim),k))
          names(s) <- NULL
          take <- c(take,lapply(s,function(x,idx){tmp <- idx[!is.element(idx,x)]
                                                  #cat("H =",H,"\n",tmp,fill=TRUE,file=filelogs,append=TRUE)
                                                  return(tmp)},idx=tmpcv))
        }
      }
    }
    
    
    ## Computation!!
    for (H in 1:length(take)){
      ti <- take[[H]]
      v0tmp <- which(apply(d[ti,],2,sd)<1e-5)
      if (length(v0tmp)>=1){
        var0.feat[names(v0tmp)] <- var0.feat[names(v0tmp)] + 1
        #warning(paste("Pre dim ", dim(d),"\n"),.call=FALSE)
        d <- d[,-v0tmp]
      }
      if (cvlab=="montecarlo" || cvlab=="loo"){
        ## Put ifelse in the cat command for printing of Variable with 0 variance
        flag <- ifelse(length(v0tmp)==0,v0tmplab <- "",v0tmplab <- paste("Var =",v0tmp,sep=" "))
        #cat("At resampling ",H, " algs ", a, " Classes ",cl,v0tmplab,file=filelogs,append=TRUE,fill=TRUE)
        if (cvlab=="montecarlo")
          cat(ti,fill=TRUE,file=filelogs,append=TRUE)
      } else {
        flag <- ifelse(length(v0tmp)==0,v0tmplab <- "",v0tmplab <- paste("Var =",v0tmp,sep=" "))
        cat("At resampling ",(H%/%k)+1," algs ", a, " Classes ",cl, v0tmplab,file=filelogs,append=TRUE,fill=TRUE)
      }
      ADJcv[[H]] <- myfun(d[ti,],P=1,FDR=fdr)
    }
}

