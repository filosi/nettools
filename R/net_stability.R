stab.core <- function(ADJcv,com=NULL,allcomb=FALSE,ADJ=NULL){
  ## NB ADJ should be already cutted for algorithm, class and pathway
  ##com <- combn(seq(H), 2)
  dist <- list()
  if (allcomb){
    for (cc in seq(dim(com)[2])){
      dist[[cc]] <- him(ADJcv[[com[,cc][1]]],ADJcv[[com[,cc][2]]])
    }
  } else {
        for (cc in 1:length(ADJcv)){
      dist[[cc]] <- him(ADJ,ADJcv[[cc]])
    }
  }
  return(dist)
}

net.stability <- function(i,wholelist,com=NULL,allcomb=FALSE,ADJ=NULL){
  ADJcv <- wholelist[[i]]
  nome <- names(wholelist)[i]
  algs <- names(ADJcv)
  dist <- list()
  for (a in algs){
    classes <- names(ADJcv[[a]])
    dist[[a]] <- list()
    for (cl in classes){
      if (allcomb){
        dist[[a]][[cl]] <- stab.core(ADJcv[[a]][[cl]],com=com,allcomb=allcomb,ADJ=ADJ)        
      }else{
        dist[[a]][[cl]] <- stab.core(ADJcv[[a]][[cl]],com=com,allcomb=allcomb,ADJ=ADJ[[nome]][[a]][[cl]])
      }
    }
  }
  return(dist)
}

#load('data.RData')
rdatafile <- '../results/adj_thr.RData'
load(rdatafile)
#load("ADJcv.RData")
library(combinat)
library(multicore)
library(WGCNA)
library(minet)
source('NetTools_v1.R')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## CROSS VALIDATION
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filelogs <- file("net_stability.log", "w")
cat("Starting Cross Validation!!",file=filelogs,append=TRUE,sep="\n")
cat("Inferring the weighted networks with CROSS VALIDATION...\n",file=filelogs,append=TRUE)

## Cross Validation Parameters
## cvlab: string which identifies the cross validation schema to be used for stability
##        By default it is set to <montecarlo>, possible values are <loo> for Leave One Out or
##        <kfold> for K-fold cross validation
cvlab <- "loo"

## k: numeric value indicates how many split on the data.
##    By default it is set to 3
##    In case cvlab=="montecarlo" the 1-1/cv samples are taken for the computation
##    In case cvlab=="loo" it should be set to 1
##    In case cvlab=="kcv" cv iteration are performed, each dividing the data into 1-1/cv
k <- 2

## h: numeric value indicates the number of iterations (cross validations)
##    By default it is set to 20
h <- 20

## How to compute distances:
## allcomb: logical defining how the distances should be computed.
##          By default it is set to FALSE, and the distance is computed
##          between each resampling and the network inferred using all the samples
##          If <allcomb>=TRUE the distance is computed between all the combination
##          in the resampling scheme.
allcomb <- FALSE

## Pathway selection
## pth.sel: By default = NULL (compute distance on all pathway in pat.list)
##          if pth.sel = c("pth1","pth2") compute the distance only
##          on the selected list.
pth.sel <- NULL

## ncores: numeric, number of cores to be used for the computation.
##         Warning: Final use of memory is ncores*use of memory!!!!
ncores <- 5

## FDR: the False Discovery Rate threshold above which the weight of a link is
##      considered reliable and not set to 0
fdr <- 1e-4

## exp.flag: A flag to be added to the result files to better identify
##           the experiment (example: "_FDR1e-4")
exp.flag <- "_AllcombFALSE"

## Write information to log file
cat("Resampling scheme:",cvlab,file=filelogs,append=TRUE,fill=TRUE)
cat("k: ",k,"\nh: ",h,file=filelogs,append=TRUE,fill=TRUE)
cat("Cores: ",ncores,file=filelogs,append=TRUE,fill=TRUE)
cat("Selected Pathways:\n",paste(pth.sel,collapse="\n"),file=filelogs,append=TRUE,fill=TRUE)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
## DO NOT MODIFY THE CODE BELOW!!!!!!!
##
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##--------------------------------------------------
## Start Computation
##--------------------------------------------------

  cat("Starting computing ADJ matrices!!\n")
  sseed <- 0
  set.seed(sseed)
  
  ADJcv <- list()
  var0.feat <- vector("double",length=length(features))
  names(var0.feat) <- features
  
  for(p in pat.list){
    cat(paste("..on pathway",p,"\n"), file=filelogs,append=TRUE)
    ##select the subset of features belonging to the pathway
    idx <- which(is.element(features,pat2[pat2[,2]==p,1]))
    ADJcv[[p]] <- list()
    for (a in algs){
    ##for (a in "WGCNAFDR"){
      ADJcv[[p]][[a]] <- list()
      myalg <- paste("Adj",toupper(a),sep="")
      myfun <- get(myalg,mode="function")
      for (cl in levels(classes)){
        ##select the subset corresponding to a class	    
        d <- states.data[classes==cl,idx]
        ddim <- dim(d)[1]
        ADJcv[[p]][[a]][[cl]] <- list()
        
        ## Montecarlo
        if (cvlab=="montecarlo")
          take <- lapply(1:h,function(x){tmp <- sample(seq(ddim),floor(ddim*(1-(1/k))))
                                         return(tmp)})
        
        ## Leave One Out
        if (cvlab=="loo"){
          if (h!=1)
            warning("h is different than 1 but you asked for Leave One Out cross-validation scheme.\nh will be ignored!!!")
          h <- ddim
          take <- lapply(1:ddim,function(x,allid){return(allid[which(allid!=x)])},allid=1:ddim)
        }
        
        ## K-fold cross-validation
        if (cvlab=="kfold"){
          if (k>=dim(d)[1])
            stop("Number of fold bigger than samples in the dataset!!!")
          else{
            take <- list()
            for (H in seq(h)){
              tmpcv <- sample(seq(ddim),ddim)
              s <- split(tmpcv,cut(seq(ddim),k))
              names(s) <- NULL
              take <- c(take,lapply(s,function(x,idx){tmp <- idx[!is.element(idx,x)]
                                                      cat("H =",H,"\n",tmp,fill=TRUE,file=filelogs,append=TRUE)
                                                      return(tmp)},idx=tmpcv))
            }
          }
        }
        
        H <- 0
        ## Computation!!
        for (ti in take){
          v0tmp <- which(apply(d[ti,],2,sd)<1e-5)
          if (length(v0tmp)>=1){
            var0.feat[names(v0tmp)] <- var0.feat[names(v0tmp)] + 1
            cat(paste("Pre dim ", dim(d),"\n"),file=filelogs,append=TRUE)
            d <- d[,-v0tmp]
          }
          if (cvlab=="montecarlo" || cvlab=="loo"){
            ## Put ifelse in the cat command for printing of Variable with 0 variance
            flag <- ifelse(length(v0tmp)==0,v0tmplab <- "",v0tmplab <- paste("Var =",v0tmp,sep=" "))
            cat("At resampling ",H, " algs ", a, " Classes ",cl,v0tmplab,file=filelogs,append=TRUE,fill=TRUE)
            if (cvlab=="montecarlo")
              cat(ti,fill=TRUE,file=filelogs,append=TRUE)
          } else {
            flag <- ifelse(length(v0tmp)==0,v0tmplab <- "",v0tmplab <- paste("Var =",v0tmp,sep=" "))
            cat("At resampling ",(H%/%k)+1," algs ", a, " Classes ",cl, v0tmplab,file=filelogs,append=TRUE,fill=TRUE)
          }
          ADJcv[[p]][[a]][[cl]][[(H+1)]] <- myfun(d[ti,],P=1,FDR=fdr)
          H <- H + 1
        }
      }
    }
  }

  save(list=c("ADJcv","allcomb","cvlab","k","h","pth.sel"),
       file=paste("../results/ADJcv_",cvlab,"_h",h,"_k",k,"_FDR",fdr,exp.flag,".RData",sep=""))

  cat("Saved ADJcv!\n")

if (allcomb){
  cat("MErda\n")
  flag <- ifelse(cvlab=="kfold",com <- combn(seq(h*k), 2),com <- combn(seq(h), 2))
  if (!is.null(pth.sel)){
    adjidx <- which(is.element(names(ADJcv),pth.sel))
    adjnames <- pth.sel
  } else {
    adjidx <- seq(length(ADJcv))
    adjnames <- names(ADJcv)
  }
  DistCV_allcomb<- mclapply(adjidx,FUN=net.stability,
                     wholelist=ADJcv,
                     com=com,allcomb=TRUE,ADJ=NULL,
                     mc.cores=ncores)
  names(DistCV_allcomb) <- adjnames
} else {
  if (!is.null(pth.sel)){
    adjidx <- which(is.element(names(ADJcv),pth.sel))
    adjnames <- pth.sel
  } else {
    adjidx <- seq(length(ADJcv))
    adjnames <- names(ADJcv)
  }
  DistCV <- mclapply(adjidx,FUN=net.stability,
                     wholelist=ADJcv,
                     com=NULL,allcomb=FALSE,ADJ=ADJ,
                     mc.cores=ncores)
  names(DistCV) <- adjnames
}

##save(list=c("DistCV","DistCV_allcomb","ADJcv","allcomb","cvlab","k","h","pth.sel"),
##      file=paste("../results/distances_overCV_",cvlab,"_h",h,"_k",k,"_FDR",fdr,exp.flag,".RData",sep=""))
save(list=c("DistCV","ADJcv","allcomb","cvlab","k","h","pth.sel"),
      file=paste("../results/distances_overCV_",cvlab,"_h",h,"_k",k,"_FDR",fdr,exp.flag,".RData",sep=""))
##save.image("net_stability.RData")
cat("Done!",file=filelogs,append=TRUE)
close(filelogs)
