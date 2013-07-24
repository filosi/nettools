spec <- function(matrix){
  sort(eigen(matrix)$values)
}

## D2 - ipsen02evolutionary
lorentz <- function(omega,mygamma,given_omega){
  l <-0
  for(i in 2:length(given_omega)){
    l = l + mygamma/( (omega-given_omega[i])**2+mygamma**2)                          }
  return(l)
}

K <- function(mygamma,given_omega){
  return(1/integrate(lorentz,lower=0,upper=Inf,mygamma=mygamma,given_omega=given_omega, stop.on.error = FALSE)$value)
}

rho <- function(omega, mygamma, ll){
  ll[[2]]*lorentz(omega,mygamma,ll[[1]])
}

ipsen_minus_one <- function(g,n){
  return(sqrt(
    1/(pi*g) +
      1/(2*g*(atan(sqrt(n)/g)+pi/2)**2)*(pi/2+ (sqrt(n)/g)/(1+(sqrt(n)/g)**2)+atan(sqrt(n)/g))
    -4*(pi-(g/sqrt(n))*log(1/(1+(sqrt(n)/g)**2))+atan(sqrt(n)/g))/
      (pi*g*(4+(sqrt(n)/g)**2)*(atan(sqrt(n)/g)+pi/2))
  )
         -1)
}

optimal_gamma <- function(n){
  return(uniroot(ipsen_minus_one,c(0.01,1),n=n,maxiter=100000,tol=.Machine$double.eps)$root)
}


##--------------------------------------------------
## DIRECTED GRAPH DISTANCE
##--------------------------------------------------

ipsen_minus_one_dir  <- function(g,n){
  return(ZZ(g)^2*MM(0,g)+WW(n,g)^2*MM(n-2,g)+WW(n,g)^2*MM(n,g)+WWp(n,g)^2*MM(2*n-2,g)
         -2*ZZ(g)*WW(n,g)*LL(0,n-2,g)-2*ZZ(g)*WW(n,g)*LL(0,n,g)-2*ZZ(g)*WWp(n,g)*LL(0,2*n-2,g)
         +2*WW(n,g)*WW(n,g)*LL(n-2,n,g) +2*WW(n,g)*WWp(n,g)*LL(n-2,2*n-2,g) +2*WW(n,g)*WWp(n,g)*LL(n,2*n-2,g)-1)
}

optimal_gamma_dir <- function(n){
  return(uniroot(ipsen_minus_one_dir,c(0.01,1),n=n,maxiter=100000,tol=.Machine$double.eps)$root)
}


MM <- function(T,g){
  return(1/2*(g^2*atan(sqrt(T)/g) + T*atan(sqrt(T)/g) + sqrt(T)*g)/(g^5 +T*g^3) + 1/4*pi/g^3)
}

LL  <- function(T,U,g){
  return(-log(g^2 + U)/((4*g^2 + T + 3*U)*sqrt(T) - (4*g^2 + 3*T + U)*sqrt(U)) + log(g^2 + T)/((4*g^2 + T + 3*U)*sqrt(T) - (4*g^2 + 3*T + U)*sqrt(U)) + pi/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g) + atan(sqrt(T)/g)/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g) + atan(sqrt(U)/g)/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g))
}

ZZ  <- function(g){
  return(2*g/pi)
}

WW  <- function(N,g){
  WWW  <- (2*N-1)*(pi/2)+(N-1)*(atan(sqrt(N-2)/g) +atan(sqrt(N)/g)) + atan(sqrt(2*N-2)/g)
  return(g*(N-1)/WWW)
}

WWp  <- function(N,g){
  return(WW(N,g)/(N-1))
}










## NB to develop

## HIM Distance for directed graphs
## undir  <- function(A){
##   n  <- dim(A)
##   zero  <- matrix(0,nrow=n,ncol=n)
##   return(rbind(cbind(zero,t(A)),cbind(A,zero)))
## }

## d2w_dir <- function(G,H,gamma){
##   return(d2w(undir(G),undir(H),gamma))
## }

## hamming_as_edit_dir <- function(G,H){
##   ## for weighted networks, weights must be in [0,1]
##   n=dim(G)[1]
##   return(sum(abs(undir(G)-undir(H)))/(2*n*(n-1)))
## }


him_dir <- function(G,H){
  n<-dim(G)[1]
  ipsen <- d2w_dir(G,H,optimal_gamma_dir(n))
  edit <- hamming_as_edit_dir(G,H)
  glocal_dist <- sqrt(edit**2/2+ipsen**2/2)
  return(c("H"=edit,"IM"=ipsen,"HIM"=glocal_dist))
}





