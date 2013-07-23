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


