\name{mpnetdist}
\alias{mpnetdist}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distances between Multiplex-networks}
\description{
  This function computes the Ipsen-Mikhailov distance between two multiplex networks.
}

%- maybe also 'usage' for other objects documented here.
\usage{
mpnetdist(x,y, d="HIM", ga=NULL, components=TRUE, \dots)
}

\arguments{
  \item{x}{A list object, see Description for more details.}
  \item{y}{A list object, see Description for more details.}
  \item{d}{\code{HIM} (default), character string containing a valid method.
    Accepted values are: \code{"IM"}, \code{"ipsen"}, \code{"Ipsen"}, 
    \code{"IpsenMikhailov"}, or \code{"Ipsen-Mikhailov"}, for Ipsen-Mikhailov distance
    and \code{"H"}, \code{"hamming"}, or \code{"Hamming"} for Hamming distance.}
  \item{ga}{\code{NULL} (default), a numeric value for the gamma
	parameter.  If \code{ga=NULL} an optimal value based on the number of
	nodes in the adjacency matrix is computed.}
  \item{components}{\code{TRUE} (default), to obtain the components 
    (Hamming and Ipsen-Mikhailov) for \code{method=HIM}; 
    otherwise the function returns only the HIM value.
    If \code{method} is set to \code{hamming} or \code{ipsen} it will be ignored.}
  \item{\dots}{Additional
	arguments to be passed to the downstream functions. Normally the
	argument passed through \dots are processed by the functions which
	compute the distance. Not all parameters are used by all functions.  }
}

\details{
  \code{mpnetdist} is a high level function which provides an interface to
  the Ipsen-Mikhailov distance for multiplex-network.
  The network can have 1 or more layers with the same number of nodes
  per layer.
  Multiplex network can be represented using a list object, with two
  element I (interlayer adjacency matrix) and L (intralayer adjacency
  matrix).
  Given a network with l layers and n nodes per layer the object I of
  the list should be a matrix of dimension \eqn{l \times l} with no self-loops
  hence the diagonal should be a vector of zeros.
  The second element of the list (L) is a 3d array. The array should
  contain the adjacency matrix for each layer and should be a \eqn{n \times
  n \times l} array.
  
  Parallel computation is provided automatically through the
  \pkg{parallel} package included by default from \R 2.15.
  The computation can be automatically parallelized on a multi-cpu
  computer using the parameter \code{n.cores}.
  
  We suggest not to change \code{ga} and use the automatic
  computation of the parameter based on the number of nodes in the
  network.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{The Ipsen-Mikhailov distance value between the input network.
}
\references{
  A. Sole-Ribalta, M. De Domenico, N. E. Kouvaris, A. Diaz-Guilera,
  S. gomez, A. Arenas \emph{Spectral Properties of the Laplacian of
	multiplex networks, 2013, arXiv:1307.2090v1}\cr
  
  G. Jurman, R. Visintainer, M. Filosi, S. Riccadonna, C. Furlanello
  \emph{The HIM glocal metric and kernel for network comparison and classification arXiv 2013, arXiv:1201.2931v3}\cr
%% ~put references to the literature/web site here ~
}
\author{
  M. Filosi
}
%% \note{
%%   %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%%   %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{


## Create toy dataset
a <- matrix(rnorm(1000),ncol=100)
b <- matrix(rnorm(1000),ncol=100)
cc <- matrix(rnorm(1000),ncol=100)

## Compute adjacency matrix!
## For different methods check the mat2adj help page!
aadj <- mat2adj(a, method="MINE", n.cores=1)
badj <- mat2adj(b, method="MINE", n.cores=1)
cadj <- mat2adj(cc, method="MINE", n.cores=1)

myarr <- array(data=NA, c(100,100,3))
myarr[,,1] <- aadj
myarr[,,2] <- badj
myarr[,,3] <- cadj

iadj <- matrix(c(0,1,0,1,0,1,0,1,0), ncol=3, byrow=TRUE)

net1 <- list(I=iadj, L=myarr)
net2 <- net1

dd <- mpnetdist(net1, net2, n.cores=1)
## IM
## 0

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dist}
\keyword{network}% __ONLY ONE__ keyword per line
