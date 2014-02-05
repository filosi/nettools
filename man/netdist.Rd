\name{netdist}
\alias{netdist}
\alias{HIM}
\alias{ipsen}
\alias{Ipsen}
\alias{hamming}
\alias{Hamming}
\alias{netdist,matrix-method}
\alias{netdist,data.frame-method}
\alias{netdist,Matrix-method}
\alias{netdist,list-method}
\alias{netdist.matrix}
\alias{netdist.Matrix}
\alias{netdist.list}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distances between network}
\description{
  This function computes the distance between two adjacency matrix given
  as matrices or igraph objects.
}

%- maybe also 'usage' for other objects documented here.
\usage{
netdist(x,\dots)

\method{netdist}{matrix}(x, h, d = "HIM", ga = NULL, components=TRUE, ...)

\method{netdist}{Matrix}(x, h, d = "HIM", ga = NULL, components=TRUE, ...)

\method{netdist}{list}(x, d = "HIM", ga = NULL, components=TRUE, ...)

}
\arguments{
  \item{x}{Adjacency matrix, igraph object or list of adjacency
    matrices/igraph objects}
  \item{h}{Adjacency matrix or igraph object. Only when \code{x} is a matrix object.}
  \item{d}{\code{HIM} (default), character string containing a valid method.
    Accepted values are: \code{"IM"}, \code{"ipsen"}, \code{"Ipsen"}, 
    \code{"IpsenMikhailov"}, or \code{"Ipsen-Mikhailov"}, for Ipsen-Mikhailov distance
    and \code{"H"}, \code{"hamming"}, or \code{"Hamming"} for Hamming distance.}
  \item{ga}{\code{NULL} (default), a numeric value for the gamma parameter used
    only if method \code{ipsen} is called.  If \code{ga=NULL} an optimal
    value based on the number of nodes in the adjacency matrix is
    computed.}
  \item{components}{\code{TRUE} (default), to obtain the components 
    (Hamming and Ipsen-Mikhailov) for \code{method=HIM}; 
    otherwise the function returns only the HIM value.
    If \code{method} is set to \code{hamming} or \code{ipsen} it will be ignored.
  }
  \item{\dots}{Additional arguments to be passed to the
    downstream functions. Normally the argument passed through \dots are
    processed by the functions which compute the distance. Not all
    parameters are used by all functions.
    \describe{
      \item{n.cores}{parameter used only when method is \code{ipsen} or
    	\code{HIM}. Very useful if the number of nodes in the adjacency
	matrix is high and when method for object \code{list} is called.
	%% Insert the following code in the details
	If method \code{list} is called, each core take care of the
	computation of eigenvalues for a matrix in the input list. If
	\code{n.cores} is bigger than the number of available cores on the
	specific machine, than \code{n.cores} is set to the number of
	available cores - 1.}
      \item{verbose}{\code{FALSE} (default), print information about the
	state of the computation}
    }
  }
}

\details{
  \code{netdist} is a high level function which provide an interface to
  different network distance methods. Distance can be computed for
  binary or weighted networks either directed and undirected.
  Each network should have the same number of nodes (the matrices
  provided as input should have the same dimensions), otherwise no
  computation will be performed and an errore message will be output.
  
  Parallel computation is provided automatically through the functions
  available in the \pkg{parallel} package included by default from \R
  2.15, only for methods \code{"ipsen"} and \code{"HIM"}.
  The computation can be automatically parallelized on a multi-cpu
  computer using the parameter \code{n.cores}.
  If \code{n.cores} is set to \code{NULL} it searches for multiple cores
  on the machine, if it finds more than 2 cores the computation is
  parallelized using the maximum number of cores available - 1.
  If \code{n.cores} is a numeric value it checks if the number of cores provided
  is suitable for the machine where \R is running. If the number of
  cores requested is bigger than the available cores \code{n.cores} will
  be set to \code{NULL}, otherwise the computation will be equally
  distributed on the number of cores provided.
  
  The input can be either an adjacency matrix, a graph object
  from \pkg{igraph} package or a list of adjacency matrix or graph
  objects.
  If the 'weight' attribute is given within the \code{graph} object it
  will be used as edge weight, otherwise binary representation of the network will be considered.
  In case edge weights are provided through the adjacency matrix, they
  should be scaled between 0-1. If some of the weights are outside the
  interval a scaling function is automatically applied.
  
  The parameter \code{ga} is used only in the \emph{Ipsen} and
  \emph{HIM} distance. We suggest not to change it and use the automatic
  computation of the parameter based on the number of nodes in the
  network.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{For method \code{HIM} and if parameter \code{components} is set
  to \code{TRUE} a named vector with 3 values is returned, respectively
  the \emph{Hamming} distance, the \emph{Ipsen-Mikhailov} distance and
  the \emph{HIM} distance.
  For method \code{ipsen} the "Ipsen-Mikhailov" distance is returned.
  For method \code{hamming} the "Hamming" distance is returned.
  If a \code{list} object is provided as input a distance matrix is
  returned.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\references{
  G. Jurman, R. Visintainer, M. Filosi, S. Riccadonna, C. Furlanello
  \emph{The HIM glocal metric and kernel for network comparison and classification arXiv 2013, arXiv:1201.2931v3}\cr
%% ~put references to the literature/web site here ~
}
\author{
  M. Filosi, S. Riccadonna
}
%% \note{
%%   %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%%   %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
a <- matrix(rnorm(2000), nrow=200)
b <- matrix(rnorm(2000), nrow=200)

aa <- mat2adj(a, method="cor")
bb <- mat2adj(b, method="cor")

dd <- netdist(aa,bb, d="HIM")
print(dd)

## Change the gamma parameter
dd <- netdist(aa,bb, d="HIM", ga=0.5)
print(dd)

alist <- list(aa,bb)
netdist(alist, d="HIM")

## Use multiple cores
\dontrun{
dd <- netdist(aa,bb,n.cores=4)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dist}
\keyword{network}% __ONLY ONE__ keyword per line
