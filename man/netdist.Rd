\name{netdist}
\alias{netdist}
\alias{HIM}
\alias{ipsen}
\alias{Ipsen}
\alias{hamming}
\alias{Hamming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distances between network}
\description{
  Function that computes the distance between two adjacency matrix given as
  matrix or igraph objects.
}

%- maybe also 'usage' for other objects documented here.
\usage{
netdist(g, h, d = "HIM", ga = NULL, components=TRUE, ...)
}

\arguments{
  \item{g}{Adjacency matrix or igraph object}
  \item{h}{Adjacency matrix or igraph object}
  \item{d}{\code{HIM}
    (default), or character string containing a valid method
    between the follwing: \code{IM}, \code{ipsen}, \code{Ipsen}, 
    \code{IpsenMikhailov}, or \code{Ipsen-Mikhailov}, for Ipsen-Mikhailov distance
    and \code{H}, \code{hamming}, or \code{Hamming} for Hamming distance.}
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
      \item{n.nodes:}{parameter used only when method is \code{ipsen} or
    	\code{HIM}.  If the number of nodes of \code{g} and \code{h} is
    	bigger than \code{n.nodes} the computation of the distance is
    	parallelized on 2 cores, otherwise only one core will be used.
      }
    }
  }
}

\details{
  Distance between networks can be computed using three methods through
  the \code{netdist} function.
  The input can be either adjacency symmetric matrix or graph object
  from \pkg{igraph} package.  
  If the 'weight' attribute is given within the \code{graph} object it
  will be used as edge weight, otherwise binary representation of the network will be considered.
  In case edge weights are provided through the adjacency matrix, they
  should be scaled between 0-1. If some of the weights are outside the
  interval a scaling function is automatically applied.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{For method \code{"HIM"} a named vector with 3 values is returned. 
  Respectively, \code{H} for "hamming", \code{IM} for "Ipsen-Mikhailov"
  and \code{HIM} for "him" distance.
  The parameter \code{ga} is set to \code{NULL}.
  For method \code{ipsen} the "Ipsen-Mikhailov" distance is returned.
  For method \code{hamming} the "Hamming" distance is returned.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\references{
  
  G. Jurman, R. Visintainer, S. Riccadonna, M. Filosi, C. Furlanello
  \emph{A glocal distance for network comparison
  arXiv 2012, arXiv:1201.2931}\cr
  
  %% ~put references to the literature/web site here ~
}
\author{
  Michele Filosi
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
