\name{netdist}
\alias{netdist}
\alias{HIM}
\alias{ipsen}
\alias{hamming}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distance between network}
\description{
  The function compute the distance between two adjacency matrix given as
  matrix or igraph objects
}

%- maybe also 'usage' for other objects documented here.
\usage{
netdist(g, h, method = "HIM", gamma = NULL)
}

\arguments{
  \item{g}{Adjacency matrix or igraph object}
  \item{h}{Adjacency matrix or igraph object}
  \item{method}{\code{HIM} (default), or character string containing a
    valid method or unique abbreviation between \code{ipsen} and
    \code{hamming}.
  }
  \item{gamma}{NULL (default), a numeric value for the gamma parameter,
    only if method \code{"ipsen"} is called, otherwise it will be
    ignored. If \code{"gamma=NULL"} an optimal value based on the number
    of columns in the adjacency matrix is computed. 
    See function \code{"HIM"} for further details.
  }
}
\details{
  Different methods for network distance can be callde using
  \code{netdist} function.
  The input can be either adjacency symmetric matrix or undirected graph
  object from \pkg{igraph} package.
  If the 'weight' attribute is given within the \pkg{graph} object it
  will be use ad edge weight, otherwise binary representation of the
  network will be considered.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{For method \code{"HIM"} a named vector with 3 values is returned. Respectively, "hamming", "Ipsen-Michailov" and "him" distance. The parameter \code{"gamma"} is set to "NA".
For method \code{"ipsen"} the "Ipsen-Michailov" distance is returned.
For method \code{"hamming"} the "Hamming" distance is returned.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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

## The function is currently defined as
function (g, h, method = "HIM", gamma = 0.08) 
{
    METHODS <- c("HIM", "ipsen", "hamming")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    g1 <- g2adj(g)
    g2 <- g2adj(h)
    myadj <- new("adj2", G1 = g1, G2 = g2, N1 = ncol(g1), N2 = ncol(g2), 
        L1 = Lap(g1), L2 = Lap(g2), method = METHODS[method])
    if (myadj@method == "HIM") 
        dd <- him(myadj)
    if (myadj@method == "ipsen") 
        dd <- ipsen(myadj, gamma = gamma)
    if (myadj@method == "hamming") 
        dd <- hamming(myadj)
    return(dd)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
