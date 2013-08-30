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
  The function computes the distance between two adjacency matrix given as
  matrix or igraph objects
}

%- maybe also 'usage' for other objects documented here.
\usage{
netdist(g, h, method = "HIM", g = NULL)
}

\arguments{
  \item{g}{Adjacency matrix or igraph object}
  \item{h}{Adjacency matrix or igraph object}
  \item{method}{\code{HIM} (default), or character string containing a
    valid method or unique abbreviation between \code{ipsen} and
    \code{hamming}.
  }
  \item{g}{NULL (default), a numeric value for the gamma parameter,
    only if method \code{"ipsen"} is called, otherwise it will be
    ignored. If \code{"g=NULL"} an optimal value based on the number
    of columns in the adjacency matrix is computed. 
    See function \code{"HIM"} for further details.
  }
}
\details{
  Different methods for network distance can be called using
  \code{netdist} function.
  The input can be either adjacency symmetric matrix or undirected graph
  object from \pkg{igraph} package.
  If the 'weight' attribute is given within the \code{graph} object it
  will be used as edge weight, otherwise binary representation of the
  network will be considered.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{For method \code{"HIM"} a named vector with 3 values is returned. 
Respectively, "hamming", "Ipsen-Mikhailov" and "him" distance. 
The parameter \code{"g"} is set to "NA".
For method \code{"ipsen"} the "Ipsen-Mikhailov" distance is returned.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
