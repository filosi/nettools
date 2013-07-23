\name{ipsen}
\alias{ipsen}
\alias{ipsen.default}
\alias{ipsen.adj2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ipsen-Mikhailov distance
}
\description{
It computes the Ipsen-Mikhailov for an object of class \code{"adj2"}
}
\usage{
ipsen(object, ...)
\method{ipsen}{default}(object,\dots)
\method{ipsen}{list}(object,\dots,gamma=0.08)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class "list".
    The list should have at least 2 elements,
    each of them is a laplacian matrix.
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{further parameters to be passed to the function.
%%     ~~Describe \code{\dots} here~~
}
\item{gamma}{Gamma parameter, if gamma="NA" an optimal gamma is computed
according to the number of nodes of the network}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Something
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
structure(function (object, ...) 
standardGeneric("ipsen"), generic = structure("ipsen", package = "nettools"), package = "nettools", group = list(), valueClass = character(0), signature = "object", default = structure(function (object, 
    ...) 
UseMethod("ipsen"), target = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), defined = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), generic = structure("ipsen", package = "nettools"), class = structure("derivedDefaultMethod", package = "methods")), skeleton = structure(function (object, 
    ...) 
UseMethod("ipsen"), target = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), defined = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), generic = structure("ipsen", package = "nettools"), class = structure("derivedDefaultMethod", package = "methods"))(object, 
    ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
