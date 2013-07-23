\name{hamming}
\alias{hamming}
\alias{hamming.default}
\alias{hamming.adj2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hamming distance
%%  ~~function to do ... ~~
}
\description{Compute the hamming distance
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hamming(object, ...)
\method{hamming}{default}(object,\dots)
\method{hamming}{adj2}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class "adj2"
}
  \item{\dots}{
further argument to be pass
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a value
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
standardGeneric("hamming"), generic = structure("hamming", package = "nettools"), package = "nettools", group = list(), valueClass = character(0), signature = "object", default = structure(function (object, 
    ...) 
UseMethod("hamming"), target = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), defined = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), generic = structure("hamming", package = "nettools"), class = structure("derivedDefaultMethod", package = "methods")), skeleton = structure(function (object, 
    ...) 
UseMethod("hamming"), target = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), defined = structure("ANY", class = structure("signature", package = "methods"), .Names = "object", package = "methods"), generic = structure("hamming", package = "nettools"), class = structure("derivedDefaultMethod", package = "methods"))(object, 
    ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
