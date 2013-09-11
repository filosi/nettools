\name{netSI}
\alias{netSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network Stability Indicators
}
\description{
Family NetSI of stability indicators for network inference methods
}
\usage{
netSI(d,indicator="all", dist='HIM', adj.method='cor',
method="montecarlo", k=3, h=20, n.cores=NULL, save=TRUE, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{dataset used for network inference and stability indicators
    computation}
  \item{indicator}{\code{all} (default) character string denoting all or
    one of the four indicators available. This must be on of the strings
    \code{"S"},\code{"SI"}, \code{"Sw"} or \code{"Sd"}}
  \item{dist}{\code{HIM} (default), character string indicating the
    distance to be used for the stability computation in case of
    indicators S and SI.}
  \item{adj.method}{\code{"cor"} (default), optional character string
    giving one of the methods available for network inference. For
    details see \code{\link{mat2adj}}.
  }
  \item{method}{\code{"montecarlo"} (default), a string identifying the
    cross validation schema to be used for the stability computation.
    This must be one of the strings \code{"montecarlo"} (default),
    \code{"LOO"} for Leave One Out or \code{"kCV"} for k-fold Cross
    Validation.}
  \item{k}{3 (default) numeric value indicating the number of subset the
    data will be splitted in.  In case method=="montecarlo" the 1-1/cv
    samples are taken for the computation In case cvlab=="loo" it should
    be set to 1 In case cvlab=="kcv" cv iteration are performed, each
    dividing the data into 1-1/cv }
  \item{h}{ 20 (default), numeric value indicates the number of iterations (cross validations).}
  \item{n.cores}{\code{NULL} (default), optional numeric value giving
    the number of cores to be used for computation. The default
    behaviour try to asses the maximum number of available cores and set
    a parallel computation using the maximum number of cores - 1 for the
    computation.}
  \item{save}{\code{TRUE} (default), logical. Should adjacency matrices
    from resampling scheme be saved?} 
  \item{verbose}{\code{TRUE} (default), logical. Print all the
    information about the computation. If \code{FALSE} messages will be suppressed.}
  \item{\dots}{Additional arguments to be passed to the downstream
    functions. Not all parameters are used by all functions. See
    function \code{\link{netdist}} and \code{\link{mat2adj}}.}
}
  
  \details{
  For different network stability indicator are computed. This function
  is based on a resampling method functions \code{\link{mat2adj}} and \code{\link{netdist}}

  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with the followind elements:
  If \code{save=TRUE}:
  \item{call}{string. Function call}
  \item{ADJlist}{A list with the adjacency matrices computed in the
    resampling scheme.}
  \item{S}{A numerical vector with the value of the stability indicator S
    for each resampling.}
  \item{SI}{A numerical vector with the value of the stability indicator SI
    for each combination of pairs within all the resampling.}
  \item{Sw}{A numerical matrix. Each row contains the value of the
    weight of a particular edge on different resampling splits.}
  \item{Sd}{A numerical matrix number_of_resampling X number_of_nodes.
    Each column contain the degree of a node on different resampling splits.}

  If \code{save=FALSE} elements \code{"call"} and \code{"ADJlist"} will
  not be provided. For the elements \code{S} and \code{SI} the mean over the
  resampling splits is computed, instead for the element \code{Sw} the
  mean is computed for each edge and for the element \code{Sd} the mean
  is computed for each node.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\references{
G. Jurman, M. Filosi, R. Visintainer, S. Riccadonna, C. Furlanello
\emph{Stability indicators in network reconstruction
arXiv 2013, arXiv:1209.1654}\cr
%% ~put references to the literature/web site here ~
}
\author{
  Michele Filosi, Samantha Riccadonna
  %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{netdist}}, \code{\link{mat2adj}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
