\name{netSI}
\alias{netSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network Stability Indicators
}
\description{
Family NetSI of stability indicators for network inference methods. The
function computes a suite of 4 stability indicators of data \code{x}
based on network distances.
}
\usage{
netSI(x,indicator="all", d='HIM', adj.method='cor',
method="montecarlo", k=3, h=20, n.cores, save=FALSE, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numerical matrix or data.frame to be used for network inference.}
  \item{indicator}{\code{all} (default) character string. Indicating all or
    one of the four indicators available. This must be on of the strings
    \code{"S"},\code{"SI"}, \code{"Sw"} or \code{"Sd"}.}
  \item{d}{\code{HIM} (default), character string indicating the
    distance to be used for the stability computation in case of
    indicators S and SI.}
  \item{adj.method}{\code{"cor"} (default), optional character string
    giving one of the methods available for network inference. For
    details see \code{\link{mat2adj}}.
  }
  \item{method}{\code{"montecarlo"} (default), a string identifying the
    cross validation schema to be used for the stability computation.
    This must be one of the strings \code{montecarlo} (default),
    \code{LOO} for Leave One Out or \code{kCV} for k-fold Cross
    Validation.}
  \item{k}{3 (default) numeric value indicating the number of subset the
    data will be splitted in. In case \code{method} is set to
    \code{montecarlo} the 1-1/\code{k} samples are taken for the computation.
    In case \code{cvlab} is set to \code{loo} it should be set to 1,
    otherwise it will be automatically set to 1.
    In case \code{cvlab} is set to \code{kCV} \code{h} iteration are
    performed, each dividing the data into randomly chosen \code{k} groups}
  \item{h}{ 20 (default), numeric value indicates the number of
    iterations the \code{montecarlo} and \code{kCV} procedure should be repeated(cross validations).}
  \item{n.cores}{\code{NULL} (default), optional numeric value giving
    the number of cores to be used for parallel computation. The default
    behaviour try to asses the maximum number of available cores and set
    a parallel computation using the maximum number of cores - 1 for the
    computation.}
  \item{save}{\code{FALSE} (default), logical. Should adjacency matrices
    and indicator's raw value from resampling scheme be saved?} 
  \item{verbose}{\code{TRUE} (default), logical. Print all the
    information about the computation. If \code{FALSE} messages will be suppressed.}
  \item{\dots}{Additional arguments to be passed to the downstream
    functions. Not all parameters are used by all functions. See
    function \code{\link{netdist}} and \code{\link{mat2adj}}.
    \describe{
      \item{sseed:}{random number generation seed}
    }
  }
}

\details{
  This function computes a suite of 4 stability indicators. It is based
  on network inference method function \code{\link{mat2adj}} and network
  distance function \code{\link{netdist}}.
  Parameters passed through \dots will be used by the low-level
  functions \code{\link{mat2adj}} and \code{\link{netdist}}.

  Indicators:
  \itemize{
    \item{\code{S}: is the global stability indicator. It asses the
      pertubations of the network given by the resampling: in particular it
      account for the distance between the network inferred using all
      samples and the network inferred on a subsampling.}
    \item{\code{SI}: is the local stability indicator. It provides a
      measure of the perturbations between different subsamplings. It
      computes the pairwise distances between networks inferred on different subsamplings.}
    \item{\code{Sw}: is the edge stability indicator. It asses the
      stability of the edges inferred in different subsamplings. In case
      of a binary network it take into account only the presence/absence of
      the link, if a "weight" attribute is provided it asses the
      perturbation of the edge weight given the presence/absence on the
      given subsampling.}
    \item{\code{Sd}: is the degree stability indicator. It asses the
      variations of the node degree given a subsampling perturbation on
      the data.}
  }
 
  
  Parallel computation is provided automatically through the functions
  available in the \pkg{parallel} package included by default from \R
  2.15. Support for parallel computation is available only on multi-core
  machines not on high performace computing facility.
  If \code{n.cores} is set to \code{NULL} it searches for multiple cores
  on the machine, if it finds more than 2 cores the computation is
  parallelized using the maximum number of cores available - 1.
  If \code{n.cores} is a numeric value it checks if the number of cores provided
  is suitable for the machine where \R is running. If the number of
  cores requested is bigger than the available cores \code{n.cores} will
  be set to \code{NULL}, otherwise the computation will be equally
  distributed on the number of cores provided.
  In case \code{indicator} is set to \code{all} or \code{SI}, parameter
  \code{n.cores} will pass to the low level function which compute
  the mutual distances, \link{netdist}.
  %%  ~~ If. necessary, more details than the description above ~~
}
\value{
  A list with the followind elements:
  If \code{save=TRUE}:
  \item{call}{character. Function call}
  \item{ADJlist}{A list with the adjacency matrices computed in the
    resampling scheme.}
  \item{S}{The value of the global indicator "S"}
  \item{SI}{The value of the global internal indicator "SI"}
  \item{Sw}{The value of the edge weight stability indicator "Sw"}
  \item{Sd}{The value of the node stability indicator "Sd"}
  \item{S_boot}{A numerical vector with the value of the stability indicator S
    for each resampling.}
  \item{SI_boot}{A numerical vector with the value of the stability indicator SI
    for each combination of pairs within all the resampling.}
  \item{Sw_boot}{A numerical matrix. Each row contains the value of the
    weight of a particular edge on different resampling splits.}
  \item{Sd_boot}{A numerical matrix number_of_resampling X number_of_nodes.
    Each column contain the degree of a node on different resampling splits.}

  If \code{save=FALSE} elements \code{call}, \code{ADJlist}, \code{ADJ},
  \code{S_boot}, \code{SI_boot}, \code{Sw_boot} and \code{Sd_boot} will
  not be provided. For \code{S} and \code{SI} the mean over the
  resampling splits is computed, instead for the element \code{Sw} and
  \code{Sd} the range/mean is computed over resamplings.
}
\references{
  M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello
  (2014)\emph{Stability Indicators in Network Reconstruction,
    PLOSONE}\cr
  
%% ~put references to the literature/web site here ~
}
\author{
  Michele Filosi, Samantha Riccadonna
  %%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{netdist}}, \code{\link{mat2adj}}, \code{\link{WGCNA}}, \code{\link{minet}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
myrawdata <- matrix(rnorm(500),ncol=20)

sstab <- netSI(myrawdata, d="HIM", n.cores=1, save=FALSE)
print(sstab$S)
print(sstab$SI)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stability }
\keyword{ network }% __ONLY ONE__ keyword per line
\keyword{ reproducibility }
