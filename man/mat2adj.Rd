\encoding{UTF-8}
\name{mat2adj}
\alias{mat2adj}
\alias{mat2adj.default}
\alias{mat2adj.data.frame}
\title{Data Matrix to Adjacency Matrix}
\usage{
mat2adj(x,\dots)

\method{mat2adj}{default}(x,method='cor', FDR=1e-3, P=6, measure='MIC', alpha=0.6, C=15,DP=1, \dots)

\method{mat2adj}{data.frame}(x,method='cor', FDR=1e-3, P=6, measure='MIC', alpha=0.6, C=15,DP=1, \dots)
}

\arguments{
  \item{x}{A matrix or data.frame of numerical values of N rows and
  P columns}
  \item{method}{\code{cor} (default), or a character string indicating which
    method to use for inferring a relation between two variables. This must
    be (an abbreviation of) one of the strings \code{WGCNA}, \code{WGCNAFDR},
    \code{bicor}, \code{bicorFDR}, \code{TOM}, \code{ARACNE},
    \code{CLR}, \code{MINE}, \code{MINEFDR}, \code{DTWMIC}}
  \item{P}{6 (default), integer which could be interpreted as the
    power degree of the correlation used by the call of a \code{WGCNA}
    method.}
  \item{FDR}{1e-3 (default), a number which indicate the number of
    values generated to compute the NULL hypothesis. To be use for
    methods \code{WGCNAFDR} and \code{bicorFDR}}
  \item{measure}{\code{MIC} (default), a valid string indicating
    the measure of the MINE suite to compute. One of \code{MIC},
    \code{MCN}, \code{MEV}, \code{MAS} or \code{MIC-R2}.}
  \item{alpha}{0.6 (default), the \code{alpha} argument to be
    passed to the function \code{mine}. See also
    \code{\link[minerva]{mine}}}
  \item{C}{15 (default), an integer value to be passed at the mine function
    main. Only for methods \code{MINE}  and \code{MINEFDR}.}
  \item{DP}{1 (default), only for method \code{DTWMIC}.}
  \item{\dots}{Additional argument to be pass to the downstream
    functions. Normally the argument passed through \dots are processed
    by the functions which compute the inference. Not all parameters are
    used by all functions.
  }
}
\description{\code{mat2adj} is a high level function which collects
  different network inference methods. The function take in input a data
  matrix N by P, with the samples on the rows and the variables on the
  columns.  An adjacency matrix P by P will be computed using N samples to
  infer the interactions.
}

\details{\code{mat2adj} function is a high-level function which includes
  different methods for inferring the relation between variables in a
  dataset.
  If the input is a \code{data.frame} the columns will be casted into
  numerical values.
  Given a N by P numerical matrix as input, the function infers the
  relation between each PxP pairs of variables.
  The default method \code{cor}, compute the interaction using the
  'Pearson' correlation coefficient. Different correlation method, such
  as \code{Spearman} could be passed to the function using \dots.\cr
  \code{MINE} methods use a non-linear mutual information coefficient to
  asses statistical dependency between variables. This method is based
  on the \pkg{minerva} package implementation.
  For some methods a \code{FDR} correction is available, based the
  estimation of the null hypothesis through permutations.
}

\value{A P by P symmetric adjacency matrix with the diagonal set to
  0. The values are in range [0 - 1]}

\seealso{
  \code{\link{WGCNA}}, \code{\link{minerva}}, \code{\link{minet}}, \code{\link{cor}}
}

\references{
  P. Langfelder, S. Horvath (2008) \emph{WGCNA: an R package for
    weighted correlation network analysis. BMC Bioinformatics 2008,
    9:559}\cr

  P. E. Meyer, F. Lafitte, G. Bontempi (2008).
  \emph{MINET: An open source R/Bioconductor Package for Mutual Information
    based Network Inference. BMC Bioinformatics}

  \url{http://www.biomedcentral.com/1471-2105/9/461}\cr

  D. Albanese, M.Filosi, R. Visintainer, S. Riccadonna, G. Jurman,
  C. Furlanello (2012). \emph{minerva and minepy: a C engine for the
    MINE suite and its R, Python and MATLAB wrappers, Bioinformatics 2012}

  \url{http://mpba.fbk.eu/cmine}\cr

  D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P.
  Turnbaugh, E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
  \emph{Detecting novel associations in large datasets}  Science\cr

  \url{http://www.exploredata.net}\cr
  
  (SOM: Supplementary Online Material at
  \url{http://www.sciencemag.org/content/suppl/2011/12/14/334.6062.1518.DC1})
}

\author{Michele Filosi\cr Roberto Visintainer\cr
  Special thanks to:
Samantha Riccadonna, Giuseppe Jurman, Davide Albanese and Cesare
Furlanello
}

\examples{
\dontrun{
data(Spellman, package="minerva")
dim(Spellman)
A <- mat2adj(Spellman)
dim(A)
}
}
