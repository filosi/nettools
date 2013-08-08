\encoding{UTF-8}
\name{mat2adj}
\alias{mat2adj}
\alias{mat2adj.default}
\alias{mat2adj.data.frame}
\alias{mat2adj.matrix}
\title{Data Matrix to Adjacency Matrix}
\usage{
mat2adj(x,\dots)

\method{mat2adj}{default}(x, \dots)

\method{mat2adj}{data.frame}(x, \dots)

\method{mat2adj}{matrix}(x,method='cor', FDR=1e-3, P=6, measure=NULL, alpha=0.6, C=15,DP=1, \dots)
}

\arguments{
  \item{x}{A matrix or data.frame of numerical values of N rows and
  P columns}
  \item{method}{\code{cor} (default), or a character string indicating which
    method to use for inferring a relation between two variables. This must
    be (an abbreviation of) one of the strings \code{WGCNA}, \code{WGCNAFDR},
    \code{bicor}, \code{bicorFDR}, \code{TOM}, \code{ARACNE},
    \code{CLR}, \code{MINE}, \code{MINEFDR}, \code{DTWMIC}}
  \item{P}{6 (default), integer used as soft-thresholding power for netwoek
    construction, used by the \code{WGCNA} and \code{TOM} methods.}
  \item{FDR}{1e-3 (default), a number which indicates the number of
    values generated to compute the NULL hypothesis. To be used for
    methods \code{WGCNAFDR}, \code{MINEFDR} and \code{bicorFDR}}
  \item{measure}{\code{MIC} (default), a valid string indicating
    the measure of the MINE suite to compute. One of \code{MIC},
    \code{MCN}, \code{MEV}, \code{MAS} or \code{MICR2}.}
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
  different network inference methods. The function takes in input a data
  matrix N by P, with N samples on the rows and P variables on the
  columns. The adjacency matrix P by P will be computed with the
  specified method, using N samples to infer the interactions between
  the variables.
}

\details{\code{mat2adj} function is a high-level function which includes
  different methods for inferring the relation between variables in a
  dataset.  If the input is a \code{data.frame} the columns will be
  converted into numerical values.  Given a N by P numerical matrix as
  input, the relation between each PxP pairs of variables is inferred.
  The available methods are the following:
  \describe{
    \item{\code{cor}}{(default) computes the interaction using the
      'Pearson' correlation coefficient. Different correlation methods, such
  as \code{Spearman} could be passed to the function using \dots.}
    
  }
  
  \code{ARACNE}
  \code{CLR}
  The method \code{TOM} (Topological Overlap Matrix) calculates the
  topological overlap matrix from a given adjacency matrix. The
  Adjacency matrix is computed with the method \code{cor} using the
  correlation coefficient and the soft thresholding P passed by the user
  (The Defaults are Pearson and P = 6)\cr
  \code{MINE} methods use a non-linear mutual
  information coefficient to asses statistical dependency between
  variables. This method is based on the \pkg{minerva} package
  implementation.\cr
  
  For some methods an \code{FDR} correction is available, based the
  estimation of the null hypothesis through permutations. Just the links
  with a P-value lower than \code{FDR} are inferred. 1/(\code{FDR})
  permutations are performed and each link is considered real just if it
  is inferred in all the permutations. The default value is 1e-3
  
  To clarify the main MINE family statistics let \eqn{D={(x,y)}} be the
  set of \emph{n} ordered pairs of elements of \code{x} and
  \code{y}. The data space is partitioned in an \emph{X-by-Y} grid,
  grouping the \emph{x} and \emph{y} values in \emph{X} and \emph{Y}
  bins respectively.\cr
  The value of \code{alpha} (default 0.6) has been
  empirically chosen by the authors of the original paper.\eqn{alpha} is
  the exponent of the search-grid size \eqn{B(n)=n^{\alpha}}. It is
  worthwhile noting that \code{alpha} and \code{C} are defined to obtain
  an heuristic approximation in a reasonable amount of time. In case of
  small sample size (\emph{n}) it is preferable to increase \code{alpha}
  to 1 to obtain a solution closer to the theoretical one.\cr
  \code{C} determines the number of starting point of the
  X-by-Y search-grid. When trying to partition the x-axis into X
  columns, the algorithm will start with at most C x X clumps. Default
  value is 15.\cr
  The \strong{Maximal Information Coefficient (MIC)} is defined as
  \deqn{\textrm{MIC}(D)=\max_{XY<B(n)} M(D)_{X,Y} = \max_{XY<B(n)}
  \frac{I^*(D,X,Y)}{log(\min{X,Y})},}{MIC(D)=max_{XY<B(n)}
  M(D)_{X,Y}=max_{XY<B(n)} I*(D,X,Y)/log(min(X,Y)),} where
  \eqn{B(n)=n^{\alpha}} is the search-grid size,
  \eqn{I^*(D,X,Y)}{I*(D,X,Y)} is the maximum mutual information over all
  grids \emph{X-by-Y}, of the distribution induced by D on a grid having
  \emph{X} and \emph{Y} bins (where the probability mass on a cell of
  the grid is the fraction of points of D falling in that cell).  The
  other statistics of the MINE family are derived from the mutual
  information matrix achieved by an \emph{X-by-Y} grid on D.
  The \strong{Maximum Asymmetry Score (MAS)} is defined as
  \deqn{\textrm{MAS}(D) = \max_{XY<B(n)} |M(D)_{X,Y} -
    M(D)_{Y,X}|.}{MAS(D) = max_{XY<B(n)} |M(D)_{X,Y} - M(D)_{Y,X}|.}
  The \strong{Maximum Edge Value (MEV)} is defined as
  \deqn{\textrm{MEV}(D) = \max_{XY<B(n)} \{M(D)_{X,Y}:
    X=2~or~Y=2\}.}{MEV(D) = max_{XY<B(n)} {M(D)_{X,Y}: X=2 or Y=2}.}
  The \strong{Minimum Cell Number (MCN)} is defined as
  \deqn{\textrm{MCN}(D,\epsilon) = \min_{XY<B(n)} \{\log(XY): M(D)_{X,Y}
    \geq (1-\epsilon)MIC(D)\}.}{MCN(D,\epsilon) = min_{XY<B(n)} {log(XY):
      M(D)_{X,Y} >= (1-\epsilon)MIC(D)}.}  More details are provided in the
  supplementary material (SOM) of the original paper.}
  
\value{A P by P symmetric adjacency matrix with the diagonal set to
  0. Self loop and direction of the edges are not taking into
  account. The values range in [0, 1].
}

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

  Jeremiah J Faith, Boris Hayete, Joshua T Thaden, Ilaria Mogno, Jamey
  Wierzbowski, Guillaume Cottarel, Simon Kasif, James J Collins, Timothy
  S Gardner. \emph{Large-Scale Mapping and Validation of Escherichia
  coli Transcriptional Regulation from a Compendium of Expression
  Profiles}
  
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
