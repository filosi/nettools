\encoding{UTF-8}
\name{mat2adj}
\alias{mat2adj}
\alias{mat2adj,data.frame-method}
\alias{mat2adj,matrix-method}
\alias{mat2adj.matrix}
\alias{mat2adj.data.frame}
\alias{mat2adj.default}
% \alias{mat2adj,mat2adj.data.frame}
% \alias{mat2adj,mat2adj.matrix}
%% \alias{mat2adj-methods}
%% \alias{mat2adj.data.frame,data.frame-method}
%% alias{mat2adj.matrix,matrix-method}

\title{From Data Matrix to Adjacency Matrix}

\description{\code{mat2adj} is a high level function providing
  different network inference methods. The function takes in input a data
  matrix \emph{N by P}, with N samples on the rows and P variables on the
  columns. The adjacency matrix \emph{P by P} will be computed with the
  specified method, using N samples to infer the interactions between
  the variables.
}

\usage{
mat2adj(x,\dots)

\method{mat2adj}{default}(x, \dots)
%% \S4method{mat2adj.default}{generic}(x, \dots)

\method{mat2adj}{data.frame}(x, \dots)
%% \S4method{mat2adj}{data.frame}(x, \dots)

\method{mat2adj}{matrix}(x,method='cor', FDR=1e-3, P=6, measure=NULL,
alpha=0.6, C=15,DP=1, \dots)

%% \S4method{mat2adj}{matrix}(x,method='cor', FDR=1e-3, P=6,
%% measure=NULL, alpha=0.6, C=15,DP=1, \dots)
}
%% \section{Methods}{
%%   \describe{
%%     \item{\code{signature(object = "ANY")}}{Generic function}
%%     \item{\code{signature(object = "data.frame")}}{Function for
%%       data.frame}
%%     \item{\code{signature(object = "matrix")}}{Function for matrix}
%%   }
%% }
\arguments{
  \item{x}{ a matrix or data.frame of numerical values of N rows and P columns}
  \item{method}{ a character string indicating which
    method will be used for inferring a relationship between two variables. This must
    be (an abbreviation of) one of \code{"cor"} (default), \code{"WGCNA"}, \code{"WGCNAFDR"},
    \code{"bicor"}, \code{"bicorFDR"}, \code{"TOM"}, \code{"ARACNE"},
    \code{"CLR"}, \code{"MINE"}, \code{"MINEFDR"}, \code{"DTWMIC"}}
  \item{P}{6 (default), integer used as soft-thresholding power for network
    construction, used by the \code{"WGCNA"} and \code{"TOM"} methods.}
  \item{FDR}{1e-3 (default), a number which indicates the number of
    values generated to compute the NULL hypothesis. To be used for
    methods \code{"WGCNAFDR"}, \code{"MINEFDR"} and \code{"bicorFDR"}}
  \item{measure}{\code{"MIC"} (default), a valid string indicating
    the measure of the MINE suite to compute. One of \code{"MIC"},
    \code{"MCN"}, \code{"MEV"}, \code{"MAS"} or \code{"MICR2"}.}
  \item{alpha}{0.6 (default), the \code{alpha} argument to be
    passed to the function \code{mine}. See also
    \code{\link[minerva]{mine}}}
  \item{C}{15 (default), an integer value to be passed at the mine function
    main. Only for methods \code{"MINE"}  and \code{"MINEFDR"}.}
  \item{DP}{1 (default), only for method \code{"DTWMIC"}.}
  \item{\dots}{Additional arguments to be passed to the downstream
    functions. Normally the argument passed through \dots are processed
    by the functions which compute the inference. Not all parameters are
    used by all functions.}
}
\details{\code{mat2adj} function is a high-level function which includes
  different methods for network inference. In particular the function
  infer the relation between all the possible pairwaise comparison
  between columns in the dataset. If the input is a \code{data.frame},
  columns were first converted into a numerical matrix. Given a N by P
  numerical matrix, the relation between each PxP pairs of
  variables is inferred with the selected method.

  The \code{"FDR"} corrected methods are based on a permutation estimate
  of the null hypothesis. A total amount of 1/(\code{"FDR"})
  permutations are performed to asses the reliability of the inferred
  link; each link is set only if it
  is inferred in all the permutations and its weight is lower then the
  value on non permuted data. The default value for \code{FDR} is 1e-3.

  All the available methods are the following:
  \describe{
        \item{\code{cor}}{(default) computes the interaction using the
	  'Pearson' correlation coefficient. Different correlation methods, such
	  as \code{Spearman} could be passed to the function using \dots.}
	
	\item{\code{ARACNE}}{Algorithm for the Reconstruction of Gene
	  Regulatory Networks, see also package \pkg{minet}}
	\item{\code{CLR}}{Context Likelihood of Relatedness see also
	  package \pkg{minet}}
	\item{\code{WGCNA}}{WeiGhted Correlation Network Analsysis. It is
	  based on a correlation measure. For further details see
	  the documentation of \pkg{WGCNA} package. The method accept
	  parameter \code{P} which is set to 6 by default}
	\item{\code{bicor}}{Biweighted Correlation method. It uses a
	  biweighted correlation as described in \link[WGCNA]{bicor}
	  package}
	\item{\code{TOM}}{Topological Overlap Measure inference method. For
	  further details see the documentation of \pkg{WGCNA}
	  package. As for \code{WGCNA} the parameter \code{P} can be
	  set(6 by default).}
	\item{\code{MINE}}{Maximum Information-based Non-parametric
	  Exploration. This method uses the \pkg{minerva} implementation of
	  the original measure. For this methods different measures are
	  available. See \pkg{minerva} for further information.
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
	\item{\code{MINEFDR}}{This calls an FDR corrected version of the
	  standard MINE method.  See the description for the \code{MINE}
	  method. Parameter \code{FDR=1e-3} (default) can
	  be set.}
	\item{\code{bicorFDR}}{This calls an FDR corrected version of the
	  \code{bicor} method. See the description for the \code{bicor}.
	  Parameter \code{FDR=1e-3} (default) can
	  be set.}
	\item{\code{WGCNAFDR}}{This calls an FDR corrected version of the
	  \code{WGCNA} method. Parameter \code{P} cannot be set for this
	  method. Parameter \code{FDR=1e-3} (default) can
	  be set.}
	\item{\code{DTWMIC}}{This method uses Dynamic Time Warping
	  transformation coupled witht the MIC statistic from the MINE
	  family. See Details for further information. Additional
	  parameters can be set with this method:
	}
      }
    }
    
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
      C. Furlanello (2013). \emph{minerva and minepy: a C engine for the
	MINE suite and its R, Python and MATLAB wrappers, Bioinformatics}
      
      \url{http://mpba.fbk.eu/cmine}\cr

      M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello
      (2014)\emph{Stability Indicators in Network Reconstruction,
      PLOSONE}
      
      D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P.
      Turnbaugh, E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
      \emph{Detecting novel associations in large datasets}  Science\cr
      
      \url{http://www.exploredata.net}\cr
      
      (SOM: Supplementary Online Material at
      \url{http://www.sciencemag.org/content/suppl/2011/12/14/334.6062.1518.DC1})
    }
    
    \author{Michele Filosi\cr 
      Special thanks to:
      Samantha Riccadonna, Giuseppe Jurman, Davide Albanese and Cesare
      Furlanello
    }
    
    \examples{
\dontrun{
data(Spellman, package="minerva")
dim(Spellman)
A <- mat2adj(Spellman,method="cor")
dim(A)
}
}
