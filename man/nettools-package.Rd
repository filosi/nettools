\name{nettools-package}
\alias{nettools-package}
\alias{nettools}
\alias{NetTools}
\alias{netTools}
\docType{package}
\title{ A Network Comparison Framework}
\description{ This package provides a collection of network inference
  methods for co-expression networks, quantitative network distances and
  a novel framework for network stability analysis.
  }

\details{
  Summary:
  \tabular{ll}{
    Package: \tab nettools\cr
    Type: \tab Package\cr
    Version: \tab 1.0.4\cr
    Date: \tab \cr
    Depends: R (>= 2.14.1), methods
    Imports: WGCNA, minet, rootSolve, dtw, Matrix, parallel, minerva, Rcpp,
    combinat, igraph, c3net, bc3net, R.utils
    License: CC BY-NC-SA 4.0
  }
}

\author{
  Michele Filosi [aut, cre],
  Roberto Visintainer [aut],
  Samantha Riccadonna [aut],
  Giuseppe Jurman [ctb],
  Cesare Furlanello [ctb]\cr
  Maintainer: Michele Filosi <filosi@fbk.eu>
}
\references{
    M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello
    (2014) \emph{Stability Indicators in Network Reconstruction,
    PLOSONE}\cr

    A. Sole-Ribalta, M. De Domenico, N. E. Kouvaris, A. Diaz-Guilera,
    S. gomez, A. Arenas \emph{Spectral Properties of the Laplacian of
	multiplex networks, 2013, arXiv:1307.2090v1}\cr

    P. Langfelder, S. Horvath (2008) \emph{WGCNA: an R package for
    weighted correlation network analysis. BMC Bioinformatics 2008,
    9:559}\cr
  
    P. E. Meyer, F. Lafitte, G. Bontempi (2008).
    \emph{MINET: An open source R/Bioconductor Package for Mutual Information
    based Network Inference. BMC Bioinformatics}.
    \url{http://www.biomedcentral.com/1471-2105/9/461}\cr

    D. Albanese, M.Filosi, R. Visintainer, S. Riccadonna, G. Jurman,
    C. Furlanello (2012). \emph{minerva and minepy: a C engine for the
    MINE suite and its R, Python and MATLAB wrappers, Bioinformatics 2012}.
    \url{http://minepy.readthedocs.io/en/latest/}\cr

    D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P.
    Turnbaugh, E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
    \emph{Detecting novel associations in large datasets}.  Science 334,
    6062\cr

    \url{http://www.exploredata.net}\cr

    (SOM: Supplementary Online Material at
    \url{http://www.sciencemag.org/content/suppl/2011/12/14/334.6062.1518.DC1})
}
\keyword{ package }
\keyword{ network }
\keyword{ reproducibility }
\keyword{ stability }
\keyword{ inference }
\keyword{ distance }
\seealso{
  Further information on specific functions for inferring the adjacency matrix could be found in \pkg{WGCNA},
  \pkg{minerva} and \pkg{minet}.
}
%% \examples{
%% ~~ simple examples of the most important functions ~~
%% }
