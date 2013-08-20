\name{netSI}
\alias{netSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network Stability Indicators
}
\description{
Family NetSI of stability indicators
}
\usage{
netSI(d,indicator="S", dist='HIM', adj.method='cor', method="montecarlo", k=3, h=20, n.cores=NULL, save=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ dataset used for network inference and stability indicators computation}
  \item{method}{ String identifying the cross validation schema to be used for the stability computation.
        One of 'montecarlo' (default), 'LOO' for Leave One Out or 'kCV'for k-fold Cross Validation.}
  \item{k}{ Numeric value indicating the number of subset the data will be splitted in.
    By default it is set to 3. 
    ##In case method=="montecarlo" the 1-1/cv samples are taken for the computation
    ##In case cvlab=="loo" it should be set to 1
    ##In case cvlab=="kcv" cv iteration are performed, each dividing the data into 1-1/cv
    }
  \item{h}{ Numeric value indicates the number of iterations (cross validations). By default it is set to 20}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
