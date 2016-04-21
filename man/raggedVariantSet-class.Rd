\name{raggedVariantSet-class}
\Rdversion{1.1}
\docType{class}
\alias{raggedVariantSet-class}
\alias{[,raggedVariantSet,ANY,ANY,ANY-method}
\alias{sampleNames,raggedVariantSet-method}
\alias{show,raggedVariantSet-method}
\alias{variantGRanges,raggedVariantSet,GRanges,missing,missing-method}
\alias{variantGRanges,raggedVariantSet,GRanges,numeric,function-method}
\alias{[,raggedVariantSet,ANY,ANY-method}

\title{Class \code{"raggedVariantSet"}}
\description{
manage information on non-aligned variant sets from multiple VCFs
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("raggedVariantSet", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{filenames}:}{files will be held in inst/rowranges, named here}
    \item{\code{sampleNames}:}{names of samples managed}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "raggedVariantSet", i = "ANY", j = "ANY", drop = "ANY")}: familiar subsetting syntax}
    \item{sampleNames}{\code{signature(object = "raggedVariantSet")}: getter }
    \item{show}{\code{signature(object = "raggedVariantSet")}: concise report }
    \item{variantGRanges}{\code{signature(rvs = "raggedVariantSet", delim = "GRanges", qthresh = "missing", applier = "missing")}: getter }
    \item{variantGRanges}{\code{signature(rvs = "raggedVariantSet", delim = "GRanges", qthresh = "numeric", applier = "function")}: getter with quality threshold }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("raggedVariantSet")
}
\keyword{classes}
