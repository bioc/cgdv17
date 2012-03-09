\name{variantGRanges}
\alias{variantGRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
acquire a list of GRanges recording variants and locations
}
\description{
acquire a list of GRanges recording variants and locations
}
\usage{
variantGRanges(rvs, delim, qthresh = 160, applier = lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rvs}{
\code{\linkS4class{raggedVariantSet}} instance
}
  \item{delim}{
GRanges instance for confinement
}
  \item{qthresh}{
lower bound on quality
}
  \item{applier}{
lapply like function
}
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
VJ Carey <stvjc@channing.harvard.edu>
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
structure(function (rvs, delim, qthresh = 160, applier = lapply) 
standardGeneric("variantGRanges"), generic = structure("variantGRanges", package = "cgdv17"), package = "cgdv17", group = list(), valueClass = character(0), signature = c("rvs", 
"delim", "qthresh", "applier"), default = `NULL`, skeleton = function (rvs, 
    delim, qthresh = 160, applier = lapply) 
stop("invalid call in method dispatch to \"variantGRanges\" (no default method)", 
    domain = NA)(rvs, delim, qthresh, applier), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
