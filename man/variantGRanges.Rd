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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
