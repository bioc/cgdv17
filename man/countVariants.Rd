\name{countVariants}
\alias{countVariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
count variants in a raggedVariantSet instance
}
\description{
count variants in a raggedVariantSet instance
}
\usage{
countVariants(rvs, delim, qthresh = 160, applier = lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rvs}{
instance of \code{\linkS4class{raggedVariantSet}}
}
  \item{delim}{
GRanges instance
}
  \item{qthresh}{
quality threshold for keeping a variant in count
}
  \item{applier}{
lapply-like function
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
