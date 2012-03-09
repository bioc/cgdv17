\name{padToReference}
\alias{padToReference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a snpStats SnpMatrix instance by padding
a ragged variant set to reference alleles wherever
a variant is not recorded
}
\description{
create a snpStats SnpMatrix instance by padding
a ragged variant set to reference alleles wherever
a variant is not recorded
}
\usage{
padToReference(rv, gr, qthresh = 160, applier = lapply)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rv}{
\code{\linkS4class{raggedVariantSet}} instance
}
  \item{gr}{
GRanges instance
}
  \item{qthresh}{
quality lower bound for retention of variant
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
