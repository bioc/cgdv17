\name{getRVS}
\alias{getRVS}
\alias{getrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
acquire data for and construct a ragged variant set instance
}
\description{
acquire data for and construct a ragged variant set instance
}
\usage{
getRVS(packname, fns2samplenames = function(x)
 gsub(".*(NA.....).*", "\\\\1", x))

getrd(x, id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{packname}{
string naming package where the resources are found
}
 \item{fns2samplenames}{function to transform filenames
to sample name tokens}
 \item{x}{instance of raggedVariantSet}
 \item{id}{character to select sample}
}
\details{
currently very specialized, as the  protocol for managing collections
of VCF files with discrepant variant sets per subject is not clear

assumes the package has inst/rowranges where row ranges of 
\code{\link[VariantAnnotation]{readVcf}} results are held
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
