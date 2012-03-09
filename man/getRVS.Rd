\name{getRVS}
\alias{getRVS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
acquire data for and construct a ragged variant set instance
}
\description{
acquire data for and construct a ragged variant set instance
}
\usage{
getRVS(packname, fns2samplenames = function(x)
 gsub(".*(NA.....).*", "\\1", x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{packname}{
string naming package where the resources are found
}
 \item{fns2samplenames}{function to transform filenames
to sample name tokens}
}
\details{
currently very specialized, as the  protocol for managing collections
of VCF files with discrepant variant sets per subject is not clear

assumes the package has inst/rowdata where row data of 
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (packname) 
{
    tmp = new("raggedVariantSet", filenames = dir(system.file("rowdata", 
        package = packname), full = TRUE))
    tmp@sampleNames = gsub(".*(NA.....).*", "\\1", tmp@filenames)
    names(tmp@filenames) = tmp@sampleNames
    tmp
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
