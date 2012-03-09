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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rv, gr, qthresh = 160, applier = lapply) 
{
    vv = variantGRanges(rv, gr, qthresh = qthresh, applier = applier)
    lens = lapply(vv, length)
    if (any(lens == 0)) 
        warning("at least one sample gave no variants with these parameters")
    nms = lapply(vv, names)
    allv = unique(unlist(nms))
    padmat = matrix("0/0", nr = length(allv), nc = length(lens))
    rownames(padmat) = allv
    colnames(padmat) = sampleNames(rv)
    fixup = function(z) {
        z = gsub("\\|", "/", z)
        if (!all(z \%in\% c("0/0", "0/1", "1/1", "1/0"))) 
            z[!(z \%in\% c("0/0", "0/1", "1/1", "1/0"))] = "0/0"
        z
    }
    refl = list()
    altl = list()
    varnl = list()
    for (i in 1:length(lens)) {
        curg = elementMetadata(vv[[i]])$geno
        REF = elementMetadata(vv[[i]])$REF
        ALTIN = elementMetadata(vv[[i]])$ALT
        DS = function(x) {
            x = IRanges::unlist(x)
            if (!all(x \%in\% c("A", "C", "G", "T"))) 
                x[which(!(x \%in\% c("A", "C", "G", "T")))] = "N"
            DNAStringSet(x)
        }
        ALT = IRanges:::newCompressedList("DNAStringSetList", 
            DS(ALTIN), end = as.integer(elementLengths(ALTIN)))
        padmat[names(curg), i] = fixup(curg)
        refl[[i]] = REF
        altl[[i]] = ALT
        varnl[[i]] = names(curg)
    }
    refv = unlist(lapply(refl, as.character))
    altv = unlist(IRanges::lapply(altl, function(x) as.character(IRanges::unlist(x))))
    names(refv) = unlist(varnl)
    names(altv) = unlist(varnl)
    ref = DNAStringSet(refv[rownames(padmat)])
    altv = altv[rownames(padmat)]
    alt = IRanges:::newCompressedList("DNAStringSetList", DNAStringSet(altv), 
        end = as.integer(cumsum(nchar(altv))))
    VariantAnnotation::MatrixToSnpMatrix(padmat, ref, alt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
