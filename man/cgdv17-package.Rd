\name{cgdv17-package}
\alias{cgdv17-package}
\alias{cgdv17}
\alias{CY17}
\alias{h1}
\alias{popvec}
\docType{package}
\title{
Complete Genomics Diversity Panel, chr17 on 46 individuals
}
\description{
Complete Genomics Diversity Panel, chr17 on 46 individuals, illustrating
subject-specific variant sets
}
\details{
\tabular{ll}{
Package: \tab cgdv17\cr
Version: \tab 0.0.9\cr
Suggests: \tab \cr
Imports: \tab Biobase, IRanges\cr
Depends: \tab R (>= 2.14), VariantAnnotation, org.Hs.eg.db, methods\cr
License: \tab Artistic 2.0\cr
LazyLoad: \tab yes\cr
biocViews: \tab genetics\cr
Built: \tab R 2.15.0; ; 2012-03-09 12:45:57 UTC; unix\cr
}

Index:
\preformatted{
countVariants           count variants in a raggedVariantSet instance
getRVS                  acquire data for and construct a ragged variant
                        set instance
padToReference          create a snpStats SnpMatrix instance by padding
                        a ragged variant set to reference alleles
                        wherever a variant is not recorded
raggedVariantSet-class
                        Class '"raggedVariantSet"'
variantGRanges          acquire a list of GRanges recording variants
                        and locations
}

see vignette; CY17 is an ExpressionSet on individuals from CEU and YRI
overlapping with the diversity set, popvec enumerates source populations,
h1 is an exemplar VCF header structure
}
\author{
VJ Carey <stvjc@channing.harvard.edu>

Maintainer: VJ Carey <stvjc@channing.harvard.edu>
}
\references{
}
\keyword{ package }
\seealso{
}
\examples{
}
