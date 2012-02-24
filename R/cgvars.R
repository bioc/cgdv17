setClass("raggedVariantSet", representation(filenames="character",
  sampleNames="character"))
setMethod("show", "raggedVariantSet", function(object){
 cat("raggedVariantSet instance with", length(object@filenames), "elements.\n", sep=" ")
 cat("some sampleNames: ", selectSome(gsub(".*(NA.....).*", "\\1", object@filenames)), "\n")
})

getRVS = function(packname)
 {
 tmp = new("raggedVariantSet", filenames=dir(system.file("rowdata", package=packname), full=TRUE))
 tmp@sampleNames = gsub(".*(NA.....).*", "\\1", tmp@filenames)
 names(tmp@filenames) = tmp@sampleNames
 tmp
 }

getrd = function(x, id) {
 ans = try(get(load(x@filenames[id])))
 if (inherits(ans, "try-error")) stop(paste("could not find filename for id ", id))
 ans
}

setMethod("sampleNames", "raggedVariantSet", function(object) {
 object@sampleNames
})

setMethod("[", "raggedVariantSet", function(x, i, j, ..., drop=TRUE) {
 if (!missing(i)) cat("[ for raggedVariantSet ignores row subsetting...\n")
 sn = sampleNames(x)
 if (missing(j)) return(x)
 if (is(j, "numeric")) {
    return(
      new("raggedVariantSet", filenames=x@filenames[j],
             sampleNames=x@sampleNames[j]) 
    )}
 if (is(j, "character")) {
    kp = match(j, sampleNames(x))
    new("raggedVariantSet", filenames=x@filenames[kp],
             sampleNames=x@sampleNames[kp]) 
 }
})

#.countVariants = function(rvs, delim, qthresh=160,
#   applier=sapply) {
#  nms = sampleNames(rvs)
#  ans = applier( nms, function(x) {
#    thisrd = getrd(rvs, x)  # GRanges
#    indelim = thisrd[ ranges(thisrd) %in% ranges(delim) ]
#    curq = elementMetadata(indelim)$QUAL
#    indelim = indelim[ which(curq >= qthresh) ]
#    length(indelim)
#  })
#  names(ans) = nms
#  ans
#}

#setGeneric("countVariants", function(rvs, delim, qthresh=160, applier=sapply)
# standardGeneric("countVariants"))
#setMethod("countVariants", c("raggedVariantSet",
#   "GRanges", "numeric", "function"), 
#   function(rvs, delim, qthresh, applier) .countVariants(rvs, delim,qthresh,applier))
#setMethod("countVariants", c("raggedVariantSet",
#   "GRanges", "missing", "missing"), 
#   function(rvs, delim, qthresh, applier) .countVariants(rvs, delim,160,sapply))
#
#.variantNames = function(rvs, delim, qthresh=160,
#   applier=lapply) {
#  nms = sampleNames(rvs)
#  ans = applier( nms, function(x) {
#    thisrd = getrd(rvs, x)  # GRanges
#    indelim = thisrd[ ranges(thisrd) %in% ranges(delim) ]
#    curq = elementMetadata(indelim)$QUAL
#    indelim = indelim[ which(curq >= qthresh) ]
#    names(indelim)
#  })
#  names(ans) = nms
#  ans
#}

#setGeneric("variantNames", function(rvs, delim, qthresh=160, applier=lapply)
# standardGeneric("variantNames"))
#setMethod("variantNames", c("raggedVariantSet",
#   "GRanges", "numeric", "function"), 
#   function(rvs, delim, qthresh, applier) .variantNames(rvs, delim,qthresh,applier))
#setMethod("variantNames", c("raggedVariantSet",
#   "GRanges", "missing", "missing"), 
#   function(rvs, delim, qthresh, applier) .variantNames(rvs, delim,160,lapply))

.variantGRanges = function(rvs, delim, qthresh=160,
   applier=lapply) {
  nms = sampleNames(rvs)
  ans = applier( nms, function(x) {
    thisrd = getrd(rvs, x)  # GRanges
    indelim = thisrd[ IRanges::"%in%"(ranges(thisrd), ranges(delim)) ]
    curq = elementMetadata(indelim)$QUAL
    indelim[ which(curq >= qthresh) ]
  })
  names(ans) = nms
  ans
}

setGeneric("variantGRanges", function(rvs, delim, qthresh=160, applier=lapply)
 standardGeneric("variantGRanges"))
setMethod("variantGRanges", c("raggedVariantSet",
   "GRanges", "numeric", "function"), 
   function(rvs, delim, qthresh, applier) .variantGRanges(rvs, delim,qthresh,applier))
setMethod("variantGRanges", c("raggedVariantSet",
   "GRanges", "missing", "missing"), 
   function(rvs, delim, qthresh, applier) 
      .variantGRanges(rvs, delim,160,lapply))

countVariants = function(rvs, delim, qthresh=160, applier=lapply) {
 sapply(variantGRanges(rvs, delim, qthresh, applier), length)
}
variantNames = function(rvs, delim, qthresh=160, applier=lapply) {
 lapply(variantGRanges(rvs, delim, qthresh, applier), names)
}
