setClass("raggedVariantSet", representation(filenames="character",
  sampleNames="character"))
setMethod("show", "raggedVariantSet", function(object){
 cat("raggedVariantSet instance with", length(object@filenames), "elements.\n", sep=" ")
 cat("some sampleNames: ", selectSome(gsub(".*(NA.....).*", "\\1", object@filenames)), "\n")
})

getRVS = function(packname="cgdv17") 
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

