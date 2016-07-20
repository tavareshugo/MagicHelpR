#' Retrive phenotypes from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return data.frame of phenotypes
#' @export
#' @rdname getPhenotypes
setGeneric("getPhenotypes", function(x) standardGeneric("getPhenotypes"))
setMethod("getPhenotypes", "MagicGen", function(x) x@phenotypes)


#' Retrive genotypes from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotype probabilities for every SNP.
#' @export
#' @rdname getGenotypes
setGeneric("getGenotypes", function(x) standardGeneric("getGenotypes"))
setMethod("getGenotypes", "MagicGen", function(x) x@genotypes)


#' Retrive markers from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a data.frame of marker names and positions.
#' @export
#' @rdname getMarkers
setGeneric("getMarkers", function(x) standardGeneric("getMarkers"))
setMethod("getMarkers", "MagicGen", function(x) x@markers)


# Set show method
setMethod("show", "MagicGen", function(object){
	phen <- getPhenotypes(object)
	gen <- getGenotypes(object)
	markers <- getMarkers(object)
	
	# Print information
	cat("Object of class", class(object), "\n")
	cat(nrow(phen), "MAGIC lines with a phenotype.\n")
	cat("There are", ncol(phen)-1, "phenotypes:\n")
	cat("Data exists for", length(gen), "markers.\n")
})
