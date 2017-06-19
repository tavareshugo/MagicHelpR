#' Retrieve phenotypes from a MagicGenPhen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return data.frame of phenotypes
#' @export
#' @rdname getPhenotypes
setGeneric("getPhenotypes", function(x) standardGeneric("getPhenotypes"))
setMethod("getPhenotypes", "MagicGenPhen", function(x) x@phenotypes)


#' Retrieve SNP genotypes of the MAGIC lines from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotypes SNP.
#' @export
#' @rdname getSnpGenotypes
setGeneric("getSnpGenotypes", function(x) standardGeneric("getSnpGenotypes"))
setMethod("getSnpGenotypes", "MagicGen", function(x) x@snp_genotypes)


#' Retrieve genotype probabilities for the MAGIC lines from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotype probabilities for every SNP.
#' @export
#' @rdname getProbGenotypes
setGeneric("getProbGenotypes", function(x) standardGeneric("getProbGenotypes"))
setMethod("getProbGenotypes", "MagicGen", function(x) x@prob_genotypes)


#' Retrieve SNP genotypes of MAGIC founders from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotypes for each founder accession.
#' @export
#' @rdname getFounderGenotypes
setGeneric("getFounderGenotypes", function(x) standardGeneric("getFounderGenotypes"))
setMethod("getFounderGenotypes", "MagicGen", function(x) x@founder_genotypes)


#' Retrieve markers from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a data.frame of marker names and positions.
#' @export
#' @rdname getMarkers
setGeneric("getMarkers", function(x) standardGeneric("getMarkers"))
setMethod("getMarkers", "MagicGen", function(x) x@markers)


# Set show methods
setMethod("show", "MagicGen", function(object){
  gen <- getProbGenotypes(object)

  # Print information
  cat("Object of class", class(object), "\n")
  cat("Using genotypes for", length(gen), "markers.\n")
})

setMethod("show", "MagicGenPhen", function(object){
  phen <- getPhenotypes(object)
  gen <- getProbGenotypes(object)

  # Print information
  cat("Object of class", class(object), "\n")
  cat(nrow(phen), "MAGIC lines with a phenotype.\n")
  cat("There are", ncol(phen)-1, "phenotypes:\n", names(phen[,-1]), "\n")
  cat("Using genotypes for", length(gen), "markers.\n")
})
