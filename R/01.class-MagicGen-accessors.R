#' Retrieve phenotypes from a MagicGenPhen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return data.frame of phenotypes
#' @export
#' @rdname getPhenotypes
setGeneric("getPhenotypes", function(x) standardGeneric("getPhenotypes"))
setMethod("getPhenotypes", "MagicGenPhen", function(x) x@phenotypes)

#' Get SNP genotypes
#'
#' This function can be used to extract the genotypes for all the MAGIC lines in
#' a MagicGen or MagicGenPhen object. The genotypes can either be the SNP genotypes
#' or the probability matrix of genotype ancestry for each MAGIC line.
#'
#' @param x an object of class MagicGen
#' @param type which type of genotype to retrieve.
#'
#' @return genotypes from MagicGenPhen object
setGeneric("getGenotypes", function(x, type = c("probability", "allele")) standardGeneric("getGenotypes"))
setMethod("getGenotypes", "MagicGen", function(x, type = c("probability", "allele")){

	# Define which method to use for genotypes
	type <- type[1]
	if(any(!(type %in% c("probability", "allele")))){
		stop("type not recognised: ", type)
	}

	if(type == "probability") return(.getProbGenotypes(x))
	if(type == "allele") return(.getSnpGenotypes(x))
})

#' Retrieve SNP genotypes of the MAGIC lines from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotypes SNP.
setGeneric(".getSnpGenotypes", function(x) standardGeneric(".getSnpGenotypes"))
setMethod(".getSnpGenotypes", "MagicGen", function(x) x@snp_genotypes)


#' Retrieve genotype probabilities for the MAGIC lines from a MagicGen object
#'
#' @param x an object of class "MagicGen".
#'
#' @return a list of genotype probabilities for every SNP.
setGeneric(".getProbGenotypes", function(x) standardGeneric(".getProbGenotypes"))
setMethod(".getProbGenotypes", "MagicGen", function(x) x@prob_genotypes)


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
