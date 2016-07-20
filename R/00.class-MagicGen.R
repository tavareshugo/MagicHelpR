#' MagicGen class
#'
#' @slot phenotypes data.frame with phenotypes for each MAGIC line. 
#' @slot markers data.frame with marker names and positions. 
#' @slot genotypes list with genotype probability matrices. 
#'
#' @rdname class-MagicGen
#' @aliases MagicGen
setClass("MagicGen",
				 representation(phenotypes = "data.frame",
				 							 markers = "data.frame", 
				 							 genotypes = "list"))
