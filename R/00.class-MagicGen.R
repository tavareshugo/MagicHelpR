#' MagicGen class
#'
#' @slot markers data.frame with marker names and positions. 
#' @slot genotypes list with genotype probability matrices. 
#'
#' @rdname class-MagicGen
#' @aliases MagicGen
setClass("MagicGen",
				 representation(markers = "data.frame", 
				 							 genotypes = "list"))


#' MagicGenPhen class
#'
#' @slot phenotypes data.frame. 
#' @slot markers data.frame. 
#' @slot genotypes list. 
#'
#' @rdname class-MagicGen
#' @aliases MagicGenPhen
setClass("MagicGenPhen",
         representation(phenotypes = "data.frame",
                        markers = "data.frame", 
                        genotypes = "list"),
         contains = "MagicGen")


