#' Calculate MAGIC founder genotype probabilities
#'
#' This function calls the function `happy.hbrem::happy()` to calculate 
#' the probability matrix of descent from each founder of the MAGIC lines
#' at each marker SNP.
#'
#' @param ped genotype file (.happy.ped).
#' @param map map file (.happy.map).
#' @param alleles alleles file (.happy.alleles).
#' @param phenotypes a data.frame with phenotypes and with one of the columns being the MAGIC line IDs. 
#' This is set to NULL by default, which returns a 'MagicGen' object with no phenotypes included.
#' @param id the name of the variable with IDs in the data.frame passed to `phenotypes` argument.
#' This is set to NULL by default, which returns a 'MagicGen' object with no phenotypes included.
#'
#' @return an object of type "MagicGen" or "MagicGenPhen".
#' @export
#'
#' @examples
#' ...
magicFounderReconstruct <- function(ped, map, alleles, phenotypes = NULL, id = NULL){
	
	# Make sure the path is expanded (no tilde and so on)
	ped <- path.expand(ped)
	map <- path.expand(map)
	alleles <- path.expand(alleles)
	
	# Make haplotype inference with happy function
	h <- happy.hbrem::happy(ped, alleles, gen = 7, file = "ped", 
														mapfile = map, haploid = TRUE)
	

	# Calculate genotype probabilities for each SNP
	geno_prob <- lapply(h$additive$genome$marker, function(marker, h){
		probs <- happy.hbrem::hdesign(h, marker)
		
		rownames(probs) <- h$subjects
		
		return(probs)
		
	}, h)
	
	names(geno_prob) <- h$additive$genome$marker
	
	# Prepare object of class "MagicGen"
	out <- new("MagicGen",
	           markers = h$additive$genome,
	           genotypes = geno_prob)
	
	# Add phenotypes if requested
	if(!is.null(phenotypes) & !is.null(id)){
	  out <- addPhenotypes(out, phenotypes, id)
	}
	
	return(out)
}

