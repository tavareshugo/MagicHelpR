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
#' @param id the name of the variable with IDs in the data.frame passed to `phenotypes` argument.
#'
#' @return an object of type "MagicGen"
#' @export
#'
#' @examples
#' ...
magicFounderReconstruct <- function(ped, map, alleles, phenotypes, id){
	
	# Get IDs of MAGIC lines
	phenotypes <- rename_(phenotypes, magic = id)
	ids <- as.character(phenotypes$magic)
	
	# Make sure the path is expanded (no tilde and so on)
	ped <- path.expand(ped)
	map <- path.expand(map)
	alleles <- path.expand(alleles)
	
	# Make haplotype inference with happy function
	h <- happy.hbrem::happy(ped, alleles, gen = 7, file = "ped", 
														mapfile = map, haploid = TRUE)
	
	# Check that all IDs exist
	if(any(!(ids %in% h$subjects))){
		stop(paste("Some IDs do not have a genotype:", ids[which(!(ids %in% h$subjects))]))
	}
	
	# Calculate genotype probabilities for each SNP
	geno_prob <- lapply(h$additive$genome$marker, function(marker, h, ids){
		probs <- happy.hbrem::hdesign(h, marker)
		
		rownames(probs) <- h$subjects
		
		return(probs[ids, ])
		
	}, h, ids)
	
	names(geno_prob) <- h$additive$genome$marker
	
	# Output an object of class "MagicGen"
	out <- new("MagicGen",
						 phenotypes = phenotypes,
						 markers = h$additive$genome,
						 genotypes = geno_prob)
	
	return(out)
}

