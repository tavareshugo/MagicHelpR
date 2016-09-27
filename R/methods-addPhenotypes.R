#' Add phenotypes to 'MagicGen' object
#'
#' This function adds phenotypes to an existing 'MagicGen' object. 
#' The genotypes are kept only for those individuals with a genotype. 
#' 
#' @param x an object of class 'MagicGen'.
#' @param phenotypes a 'data.frame' object with phenotypes.
#' @param id the column name containing the MAGIC IDs in the 'phenotypes' table.
#'
#' @rdname addPhenotypes
setGeneric("addPhenotypes", 
           function(x, phenotypes, id) 
             standardGeneric("addPhenotypes"))


#' @return a a 'MagicGen' object with phenotypes added to it.
#' @export
#' @docType methods
#' @rdname addPhenotypes
#'
#' @examples
#' ...to be added...
setMethod("addPhenotypes", "MagicGen", 
          function(x, phenotypes, id){
            
            # Stop if there are already phenotypes
            if(is(x, "MagicGenPhen")){
              stop("Phenotypes are already present. Consider re-running 'magicFounderReconstruct()'.")
            }
            
            # Get IDs of MAGIC lines
            phenotypes <- rename_(phenotypes, magic = id)
            ids <- as.character(phenotypes$magic)
            
            # Check that all IDs exist
            if(any(!(ids %in% rownames(x@genotypes[[1]])))){
              stop(paste("Some IDs do not have a genotype:", ids[which(!(ids %in% h$subjects))], "\n",
                         "Consider removing them from the phenotype table."))
            }
            
            # Subset genotype probability matrices for each SNP
            geno_prob <- lapply(getGenotypes(x), function(probs, ids){
              return(probs[ids, ])
            }, ids)
            
            # Update "MagicGen" object
            out <- new("MagicGenPhen",
                       phenotypes = phenotypes,
                       markers = getMarkers(x),
                       genotypes = geno_prob)
            
            return(out)

            })








