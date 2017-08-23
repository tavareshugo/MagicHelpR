#' Add phenotypes to 'MagicGen' object
#'
#' This function adds phenotypes to an existing 'MagicGen' object.
#' The genotypes are kept only for those individuals with a genotype.
#'
#' @param x an object of class 'MagicGen'.
#' @param phenotypes a 'data.frame' object with phenotypes.
#' @param id the column name containing the MAGIC IDs in the 'phenotypes' table.
#'
#' @export
#' @rdname addPhenotypes
setGeneric("addPhenotypes",
           function(x, phenotypes, id)
             standardGeneric("addPhenotypes"))


#' @return a 'MagicGenPhen' object with phenotypes added to initial object.
#'
#' @export
#' @rdname addPhenotypes
#'
#' @examples
#' # to be added...
setMethod("addPhenotypes", "MagicGen",
          function(x, phenotypes, id){

            # Stop if there are already phenotypes
            if(is(x, "MagicGenPhen")){
              stop("Phenotypes are already present. Consider re-running 'magicFounderReconstruct()'.")
            }

            # Ensure x is a data.frame (to ensure compatibility with tbl objects)
            phenotypes <- as.data.frame(phenotypes)

            # Get IDs of MAGIC lines
            phenotypes <- rename_(phenotypes, magic = id)
            ids <- as.character(phenotypes$magic)

            # Check that all IDs exist
            if(any(!(ids %in% rownames(x@prob_genotypes[[1]])))){
              stop(paste("Some IDs do not have a genotype:", ids[which(!(ids %in% rownames(x@prob_genotypes[[1]])))], "\n",
                         "Consider removing them from the phenotype table."))
            }

            # Subset genotypes to include only those MAGIC lines with a phenotype
            geno_prob <- lapply(getGenotypes(x), function(g, ids){
              return(g[ids, ])
            }, ids)

            geno_snp <- lapply(getGenotypes(x, type = "allele"), function(g, ids){
              return(g[ids])
            }, ids)


            # Create new "MagicGenPhen" object
            out <- new("MagicGenPhen",
                       phenotypes = phenotypes,
                       markers = getMarkers(x),
                       prob_genotypes = geno_prob,
                       snp_genotypes = geno_snp,
                       founder_genotypes = getFounderGenotypes(x))

            return(out)

            })








