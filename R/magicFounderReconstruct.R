#' Calculate MAGIC founder genotype probabilities
#'
#' This function calls the function `happy.hbrem::happy()` to calculate 
#' the probability matrix of descent from each founder of the MAGIC lines
#' at each marker SNP.
#'
#' @param snp_dir directory with all genotype files. The function will be looking 
#' for files named "all_chr.MAGIC.happy.data", "all_chr.MAGIC.happy.map", 
#' "all_chr.MAGIC.happy.alleles", "founder_genotypes.tsv" and "magic_genotypes.tsv". 
#' All these files can be created with either the \code{tidyArabMagic} or
#' \code{downloadArabMagic} functions.
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
magicFounderReconstruct <- function(snp_dir, phenotypes = NULL, id = NULL){
  
  ##### Check all files #####
  # File paths
  # Note: need to make sure the path is expanded (no tilde for happy.hbrem functions)
  ped <- path.expand(file.path(snp_dir, "all_chr.MAGIC.happy.data"))
  map <- path.expand(file.path(snp_dir, "all_chr.MAGIC.happy.map"))
  alleles <- path.expand(file.path(snp_dir, "all_chr.MAGIC.happy.alleles"))
  snps <- path.expand(file.path(snp_dir, "magic_genotypes.tsv"))
  founders <- path.expand(file.path(snp_dir, "founder_genotypes.tsv"))
  
  # Check that all files exist
  if(!file.exists(ped)) stop("all_chr.MAGIC.happy.data file was not found")
  if(!file.exists(map)) stop("all_chr.MAGIC.happy.map file was not found")
  if(!file.exists(alleles)) stop("all_chr.MAGIC.happy.alleles file was not found")
  if(!file.exists(snps)) stop("magic_genotypes.tsv file was not found")
  if(!file.exists(founders)) stop("founder_genotypes.tsv file was not found")
  
  
  ##### Haplotype inference ######
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
  
  
  ##### Get SNP genotypes #####
  # Read SNP genotypes for MAGIC lines
  geno_snp <- read.table(snps, header = TRUE, stringsAsFactors = FALSE)
  head(geno_snp)
  
  # Split the genotypes into a list
  geno_snp <- geno_snp %>% 
    filter(marker %in% names(geno_prob)) %>% 
    split(.$marker)
  
  # Convert the genotypes into a list of vectors
  geno_snp <- lapply(geno_snp, function(x){
    # Make SNP genotypes the same order as in the happy output
    x <- x[match(h$subjects, x$magic),]
    
    # Convert into a matrix
    x <- as.matrix(x)
    
    # Add names of individuals as rownames
    rownames(x) <- x[,1]
    
    # Return only the relevant genotype column
    x[,3]
  })
  
  
  ##### Get founder genotypes #####
  # Read table with genotypes
  geno_founder <- read.table(founders, header = TRUE, stringsAsFactors = FALSE)
  
  # Split into a list
  geno_founder <- split(geno_founder, geno_founder$marker)
  
  
  ##### Create MagicGen object #####
  # Prepare object of class "MagicGen"
  out <- new("MagicGen",
             markers = h$additive$genome,
             prob_genotypes = geno_prob,
             snp_genotypes = geno_snp,
             founder_genotypes = geno_founder)
  
  
  ##### Add phenotypes if requested #####
  
  # Add phenotypes if requested
  if(!is.null(phenotypes) & !is.null(id)){
    out <- addPhenotypes(out, phenotypes, id)
  }
  
  return(out)
}

