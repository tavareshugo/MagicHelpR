#' Estimate founder effect based on probabilistic sampling
#'
#' @param x an object of class MagicData.
#' @param phenotype the phenotype to get the results from.
#' @param marker the SNP marker to estimate the effects for.
#' @param n_samples number of Monte Carlo samples to draw.
#'
#' @rdname estimateFounderEffect
setGeneric("estimateFounderEffect", 
           function(x, phenotype, marker, n_samples = 500) 
             standardGeneric("estimateFounderEffect"))



#' @return a data.frame with phenotype estimates for each founder accession. 
#' @export
#' @docType methods
#' @rdname estimateFounderEffect
#'
#' @examples
#' ...
setMethod("estimateFounderEffect", "MagicGenPhen", 
          function(x, phenotype, marker, n_samples = 500){
  
  phenotype <- getPhenotypes(x)[, phenotype]
  gen_prob <- getProbGenotypes(x)[[marker]]
  
  # Make a design matrix for founders
  founder_matrix <- matrix(0, ncol = 19, nrow = 19)
  diag(founder_matrix) <- 1
  colnames(founder_matrix) <- colnames(gen_prob)
  rownames(founder_matrix) <- colnames(gen_prob)
  
  # For each MAGIC line, infer a founder based on its probability
  magic_founder <- apply(gen_prob, 1, function(x, n_samples) sample(19, n_samples, T, x), n_samples)
  
  # Make a design matrix for each inferred MAGIC founder
  dmat <- lapply(1:n_samples, function(i, magic_founder, founder_matrix){
    founder_matrix[ , magic_founder[i, ]]
  }, magic_founder, founder_matrix)
  
  # Calculate mean phenotype for each founder in each replicate sampling
  estimated_effect <- lapply(dmat, function(x, phen){
    N <- apply(x, 1, sum)
    S <- x %*% phen
    mean_pheno <- ifelse(N > 0, S/N, NA)
    
    data.frame(accession = rownames(S), n = N, mean_phenotype = mean_pheno)
    
  }, phenotype) %>%
    bind_rows()

})








