#' Estimate founder QTL effect
#'
#' This function estimates the phenotypic effect of each accession's allele at
#' a particular marker. It is based on the probabilistic assignment of each
#' MAGIC line to a single accession. The procedure is repeated several times to
#' produce an average estimate and associated uncertainty. This function is based
#' on the function \code{imputed.one.way.anova()} from the scripts available at
#' http://mtweb.cs.ucl.ac.uk/mus/www/magic/
#'
#' @param x an object of class MagicData.
#' @param phenotype the phenotype to get the results from.
#' @param marker the SNP marker to estimate the effects for.
#' @param n_samples number of Monte Carlo samples to draw.
#' @param summarised whether or not to report summarised results. The summarised
#' output gives the mean, median and 2.5% and 97.5% percentiles of the phenotype
#' imputations. The two percentiles define the 95% range of the estimate, which
#' can be used as a confidence interval to report accuracy of the estimates.
#' Default: TRUE
#' @param standardised whether or not the phenotypic values should be standardised
#' to the mean. In that case the result reflects how many standard deviations
#' the estimated effect deviates from the observed trait mean. Default: TRUE
#' @param covariates optional name of column(s) from phenotypes to use as
#' covariate(s). The effect will be estimated from the residuals of a standard
#' linear model between the trait of interest and the specified covariates.
#' Default: NULL
#'
#' @export
#' @rdname estimateFounderEffect
setGeneric("estimateFounderEffect",
           function(x, phenotype, marker, n_samples = 500,
                    summarised = TRUE, standardised = TRUE,
                    covariates = NULL)
             standardGeneric("estimateFounderEffect"))



#' @return a data.frame with imputed phenotypes for each founder accession.
#' @export
#' @rdname estimateFounderEffect
#'
#' @examples
#' \dontrun{
#' estimateFounderEffect(magic_genotypes, "MN1_29291")
#' }
setMethod("estimateFounderEffect", "MagicGenPhen",
          function(x, phenotype, marker, n_samples = 500,
                   summarised = TRUE, standardised = TRUE,
                   covariates = NULL){

  phenotype <- getPhenotypes(x)[, phenotype]
  gen_prob <- getGenotypes(x)[[marker]]

  # If covariates are specified, get residuals of linear model
  if(!is.null(covariates)){
    covariates <- as.matrix(getPhenotypes(x)[,covariates])
    phenotype <- lm(phenotype ~ covariates) %>% resid()
  }

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

    data.frame(accession = rownames(S), n = N, effect = mean_pheno, stringsAsFactors = FALSE)
    

  }, phenotype) %>%
    bind_rows(.id = "rep")

  # Standardise the trait if requested
  if(standardised){
    estimated_effect <- estimated_effect %>%
      group_by(rep) %>%
      mutate(effect_standard = (effect - mean(phenotype))/sd(phenotype)) %>%
      ungroup()
  }

  # Summarise if requested
  if(summarised & standardised){

  	# Get the founder SNP genotypes and add to output
  	founder_gen <- getFounderGenotypes(x)[[marker]]
  	
  	# Summarise the imputed phenotype
    estimated_effect <- estimated_effect %>%
    group_by(accession) %>%
      summarise(effect_mean = mean(effect_standard, na.rm = TRUE),
                effect_up = quantile(effect_standard, 0.975, na.rm = TRUE),
                effect_lo = quantile(effect_standard, 0.025, na.rm = TRUE),
                n = median(n)) %>%
    	left_join(founder_gen, by = "accession") %>%
      mutate(accession = reorder(factor(accession), effect_mean)) %>%
      ungroup()


  } else if(summarised & !standardised){

  	# Get the founder SNP genotypes and add to output
  	founder_gen <- getFounderGenotypes(x)[[marker]]

  	estimated_effect <- estimated_effect %>%
      group_by(accession) %>%
      summarise(effect_mean = mean(effect, na.rm = TRUE),
                effect_up = quantile(effect, 0.975, na.rm = TRUE),
                effect_lo = quantile(effect, 0.025, na.rm = TRUE),
                n = median(n)) %>%
  		left_join(founder_gen, by = "accession") %>%
      mutate(accession = reorder(factor(accession), effect_mean)) %>%
      ungroup()
  }

  return(estimated_effect)

})








