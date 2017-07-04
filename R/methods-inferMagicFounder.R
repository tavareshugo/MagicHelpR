#' Infer most likely founder for a given SNP
#'
#' @param x an object of class MagicData.
#' @param marker name of SNP marker to infer founder.
#' @param prob_thr probability threshold for attributing a genotype to a founder.
#'
#' @rdname inferMagicFounder
setGeneric("inferMagicFounder",
           function(x, marker, prob_thr = 0.5)
             standardGeneric("inferMagicFounder"))


#' @return a data.frame with the most likely founder for each MAGIC line.
#' @export
#' @docType methods
#' @rdname inferMagicFounder
#'
#' @examples
#' ...
setMethod("inferMagicFounder", "MagicGenPhen",
          function(x, marker, prob_thr = 0.5){

  # Get the genotype probabilities for the specified marker
  gen_prob <- getGenotypes(x)[[marker]]
  gen_snp <- getGenotypes(x, type = "allele")[[marker]]

  # Convert genotype probabilities to binary using user's threshold
  gen_prob <- ifelse(gen_prob > prob_thr, 1, 0)

  # Report number of times each MAGIC line is represented
  n_magic <- rowSums(gen_prob)
  message(sum(n_magic == 1), " MAGIC lines were attributed to exactly one founder.")
  message(sum(n_magic == 0), " MAGIC lines were not attributed a founder.")
  message(sum(n_magic > 1), " MAGIC lines were attributed more than one founder.")

  # Presumed founder allele for each magic
  magic_founder <- apply(gen_prob, 2, function(x) ifelse(x == 1, names(x), NA)) %>%
    as.data.frame(stringsAsFactors = F) %>%
    gather(founder, magic, na.rm = T)

  # Make sure that founder is a factor with levels for all founders
  magic_founder$founder <- factor(magic_founder$founder, levels = colnames(gen_prob))

  # Merge with phenotypes
  magic_founder <- merge(magic_founder, getPhenotypes(x), by = "magic")

  # Add SNP genotypes
  magic_founder$allele <- gen_snp[magic_founder$magic]

  return(magic_founder)
})


