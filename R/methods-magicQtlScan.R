#' Make a QTL scan
#'
#' @param x an object of class MagicData.
#' @param phenotype to run QTL analysis on.
#' @param covariates optional name of column(s) from phenotypes to use as covariate(s). Default: NULL
#' @param snp_cond optional ID of marker to use as covariate in the model. Default: NULL
#' @param h1 optional optional model to test. Need to provide `h0` as well. Default: NULL
#' @param h0 optional optional null model to test against. Need to provide `h1` as well. Default: NULL
#' @param perm optional number of permutations to run. If provided will permute the phenotypes 
#' the number of requested times and then run the QTL scan. It will add a new column to
#' the output with a genome-wide p-value (proportion of times where observed F < permuted F). Default: NULL
#' @param cores number of cores to use (Windows only supports 1 core). Default: 1
#'
#' @rdname magicQtlScan
setGeneric("magicQtlScan", function(x, phenotype, covariates = NULL, snp_cond = NULL, 
																		h1 = NULL, h0 = NULL, perm = NULL, cores = 1) 
	standardGeneric("magicQtlScan"))


#' @return a data.frame with the QTL results.
#' @export
#' @docType methods
#' @rdname magicQtlScan
#'
#' @examples
#' ...
setMethod("magicQtlScan", "MagicGen", 
					function(x, phenotype, covariates = NULL, snp_cond = NULL, 
									 h1 = NULL, h0 = NULL, perm = NULL, cores = 1){
	
	cat("Phenotype:", phenotype, "\n")
	
	# Get marker and genotype information
	markers <- getMarkers(x)
	genotypes <- getGenotypes(x)
	
	# Get phenotype and covariate information
	phenotype <- as.matrix(getPhenotypes(x)[, phenotype])
	
	if(!is.null(covariates)) covariates <- as.matrix(getPhenotypes(x)[, covariates])
	
	# Check if there is a SNP to add to the model
	if(!is.null(snp_cond)){
		if(snp_cond %in% names(genotypes)){
			
			if(length(snp_cond) > 2) stop("Only one SNP can be added to the model.")
			snp_cond <- genotypes[[snp_cond]]
			
		} else {
			stop("SNP was not found in genotype matrices")
		}
	}
	
	# Model formulas for `lm`, based on user inputs
	qtl_models <- .makeQtlModel(covariates, snp_cond, h1, h0)
	
	# Function to fit the linear model to all SNPs in parallel
	allSnpScan <- function(PHEN = phenotype, GEN = genotypes, COV = covariates,
												 SNP = snp_cond, H1 = qtl_models$h1, H0 = qtl_models$h0, Cores = cores){
		parallel::mclapply(GEN, .qtlFit, PHEN, COV, SNP, H1, H0, 
											 mc.cores = Cores) %>%
			bind_rows(.id = "marker")
	}
	
	# Fit linear model to all SNPs
	qtl_scan <- allSnpScan()
	
	# Merge with marker table to get coordinates of each marker
	qtl_scan <- merge(markers, qtl_scan)
	
	
	##### Permutations #####
	if(!is.null(perm) & perm > 0){
		cat("Performing", perm, "permutations. This might take a while...\n")
		
		# Shuffle phenotypes and covariates
		perm_index <- replicate(perm, sample(length(phenotype)))

		# Run QTL scan on each permuted phenotype (shows a progress bar)
		qtl_perm <- plyr::adply(perm_index, 2, function(i, phenotype, covariates){
			phenotype <- phenotype[i]
			covariates <- covariates[i, ]
			allSnpScan(PHEN = phenotype, COV = covariates)
			}, phenotype, covariates, .progress = "text", .id = "perm_rep")
		
		# Get the minimum p-value for each permutation
		## Note: before I was taking the max F-value, but degrees freedom 
		## vary when there is multicolinearity in genotype matrix
		min_p <- qtl_perm %>%
			group_by(perm_rep) %>%
			summarise(min_p = min(p)) %>%
			.$min_p

		# Get a genome-wide p-value by checking how many times 
		# the observed p is above the permuted p
		qtl_scan$p_perm <- sapply(qtl_scan$p, function(x, min_p) sum(x > min_p)+1, min_p)/(perm+1)
	}
	
	return(qtl_scan)
	
})


#' Specify the null and alternative models for QTL scan
#' 
#' @param COV covariates
#' @param SNP snp genotype probabilities
#' @param h1 alternative model
#' @param h0 null model
.makeQtlModel <- function(COV, SNP, h1, h0){
	
	# Check if user specified their own models
	if(any(!is.null(h1), !is.null(h0))){
		
		# Fail if either H1 or H0 were not provided
		if(any(is.null(h1), is.null(h0))){
			stop("Both h1 or h0 have to be specified.")
		}
		
		cat("F-test comparing the models:\nH1:", h1, "\nH0:", h0, "\n")
		
		return(list(h1 = h1, h0 = h0))
		
	} else if(is.null(COV) & is.null(SNP)){
		
		cat("F-test comparing the models:\nH1: PHEN ~ GEN\nH0: PHEN ~ 1\n")
		
		h1 <- "PHEN ~ GEN"
		h0 <- "PHEN ~ 1"

	} else {
		h1 <- paste0("PHEN ~ GEN")
		h0 <- paste0("PHEN ~ ")
		
		if(!is.null(COV)){
			h1 <- paste0(h1, " + COV")
			h0 <- paste0(h0, "COV")
		}
		if(!is.null(SNP)){
			h1 <- paste0(h1, " + SNP")
			if(!is.null(COV)){
				h0 <- paste0(h0, " + SNP")
			} else{ h0 <- paste0(h0, "SNP") }
			
		}
		
		cat("F-test comparing the models:\nH1:", h1, "\nH0:", h0, "\n")
	}
	
	return(list(h1 = h1, h0 = h0))
}


#' Fit a linear model to genotype data
#'
#' @param GEN genotype matrix.
#' @param PHEN phenotype vector.
#' @param COV covariate matrix or vector.
#' @param SNP snp ID to use as covariate.
#'
#' @return data.frame with test results.
.qtlFit <- function(GEN, PHEN, COV, SNP, h1, h0){
	
	fit1 <- lm(as.formula(h1))
	fit0 <- lm(as.formula(h0))
	
	return(data.frame(f = anova(fit0, fit1)[[5]][2],
										p = anova(fit0, fit1)[[6]][2],
										r2_h1 = summary(fit1)$adj.r.squared,
										r2_h0 = summary(fit0)$adj.r.squared))
}

