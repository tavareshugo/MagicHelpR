#' Make a QTL scan
#'
#' @param x an object of class MagicGenPhen
#' @param phenotype to run QTL analysis on.
#' @param covariates optional name of column(s) from phenotypes to use as covariate(s). Default: NULL
#' @param snp_cond optional ID of marker(s) to use as covariate(s) in the model. Default: NULL
#' @param h1 optional optional model to test. Need to provide `h0` as well. Default: NULL
#' @param h0 optional optional null model to test against. Need to provide `h1` as well. Default: NULL
#' @param perm optional number of permutations to run. If provided will permute the phenotypes 
#' the number of requested times and then run the QTL scan. It will add a new column to
#' the output with a genome-wide p-value (proportion of times where observed p > permuted p). Default: 0
#' @param cores number of cores to use (Windows only supports 1 core). Default: 1
#'
#' @rdname magicQtlScan
setGeneric("magicQtlScan", function(x, phenotype, covariates = NULL, snp_cond = NULL, 
																		h1 = NULL, h0 = NULL, perm = 0, cores = 1) 
	standardGeneric("magicQtlScan"))


#' @return a data.frame with the QTL results.
#' @export
#' @docType methods
#' @rdname magicQtlScan
#'
#' @examples
#' ...
setMethod("magicQtlScan", "MagicGenPhen", 
					function(x, phenotype, covariates = NULL, snp_cond = NULL, 
									 h1 = NULL, h0 = NULL, perm = 0, cores = 1){
	
	
	# Prepare variables for model
	model_vars <- .prepareModelVariables(x, phenotype, covariates, snp_cond)
	
	# Prepare model formulas for `lm`, based on user inputs
	qtl_models <- .makeQtlModel(phenotype, covariates, snp_cond, h1, h0)
	
	# Function to fit the linear model to all SNPs using parallel mclapply
	allSnpScan <- function(VARS = model_vars, GEN = getGenotypes(x), 
												 H1 = qtl_models$h1, H0 = qtl_models$h0, Cores = cores){
		parallel::mclapply(GEN, .qtlFit, VARS, H1, H0, 
											 mc.cores = Cores) %>%
			bind_rows(.id = "marker")
	}
	
	# Fit linear model to all SNPs
	t1 <- proc.time()  # to report time taken
	qtl_scan <- allSnpScan()
	t2 <- proc.time() - t1
	cat("Elapsed time: ", t2[3], "seconds\n")
	
	# Merge with marker table to get coordinates of each marker
	qtl_scan <- merge(getMarkers(x), qtl_scan)
	
	
	##### Permutations #####
	if(perm > 0){
		cat("Performing", perm, "permutations. This might take a while...\n")
		
		# Shuffle phenotypes and covariates
		perm_index <- replicate(perm, sample(length(model_vars[[phenotype]])))

		# Run QTL scan on each permuted phenotype (shows a progress bar)
		qtl_perm <- plyr::adply(perm_index, 2, function(i, model_vars){
			# First shuffle the co-variates
			model_vars <- lapply(model_vars, function(x){
				if(is.vector(x)) return(x[i])
				if(is.matrix(x)) return(x[i,])
			})
			
			# Now run QTL model with shuffled variables
			allSnpScan(VARS = model_vars)
			}, model_vars, .progress = "text", .id = "perm_rep")
		
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



#' Internal MagicHelpR function
#' 
#' Make a list of the variables to include in the linear model
#'
#' @param x object of class MagicGen
#' @param phenotype phenotype variable
#' @param covariates covariate variables
#' @param snp_cond covariate SNPs
.prepareModelVariables <- function(x, phenotype, covariates, snp_cond){
	
	# Check that all requested variables exist
	if(!(phenotype %in% colnames(getPhenotypes(x)))) stop("Phenotype ", phenotype, " does not exist.")
	if(any(!(covariates %in% colnames(getPhenotypes(x))))) stop("Some covariates do not exist.")
	if(any(!(snp_cond %in% names(getGenotypes(x))))) stop("Some specified snp_cond do not exist.")
	
	# Make lists of phenotype, covariates and snp matrices
	phenotype_list <- lapply(phenotype, function(i, x) getPhenotypes(x)[,i], x = x)
	covariates_list <- lapply(covariates, function(i, x) getPhenotypes(x)[,i], x = x)
	snp_list <- lapply(snp_cond, function(i, x) getGenotypes(x)[[i]], x = x)
	
	# Concatenate the lists and name elements appropriately
	model_vars <- c(phenotype_list, covariates_list, snp_list)
	names(model_vars) <- c(phenotype, covariates, snp_cond)
	
	return(model_vars)
	
}



#' Internal MagicHelpR function
#' 
#' Specify the null and alternative models for QTL scan.
#' 
#' @param phenotype phenotype variable
#' @param covariates covariate variables
#' @param snp_cond IDs of SNPs to use as covariates
#' @param h1 alternative model
#' @param h0 null model
.makeQtlModel <- function(phenotype, covariates, snp_cond, h1, h0){
	
	# Check if user specified their own models
	if(any(!is.null(h1), !is.null(h0))){
		
		# Fail if either H1 or H0 were not provided
		if(any(is.null(h1), is.null(h0))){
			stop("Both h1 or h0 have to be specified.")
		}
		
		cat("\nF-test comparing the models:\nH1: ", h1, "\nH0: ", h0, "\n")
		
		return(list(h1 = h1, h0 = h0))
		
	} else if(is.null(covariates) & is.null(snp_cond)){
		
		h1 <- paste(phenotype, "~ GEN")
		h0 <- paste(phenotype, "~ 1")
		
		cat("\nF-test comparing the models:\nH1: ", h1, "\nH0: ", h0, "\n")
		
	} else {
		h1 <- paste0(phenotype, " ~ GEN")
		h0 <- paste0(phenotype, " ~ ")
		
		if(!is.null(covariates)){
			COV <- paste(covariates, collapse = " + ")
			h1 <- paste(h1, "+", COV)
			h0 <- paste0(h0, COV)
		}
		if(!is.null(snp_cond)){
			SNP <- paste(snp_cond, collapse = " + ")
			h1 <- paste0(h1, " + ", SNP)
			if(!is.null(covariates)){
				h0 <- paste0(h0, " + ", SNP)
			} else{ h0 <- paste0(h0, SNP) }
			
		}
		
		cat("\nF-test comparing the models:\nH1:", h1, "\nH0:", h0, "\n")
	}
	
	return(list(h1 = h1, h0 = h0))
}



#' Fit a linear model to genotype data
#'
#' @param GEN genotype matrix.
#' @param model_vars list of variables to include in model.
#' @param h1 alternative model.
#' @param h0 null model.
#'
#' @return data.frame with test results.
.qtlFit <- function(GEN, model_vars, h1, h0){
	
	fit1 <- with(model_vars, lm(as.formula(h1)))
	fit0 <- with(model_vars, lm(as.formula(h0)))
	
	qtl_aov <- anova(fit0, fit1)

	return(data.frame(f = qtl_aov$F[2],
										df1 = qtl_aov$Df[2],
										df2 = qtl_aov$Res.Df[2],
										p = qtl_aov$`Pr(>F)`[2],
										r2_h1 = summary(fit1)$adj.r.squared,
										r2_h0 = summary(fit0)$adj.r.squared))
}

