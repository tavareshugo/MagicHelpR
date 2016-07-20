#' Download Arabidopsis MAGIC genotypes
#' 
#' This function downloads SNP genotypes for 1254 SNPs typed with GoldenGate 
#' technology and reported in Kover et. al (2009). Files are downloaded from
#' http://mus.well.ox.ac.uk/POOLING/ARABIDOPSIS/FOUNDER/GENOTYPES/
#'
#' @param save_dir directory to save all files.
#' @param tidy whether or not to "tidy" the retrieved files. 
#' This will call the `tidyArabMagic()` function. Default: TRUE
#'
#' @return files are downloaded to a directory of choice.
#' 
#' @rdname downloadArabMagic
#' @export
#'
#' @examples
#' downloadArabMagic('~/temp/magic_snp_1k')
downloadArabMagic <- function(save_dir, tidy = TRUE, example_data = FALSE){
	
	save_dir <- path.expand(save_dir)
	
	# Create the directory if it doesn't exist
	if(!dir.exists(save_dir)){
		cat("Target directory does not exist. It will be created.")
		dir.create(save_dir)
	}
	
	# Genotype files
	magic_server <- 'http://mus.well.ox.ac.uk/POOLING/ARABIDOPSIS/FOUNDER/GENOTYPES/'
	geno_founder <- 'founders.genotypes.txt'
	geno_magic <- 'magic.15012010.tar.gz'
	
	# Download the files
	cat("Downloading files from:\n", magic_server, "\n")
	download.file(file.path(magic_server, geno_founder), file.path(save_dir, geno_founder))
	download.file(file.path(magic_server, geno_magic), file.path(save_dir, geno_magic))
	
	# Untar the MAGIC genotypes
	untar(file.path(save_dir, geno_magic), exdir = save_dir)
	
	# Tidy the files if requested
	if(tidy) cat("Tidying files.\n"); tidyArabMagic(save_dir)
	
	# Download example phenotype data if required
	if(example_data){
		cat("Downloading example data from:\nhttp://mus.well.ox.ac.uk/magic/MAGIC.phenotype.example.12102015.txt\n")
		download.file("http://mus.well.ox.ac.uk/magic/MAGIC.phenotype.example.12102015.txt",
									file.path(save_dir, "magic_phenotype_example.txt"))
	} 
}


