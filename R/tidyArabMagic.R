#' Tidy Arabidopsis MAGIC genotypes
#'
#' @param snp_dir directory where SNP genotypes are stored. These files can be 
#' obtained with the `downloadMagicSnp()` function.
#'
#' @export
#'
#' @examples
#' downloadArabMagic('~/temp/magic_snp_1k/', tidy = FALSE)
#' tidyArabMagic('~/temp/magic_snp_1k/')
tidyArabMagic <- function(snp_dir){
	
	snp_dir <- path.expand(snp_dir)
	
	# Make a list of original file names
	ped_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.data"))
	alleles_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.alleles"))
	map_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.map"))
	
	# Check that all expected files exist
	if(!all(ped_files %in% list.files(snp_dir, pattern = ".data", full.names = T)) |
		 !all(alleles_files %in% list.files(snp_dir, pattern = ".alleles", full.names = T)) |
		 !all(map_files %in% list.files(snp_dir, pattern = ".map", full.names = T))){
		stop("Did not find all expected files. Try using the downloadArabMagic() function.")
	}
	
	
	##### Make file names reflect format #####
	file.rename(map_files, gsub("\\.map", "\\.happy\\.map", map_files))
	file.rename(ped_files, gsub("\\.data", "\\.happy\\.data", ped_files))
	file.rename(alleles_files, gsub("\\.alleles", "\\.happy\\.alleles", alleles_files))
	
	ped_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.happy.data"))
	alleles_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.happy.alleles"))
	map_files <- file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.happy.map"))

	
	##### Create plink PED and MAP files #####
	# Convert files to plink format
	plink_files <- lapply(1:5, function(i, snp_dir){
		.happy2plink(ped = file.path(snp_dir, paste0("chr", i, ".MAGIC.happy.data")),
								 map = file.path(snp_dir, paste0("chr", i, ".MAGIC.happy.map")),
								 alleles = file.path(snp_dir, paste0("chr", i, ".MAGIC.happy.alleles")),)
	}, snp_dir)
	
	# Write files
	lapply(1:5, function(i, plink_files, snp_dir){
		write.table(plink_files[[i]]$ped, 
								file.path(snp_dir, paste0("chr", i, ".MAGIC.plink.ped")),
								, quote = F, row.names = F, col.names = F, sep = "\t")
		write.table(plink_files[[i]]$map, 
								file.path(snp_dir, paste0("chr", i, ".MAGIC.plink.map")),
								, quote = F, row.names = F, col.names = F, sep = "\t")
	}, plink_files, snp_dir)
	

	##### Combine .happy.map files #####
	# Combine happy.map files
	lapply(map_files, read.table, header = T) %>%
		bind_rows() %>%
		write.table(file.path(snp_dir, "all_chr.MAGIC.happy.map"), 
								quote = F, row.names = F, sep = "\t")
	
	
	###### Combine .happy.alleles files ######
	# Read all files
	alleles_read <- lapply(alleles_files, readLines)
	
	# Extract the total number of SNPs
	tot_alleles <- sapply(alleles_read, function(x) strsplit(x[1], "\t")[[1]][2]) %>%
		as.numeric() %>% sum()
	
	# Get founder ids (second line of file)
	founder_ids <- alleles_read[[1]][2]
	
	# Combine the vectors excluding first two lines of each file
	alleles_read <- Reduce(c, lapply(alleles_read, function(x) x[c(-1, -2)]))
	
	# Write the file adding back the two lines with 
	# number of markers and accession names
	c(paste0("markers\t", tot_alleles,"\tstrains\t19"), 
		founder_ids, 
		alleles_read) %>%
		writeLines(file.path(snp_dir, "all_chr.MAGIC.happy.alleles"))
	
	
	###### Combine .happy.data files ######
	# Read ped files (excluding columns 1:6)
	# and combine by column
	ped_read <- lapply(ped_files, function(x) read.table(x)[-c(1:6)]) %>%
		bind_cols()
	
	# Output file adding the first 6 columns from one of the files
	cbind(read.table(ped_files[1])[,1:6], ped_read) %>%
		write.table(file.path(snp_dir, "all_chr.MAGIC.happy.data"), 
								quote = F, row.names = F, col.names = F)
	
	
	##### Combine plink.map files ####
	all_map <- lapply(file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.plink.map")), read.table, header = F) %>%
		bind_rows() 
	
	write.table(all_map, file.path(snp_dir, "all_chr.MAGIC.plink.map"), 
							quote = F, row.names = F, col.names = F, sep = "\t")
	
	
	##### Combine plink.ped files #####
	# Read ped files (excluding columns 1:6)
	# and combine by column
	ped_read <- lapply(plink_files, function(x) x[["ped"]][, -c(1:6)]) %>%
		bind_cols()
	
	# Output file adding the first 6 columns from one of the files
	cbind(plink_files[[1]]$ped[,1:6], ped_read) %>%
		write.table(file.path(snp_dir, "all_chr.MAGIC.plink.ped"), 
								quote = F, row.names = F, col.names = F)
	
	
	###### Convert Founder Genotypes to plink ######
	
	founders <- read.table(file.path(snp_dir, 'founders.genotypes.txt'), header = T)
	
	# Need to fix some SNP names which are wrong
	colnames(founders) <- gsub("\\.", "_", colnames(founders))
	
	# Make the .map file
	# Get SNPs which exist in the founders
	map_founders <- all_map[all_map$V2 %in% colnames(founders)[-1], ]
	
	# Make the ped file
	# First 6 columns of .ped file
	ped_founders <- data.frame(fid = founders$Founder, iid = founders$Founder,
														 dad = 0, mom = 0,
														 sex = 0, phen = 0)
	
	# Genotypes in .ped format
	ped_geno <- apply(founders[, map_founders$V2], 1, function(x){
		x[which(is.na(x))] <- "00"
		strsplit(x, "") %>% unlist()
	}) %>% t()
	
	ped_founders <- cbind(ped_founders, ped_geno)
	
	# Write files
	write.table(map_founders, file.path(snp_dir, "founders.map"), row.names = F, col.names = F, quote = F)
	write.table(ped_founders, file.path(snp_dir, "founders.ped"), row.names = F, col.names = F, quote = F)
	
}


#' Convert happy files to plink format
.happy2plink <- function(ped, map, alleles, snp_dir){
	# Read HAPPY .data file
	ped <- read.table(ped, stringsAsFactors = F)
	
	# Read HAPPY .map file
	map <- read.table(map, header = T, stringsAsFactors = F)
	
	# Read HAPPY .alleles file
	alleles <- readLines(alleles)
	alleles <- alleles[grep("marker\t", alleles)] %>% read.table(text = .) %>%
		rename(marker = V2) %>% select(marker)
	
	# If any alleles are not in common between the .alleles and .map files
	# Remove them from the ped file
	miss_allele <- which(!(alleles$marker %in% map$marker))
	
	if(length(miss_allele) > 0){
		ped <- ped[,-1*(6 + c(miss_allele*2-1, miss_allele*2))]
	}
	
	# Replace NA for '0' in .ped file
	ped[is.na(ped)] <- "0"
	
	# Merge alleles and map tables, this way 
	# exclude any missing markers from the map file
	map <- merge(alleles, map) %>%
		mutate(cm = 0) %>%
		select(chromosome, marker, cm, bp)
	
	# Return a list with .ped and .map files
	return(list(ped = ped, map = map))
}



