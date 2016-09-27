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
	
	
	##### Rename files to reflect format #####
	# Rename to include "happy" in the file  format
	file.rename(map_files, gsub("\\.map", "\\.happy\\.map", map_files))
	file.rename(ped_files, gsub("\\.data", "\\.happy\\.data", ped_files))
	file.rename(alleles_files, gsub("\\.alleles", "\\.happy\\.alleles", alleles_files))
	
	# Make variables with these file names
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
	lapply(map_files, read.table, header = T, stringsAsFactors = F) %>%
		bind_rows() %>%
		write.table(file.path(snp_dir, "all_chr.MAGIC.happy.map"), 
								quote = F, row.names = F, sep = "\t")
	
	
	###### Combine .happy.alleles files ######
	# Read all files
	alleles_read <- lapply(alleles_files, readLines)
	
	# Fix marker MN1_21908389
	## This has an extra column in the probabilities for allele C
	i <- grep("MN1_21908389", alleles_read[[1]])
	alleles_read[[1]][i+2] <- "allele\tC\t0\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0.083\t0.083\t0\t0\t0.083"
	
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
	ped_read <- lapply(ped_files, function(x) read.table(x)[, -c(1:6)]) %>%
		bind_cols()
	
	# Output file adding the first 6 columns from one of the files
	cbind(read.table(ped_files[1])[,1:6], ped_read) %>%
		write.table(file.path(snp_dir, "all_chr.MAGIC.happy.data"), 
								quote = F, row.names = F, col.names = F)
	
	
	##### Combine plink.map files ####
	all_map <- lapply(file.path(snp_dir, paste0("chr", 1:5, ".MAGIC.plink.map")), 
										read.table, header = F, stringsAsFactors = F) %>%
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
	# Founder genomes are "reconstructed" from the alleles files
	alleles_read <- readLines(file.path(snp_dir, "all_chr.MAGIC.happy.alleles"))
	
	# For each marker present in the MAGIC lines, extract the allele of each accession
	## Some are missing!
	founder_geno <- lapply(all_map$V2, function(marker, alleles_read){
		
		# Get founder names
		founder_id <- unlist(strsplit(alleles_read[2], "\t"))[-1]
		
		# Read the genotype probabilities
		alleles <- alleles_read[grep(marker, alleles_read) + 2:3]
		
		# Stop if the marker was not found
		if(length(alleles) == 0) stop(paste("Missing marker", marker))
		
		# Convert to data.frame (using read.table for this)
		alleles <- read.table(text = alleles, fill = NA, stringsAsFactors = F)[, -1]
		
		# Stop if the number of columns is not as expected
		if(ncol(alleles) > 20) stop(paste(marker, "Too many probabilities present in file"))
		
		# Name the columns with each accession ID
		names(alleles) <- c("allele", founder_id)
		
		# Get the allele from each accession
		alleles <- alleles %>%
			gather("accession", "prob", Bur:Zu) %>%
			group_by(accession) %>%
			summarise(allele = ifelse(sum(prob) == 0, "0", allele[which.max(prob)]))
		
		# Issue warning if the number of alleles is not as expected
		if(sum(grepl("o", alleles$allele)) > 19) warning(paste(marker, "too many founder alleles!"))
		if(sum(grepl("o", alleles$allele)) < 19) warning(paste(marker, "missing some founder genotypes"))
		
		# Returna  data.frame with marker name, accession and two columns with the allele
		return(data.frame(marker = marker, acc = alleles$accession, 
											allele1 = alleles$allele, 
											allele2 = alleles$allele))
		
	}, alleles_read)
	
	# Make the .ped file - first 6 columns
	ped_founders <- data.frame(fid = founder_geno[[1]]$acc,
						 iid = founder_geno[[1]]$acc,
						 dad = 0, mom = 0,
						 sex = 0, phen = 0)
	
	# Get the genotype columns from the above list
	ped_genos <- lapply(founder_geno, function(x) select(x, allele1, allele2))
	
	# Combine the two and write to file
	cbind(ped_founders, ped_genos) %>%
		write.table(file.path(snp_dir, "founders.ped"), row.names = F, col.names = F, quote = F)
	
	# Write the .map file, which is the same as for the MAGIC lines
	write.table(all_map, file.path(snp_dir, "founders.map"), row.names = F, col.names = F, quote = F)
	
	# Write genotypes in long tabular format (potentially useful for other analysis)
	founder_geno <- bind_rows(founder_geno)
	founder_geno <- merge(all_map, founder_geno, by.x = "V2", by.y = "marker")
	names(founder_geno) <- c("marker", "chr", "cm", "pos", "acc", "allele")
	
	founder_geno %>%
		select(marker, chr, pos, acc, allele) %>%
		mutate(allele = ifelse(allele == "0", NA, allele)) %>%
		write.table(., file.path(snp_dir, "founders.tsv"), row.names = F, quote = F, sep = "\t")
}


#' Convert happy files to plink format
.happy2plink <- function(ped, map, alleles, snp_dir){
	# Read HAPPY .data file
	ped <- read.table(ped, stringsAsFactors = F)
	
	# Read HAPPY .map file
	map <- read.table(map, header = T, stringsAsFactors = F)
	
	# Read HAPPY .alleles file
	markers <- readLines(alleles)
	markers <- markers[grep("marker\t", markers)] %>% read.table(text = ., stringsAsFactors = F) %>%
		.$V2
	
	# If any alleles are not in common between the .alleles and .map files
	# Remove them from the ped file
	# Some of these occur because they have no basepair coordinates
	miss_allele <- which(!(markers %in% map$marker))
	markers <- markers[-miss_allele]
	
	if(length(miss_allele) > 0){
		ped <- ped[,-1*(6 + c(miss_allele*2-1, miss_allele*2))]
	}
	
	# Replace NA for '0' in .ped file
	ped[is.na(ped)] <- "0"
	
	# Get the markers that were present in the alleles file
	# which is the right order for the ped file. This also
	# excludes any missing markers from the map file.
	map <- map[match(markers, map$marker),] %>%
		mutate(cm = 0) %>%
		select(chromosome, marker, cm, bp)
	
	# Return a list with .ped and .map files
	return(list(ped = ped, map = map))
}



