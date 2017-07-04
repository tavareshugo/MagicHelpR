#' Tidy Arabidopsis MAGIC genotypes
#'
#' @param snp_dir directory where SNP genotypes are stored. These files can be
#' obtained with the `downloadMagicSnp()` function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' downloadArabMagic('~/temp/magic_snp_1k/', tidy = FALSE)
#' tidyArabMagic('~/temp/magic_snp_1k/')
#' }
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

  # Adjust file name variables with the new names
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

  # Manually fix marker MN1_21908389
  ## This marker has an extra column in the probabilities for allele C
  i <- grep("MN1_21908389", alleles_read[[1]])
  alleles_read[[1]][i+2] <- "allele\tC\t0\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0\t0.083\t0.083\t0.083\t0.083\t0\t0\t0.083"

  # Extract the total number of SNPs
  tot_alleles <- sapply(alleles_read, function(x) strsplit(x[1], "\t")[[1]][2]) %>%
    as.numeric() %>% sum()

  # Get founder ids (second line of file)
  founder_ids <- alleles_read[[1]][2]

  # Combine the vectors excluding the first two lines of each file
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


  ###### Create long table with all MAGIC genotypes ######
  # This will be useful for easy loading into R

  # First get only one allele for each marker (as they are homozygous)
  magic_geno <- ped_read[, seq(1, ncol(ped_read), 2)]

  # Then get the marker names from the alleles file
  # (note that the .map file contains extra markers that are not present in
  # the .data file)
  magic_markers <- readLines(file.path(snp_dir, "all_chr.MAGIC.happy.alleles"))

  magic_markers <- magic_markers[grep("marker\t", magic_markers)] %>%
    strsplit("\t") %>%
    lapply(function(x) x[[2]]) %>%
    unlist()

  # Now add these as the names of the genotype table
  if(ncol(magic_geno) != length(magic_markers)) stop("Check some bug here!")
  names(magic_geno) <- magic_markers

  # Finally add the MAGIC IDs to the table
  magic_geno <- cbind(read.table(ped_files[1])[,1], magic_geno)
  names(magic_geno)[1] <- "magic"

  # Write the output in long format
  magic_geno %>%
    gather("marker", "allele", -magic) %>%
    write.table(file.path(snp_dir, "magic_genotypes.tsv"),
                row.names = F, quote = F, sep = "\t")


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
  founder_geno <- lapply(magic_markers, function(marker, alleles_read){

    # Get founder accession names (second line of the .happy.alleles file)
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

    # Issue warning if there are missing founder genotypes
    ## and there are...
    if(sum(grepl("0", alleles$allele)) > 0) warning(paste(marker, "missing some founder genotypes"))

    # Return a data.frame with marker name, accession and two columns with the allele
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
  names(ped_genos) <- magic_markers

  # Combine the two and write to file
  ## need to ensure only the markers with coordinates in the .map file are included
  bind_cols(ped_founders, ped_genos[all_map$V2]) %>%
    write.table(file.path(snp_dir, "founders.ped"), row.names = F, col.names = F, quote = F)

  # Write the .map file, which is the same as for the MAGIC lines
  write.table(all_map, file.path(snp_dir, "founders.map"), row.names = F, col.names = F, quote = F)

  # Write genotypes in long tabular format (potentially useful for other analysis)
  founder_geno <- bind_rows(founder_geno)

  founder_geno %>%
    select(acc, marker, allele1) %>%
    rename(allele = allele1, accession = acc) %>%
    mutate(allele = ifelse(allele == "0", NA, allele)) %>%
    write.table(., file.path(snp_dir, "founder_genotypes.tsv"), row.names = F, quote = F, sep = "\t")
}


#' Convert happy files to plink format
.happy2plink <- function(ped, map, alleles, snp_dir){
  # Read HAPPY .data file
  ped <- read.table(ped, stringsAsFactors = F)

  # Read HAPPY .map file
  map <- read.table(map, header = T, stringsAsFactors = F)

  # Read HAPPY .alleles file
  markers <- readLines(alleles)
  markers <- markers[grep("marker\t", markers)] %>%
    read.table(text = ., stringsAsFactors = F) %>%
    .$V2

  # If any markers are not in common between the .alleles and .map files
  # Remove them from the ped file
  # Some of these occur because they have no basepair coordinates
  miss_allele <- which(!(markers %in% map$marker))
  markers <- markers[-miss_allele]

  if(length(miss_allele) > 0){
    ped <- ped[,-1*(6 + c(miss_allele*2-1, miss_allele*2))]
  }

  # Replace NA for '0' in .ped file (this is the encoding used in PLINK)
  ped[is.na(ped)] <- "0"

  # Get the markers that were present in the alleles file,
  # which is the same order as the genotypes in the ped file.
  # This also removes missing markers from the map file.
  # I also add a new column with cM positions (set to 0, which is the default
  # missing in PLINK format).
  map <- map[match(markers, map$marker),] %>%
    mutate(cm = 0) %>%
    select(chromosome, marker, cm, bp)

  # Return a list with .ped and .map files in PLINK format
  return(list(ped = ped, map = map))
}



