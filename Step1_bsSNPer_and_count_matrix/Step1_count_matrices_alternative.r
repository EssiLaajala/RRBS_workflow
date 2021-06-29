# Works at least with R versions between 3.6.1 and 4.0.4.

# The easiest option would be to immediately filter out all measurements that are under 10 reads. However, that way we would lose information and PQLseq would converge less often. Here we keep measurements with 1-9 reads.
# This performs the same task as Step1_count_matrices.R. This version is uglier but needs less memory.

# Time and memory for 173 samples, on average approx. 5-6 million CpG sites detected in each sample with at least 1 read, 2-3 million with at least 10 reads
# Memory 40G, time 3 hours 15 minutes (this version)
# The other version (Step1_count_matrices.R) needed 80G memory, time 2 hours 45 minutes
# These requirements also depend on the minimum_coverage and minimum_nbr_samples you choose. I used 10 and 10. I do not recommend values under 3. 
# I tested these first interactively with just 3 samples. Five minutes with 16G was enough for that.

# Input: cov files produced by coverage2cytosine (where strand information has been merged). 
# The cov file has 6 columns: chr, start, end, methylation percentage, count methylated, count non-methylated

###################################################################################
# Choose these before you start
###################################################################################

# Set the working directory to the location of the cov files you want to include 
files <- # list.files(pattern="merged")
# files is a character vector of filenames such as "Subject1.CpG_report.merged_CpG_evidence.cov"

# Create a vector called sample_names. Must be in the same order and of the same length as files 
sample_names <- # sapply(files, function(z) strsplit(z, split=".CpG")[[1]][1])
# sample_names is a character vector of IDs such as Subject1 in this example. 

minimum_coverage <- # 10

# Choose a minimum number of samples the CpG site should be detected in (coverage >= minimum_coverage) to be included. You can save a lot of memory and time with some preliminary filtering, even though the actual coverage filtering might be done later. In my case, I know already that if some CpG site is detected in less than 10 samples, it has no way of fulfilling the actual coverage criteria (applied later, not as part of this script)
minimum_nbr_samples <- # 10 ##Choose this number! Must be an integer between 2 and length(files)

# Add filenames to write to (as characters). Full paths or paths from the working directory.
path_to_write_total <- # "Step1_count_matrices/Results/total_prefiltered.txt"
path_to_write_meth <- # "Step1_count_matrices/Results/meth_prefiltered.txt"

# Do you want to sort the count matrices by chromosome and location? (TRUE/FALSE). If this is set to FALSE, they will end up in aplhabetical order.
sort_the_matrix <- # TRUE
# If this is set to TRUE, don't worry about the warning "NAs introduced by coercion"

###################################################################################
# Create a list of coverage-filtered counts in each sample
###################################################################################
 
# lista is a list of matrices with two columns: methylated and total read counts. The rownames will be such as chr1:12345. The names of the elements of the list will be the sample names (such as "Subject1" in this example).
lista <- list()
for(i in 1:length(files)){
        a <- read.table(files[i], sep="\t", header=F)
        colnames(a) <- c("chr", "start", "end", "methylation percentage",  "count methylated", "count non-methylated") # Check these!
        # Replace count non-methylated with total counts
        total <- rowSums(a[,c("count methylated", "count non-methylated")])
        a <- cbind(a, total)
        a <- a[,setdiff(colnames(a), "count non-methylated")]
        # Remove sites with coverage < minimum_coverage
        filtered <- a[a[,"total"] >= minimum_coverage,]
        rownames(filtered) <- paste(filtered[,"chr"], filtered[,"start"], sep=":")
        # Mehtylated and total counts are enough information
        filtered <- filtered[,c("count methylated", "total")]
        lista[[i]] <- filtered
}
names(lista) <- sample_names

# Each matrix contains millions of rows and two columns, and the format should be this:
#> lista[[1]][1:4,]
#            count methylated total
#chr5:12155               14    45
#chr5:12177               73    89
#chr5:12182               87    90
#chr5:12188               84    90 

##########################################
# Frequently detected CpG sites
##########################################

all_sites <- vector()
for(i in 1:length(lista)) all_sites <- c(all_sites, rownames(lista[[i]]))
a <- duplicated(all_sites) # This is done first because all_sites is a huge vector and the next step (table) would take too long. Function "duplicated" was 57 times faster than the next step (function "table") in my toy data of 3 samples and 17 million sites
non_unique <- all_sites[a]
b <- table(non_unique)
# Pick CpG sites that are detected (coverage >= minimum_coverage) in at least some number of samples (the number you inserted in minimum_nbr_samples):
frequent_sites <- names(b[b>=(minimum_nbr_samples-1)]) #minimum_nbr_samples-1 because duplicated already removed one of each

############################################################################################################
# Create a list of counts in each sample, this time keeping values below minimum_coverage
############################################################################################################

# Use frequent_sites to get this done within reasonable time

lista <- list()
for(i in 1:length(files)){
        a <- read.table(files[i], sep="\t", header=F)
        colnames(a) <- c("chr", "start", "end", "methylation percentage",  "count methylated", "count non-methylated") # Check these!
        # Replace count non-methylated with total counts
        total <- rowSums(a[,c("count methylated", "count non-methylated")])
        a <- cbind(a, total)
        a <- a[,setdiff(colnames(a), "count non-methylated")]
        # Remove sites with coverage < minimum_coverage
		rownames(a) <- paste(a[,"chr"], a[,"start"], sep=":")
        filtered <- a[intersect(frequent_sites, rownames(a)),]
        # Mehtylated and total counts are enough information
        filtered <- filtered[,c("count methylated", "total")]
        lista[[i]] <- filtered
}
names(lista) <- sample_names

##########################################
# Count matrix
##########################################

if(sort_the_matrix == FALSE){
	# nrow is the number of unique sites (detected in at least the minimum number of samples, which is 10 in this example), ncol is the nbr of samples
	total <- matrix(nrow=length(frequent_sites), ncol=length(lista))
	meth <- matrix(nrow=length(frequent_sites), ncol=length(lista))
	rownames(total) <- frequent_sites
	rownames(meth) <- frequent_sites
	colnames(total) <- names(lista)
	colnames(meth) <- names(lista)

	for(i in names(lista)){
			total[rownames(lista[[i]]),i] <- lista[[i]][,"total"]
			meth[rownames(lista[[i]]),i] <- lista[[i]][,"count methylated"]
	}
}else{
	
	##########################################
	# Sort sites by chromosome and location
	##########################################	
		
	# This will sort numbers as numbers (e.g. chr2 will be before chr10 and chr1:123 will be before chr1:1112)
	# The CpG sites in vector frequent_sites must be of the format chr1:12345
		
	# Extract the chromosome and separate numbers 1-22 ("number chromosomes") from others ("letter chromosomes" such as X and Y)
	chr <- sapply(frequent_sites, function(z) strsplit(z, split=":")[[1]][1])
	a <- is.na(as.numeric(sapply(chr, function(z) strsplit(z, split="r")[[1]][2]))) # Chromosomes such as X,Y produce NAs so they will be TRUE
	# It will warn that NAs are introduced by coercion (that's supposed to happen)
	nbr <- frequent_sites[!a] # Locations in chromosomes 1-22
	letter <- frequent_sites[a] # locations in other chromosomes (X,Y etc.)
	chr_nbr <- chr[!a]
	chr_letter <- chr[a]

	### Separately sort chromosomes with numbers (1-22) and the others with letters

	## Numbers 1-22
	chr_nbr <- as.numeric(sapply(chr_nbr, function(z) strsplit(z, split="r")[[1]][2]))
	for_sorting_nbr <- cbind(chr_nbr, as.numeric(sapply(nbr, function(z) strsplit(z, split=":")[[1]][2])))
	colnames(for_sorting_nbr) <- c("chr", "loc")
	rownames(for_sorting_nbr) <- nbr
	# Sort the matrix
	for_sorting_nbr <- for_sorting_nbr[order(for_sorting_nbr[,"loc"]),]
	for_sorting_nbr <- for_sorting_nbr[order(for_sorting_nbr[,"chr"]),]

	## Other Chromosomes
	for_sorting_letter <- cbind(chr_letter, sapply(letter, function(z) strsplit(z, split=":")[[1]][2]))
	colnames(for_sorting_letter) <- c("chr", "loc")
	rownames(for_sorting_letter) <- letter
	# Sort the matrix
	for_sorting_letter <- for_sorting_letter[order(as.numeric(for_sorting_letter[,"loc"])),]
	for_sorting_letter <- for_sorting_letter[order(for_sorting_letter[,"chr"]),]

	for_sorting <- rbind(for_sorting_nbr, for_sorting_letter)

	##########################################
	# Count matrix
	##########################################

	# nrow is the number of unique sites (detected in at least the minimum number of samples, which is 10 in this example), ncol is the nbr of samples
	total <- matrix(nrow=nrow(for_sorting), ncol=length(lista))
	meth <- matrix(nrow=nrow(for_sorting), ncol=length(lista))
	rownames(total) <- rownames(for_sorting)
	rownames(meth) <- rownames(for_sorting)
	colnames(total) <- names(lista)
	colnames(meth) <- names(lista)

	for(i in names(lista)){
			total[rownames(lista[[i]]),i] <- lista[[i]][,"total"]
			meth[rownames(lista[[i]]),i] <- lista[[i]][,"count methylated"]
	}
}

############################################
# Write
############################################

write.table(total, path_to_write_total, quote=F, sep="\t")
write.table(meth, path_to_write_meth, quote=F, sep="\t")












