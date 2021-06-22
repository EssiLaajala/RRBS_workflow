# Works at least with R versions between 3.6.1 and 4.0.4.

# The other (easier) option would be to immediately filter out all measurements that are under 10 reads. That way we would lose information and PQLseq would converge less often. Here we keep measurements with 1-9 reads.
# This script needs a lot of memory (70G was barely enough). I tested it first interactively with just 3 samples, 16G was enough for that
# Used cov files produced by coverage2cytosine (where strand information has been merged). 
# The cov file has 6 columns: chr, start, end, methylation percentage, count methylated, count non-methylated

###################################################################################
# Choose these before you start
###################################################################################

# Set the working directory to the location of the cov files you want to include 
files <- # list.files(pattern="merged")
# files is a character vector of filenames such as "Subject1.CpG_report.merged_CpG_evidence.cov"

# Create a vector called sample_names. Must be in the same order and of the same length as files 
sample_names # <- sapply(files, function(z) strsplit(z, split=".CpG")[[1]][1])
# sample_names is a character vector of IDs such as Subject1 in this example. 

# Choose a minimum number of samples the CpG site should be detected in (coverage >= 1) to be included. You can save a lot of memory and time with some preliminary filtering, even though the actual coverage filtering might be done later. In my case, I know already that if some CpG site is detected in less than 10 samples, it has no way of fulilling the actual coverage criteria (applied later, not as part of this script)
limit <- # 10 ##Choose this number! Must be an integer between 2 and length(files)

# Add filenames to write to (as characters). Full paths or paths from the working directory.
path_to_write_total <- # "Step1_count_matrices/Results/total_prefiltered.txt"
path_to_write_meth <- # "Step1_count_matrices/Results/meth_prefiltered.txt"

# Do you want to sort the count matrices by chromosome and location? (TRUE/FALSE). If this is set to FALSE, they will end up in aplhabetical order.
sort_the_matrix <- # TRUE
# If this is set to TRUE, don't worry about the warning "NAs introduced by coercion"

###################################################################################
# Make a list of (modified) cov files
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
        rownames(a) <- paste(a[,"chr"], a[,"start"], sep=":")
        # Mehtylated and total counts are enough information
        a <- a[,c("count methylated", "total")]
        lista[[i]] <- a
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
# Filter by frequency of detection
##########################################

all_sites <- vector()
for(i in 1:length(lista)) all_sites <- c(all_sites, rownames(lista[[i]]))
a <- duplicated(all_sites) # This is done first because all_sites is a huge vector and the next step (table) would take too long. Function "duplicated" was 57 times faster than the next step (function "table") in my toy data of 3 samples and 17 million sites
non_unique <- all_sites[a]
b <- table(non_unique)
# Pick CpG sites that are detected (coverage >= 1) in at least some number of samples (the number you inserted in limit):
frequent_sites <- names(b[b>=(limit-1)]) #limit-1 because duplicated already removed one of each

# Filter each matrix (include only frequent_sites):
for(i in 1:length(lista)) lista[[i]] <- lista[[i]][intersect(rownames(lista[[i]]), frequent_sites),]

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

## I did some more preliminary filtering to save a bit of time/memory at later steps (minimum 10 samples with coverage >=10).
#a <- rowSums(total>=10, na.rm=T)
#keep <- names(a[a>=10])
#total <- total[keep,]
#meth <- meth[keep,]

write.table(total, path_to_write_total, quote=F, sep="\t")
write.table(meth, path_to_write_meth, quote=F, sep="\t")












