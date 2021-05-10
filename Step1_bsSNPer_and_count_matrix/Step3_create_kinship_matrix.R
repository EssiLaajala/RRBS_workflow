# Works at least with R versions between 3.6.3. and 4.0.4.

library(stringr)
# Developed under stringr version 1.4.0.

# This creates a kinship matrix based on correlations between SNP profiles, as described in the paper on the earlier version of PQLseq (MACAU). Lea AJ, Tung J, and Zhou X 2015: “A Flexible, Efficient Binomial Mixed Model for Identifying Differential DNA Methylation in Bisulfite Sequencing Data”, PLoS Genetics 11 (11): e1005650. The SNPs were detected using BS-SNPer.

# Digging out the genotype info from a vcf file is nasty work. Vcf files come in various formats. This will work for the output format of BS-SNPer, assuming it remains unchanged. Check column "FORMAT", which describes the format of the column "NA00001". In these vcf files genotype info ("GT") is in the very beginning of the :-delimited NA00001 columns. The FORMAT column says the format is GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR for all rows.

###################################################################
# Choose these before you start
###################################################################

## Choose an output directory
#output_directory <- "path/kinship_files/"

## Read a couple lines from a count matrix just to get the order of the columns (kinship matrix should be in the same order)
#total <- as.matrix(read.table("total_noSNP_filtered.txt", sep="\t",nrow=2))

## The vcf files:
## Set working directory to the output directory you used for BS-SNPer. It contains SNP out files named such as "SNP_Subject1.out" + temporary files. Check that all SNP out files end with "SNP finished", indicating that BS-SNPer wasn't interrupted.
## SNP out files are vcf files that contain 10 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NA00001. Please check that this is the format you have.
colnames_vcf <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA00001")
#files <- list.files(pattern="SNP")
#a <- sapply(files, function(z) strsplit(z, split="SNP_")[[1]][2])
#names(files) <- sapply(a, function(z) strsplit(z, split="\\.")[[1]][1])
## files is a character vector containing the names of SNP files ("SNP_Subject1.out" etc.) and names(files) should correspond to the column names of the count matrix total (sample names such as "Subject1" in this example).

## Choose the minor allele freq limit, minimum number of samples 
#minor_allele_freq_limit <- 9 # Choose this number! In this example, 9 is approx 5 % of the total number of samples (173)

##########################################################################
# Create a list of our CpGs with a SNP marked "PASS" (by bsSNPer)
##########################################################################

# Digging out the genotype info from a vcf file is nasty work. Vcf files come in various formats. This will work for the output format of BS-SNPer, assuming it remains unchanged. Check column "FORMAT", which describes the format of the column "NA00001". In these vcf files genotype info ("GT") is in the very beginning of the :-delimited NA00001 columns. .
# In the end you should have a list of matrices (one for each sample) like this (nbr is the number of reference alleles):
#            REF GT   nbr
#chr1:762787 "A" "CC" "0"
#chr1:763080 "G" "GT" "1"
#chr1:834404 "C" "CT" "1"
#chr1:842362 "C" "CT" "1"

# Ran in approx 20-30 seconds per sample
snp_list <- vector("list", length(files))
for(sample in 1:length(files)){

	# Read a SNP.out file (for one sample)
	snp_out <- as.matrix(read.table(files[sample], header=F, sep="\t"))
	colnames(snp_out) <- colnames_vcf

	# Rownames of the snp_out matrix need to be of the format chr1:12345 (same as in total_filt and meth_filt)
	rownames(snp_out) <- paste(snp_out[,"CHROM"], snp_out[,"POS"], sep=":")
	rownames(snp_out) <- gsub(" ", "", rownames(snp_out)) # just to make sure
	
	# Extract PASS SNPs:
	snp_out <- snp_out[snp_out[,"FILTER"]=="PASS",]
	
	# Add genotype
	# Specific only for the current output format of BS-SNPer! The FORMAT should be GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR for all rows (at least GT should be the first element). Please check this from the FORMAT column
	GT <- sapply(snp_out[,"NA00001"], function(z) strsplit(z, split=":")[[1]][1]) # Careful! Hard-coded!
	genotype <- cbind(snp_out[,c("REF", "ALT")], GT) # Careful! Hard-coded!
	gt2 <- apply(genotype, 1, function(z) gsub("1", z[2], z["GT"]))# Careful! Hard-coded!
	genotype[,"GT"] <- gt2 # Careful! Hard-coded!
	gt2 <- apply(genotype, 1, function(z) gsub("0", z[1], z["GT"]))# Careful! Hard-coded!
	GT <- gsub("\\/", "", gt2)# Careful! Hard-coded!
	snp_out <- cbind(snp_out, GT)
	alleles <- snp_out[,c("REF", "GT")] # Careful! Hard-coded!
	# Count the number of times the reference allele ("REF") appears in the genotype ("GT"). It's either 0 or 1 (because snp_out only contains sites where some deviation from the reference was detected in this sample)
	nbr <- apply(alleles, 1, function(z) str_count(z[2],z[1])) # Careful! Hard-coded!
	snp_list[[sample]] <- cbind(alleles, nbr)
	
	print(sample)
}
names(snp_list) <- files

################################
# Create the SNP matrix
################################

all_snps <- vector()
for(i in 1:length(snp_list)) all_snps <- c(all_snps, rownames(snp_list[[i]]))
all_snps <- unique(all_snps)

# snp_matrix will contain the numbers of reference alleles for all snps (that were detected in at least one sample). The dimensions (rows, columns) are length(all_snps), number of samples
# Please keep in mind, this assumes the number of reference alleles to be 2 unless a SNP was detected. Sometimes a SNP was not detected due to limited coverage. In other words, we are imputing missing values with 2 reference alleles.
snp_matrix <- vector()
for(i in 1:length(snp_list)){
	a <- rep(2, length(all_snps))
	names(a) <- all_snps
	a[rownames(snp_list[[i]])] <- as.numeric(snp_list[[i]][,"nbr"])
	snp_matrix <- cbind(snp_matrix, a)
}
colnames(snp_matrix) <- names(files)

######################################
# Filter for minor allele frequency
######################################

a <- rowSums(snp_matrix == 2)
limit <- ncol(snp_matrix)-minor_allele_freq_limit
keep <- names(a[a<=limit])
print(length(setdiff(names(a), keep)))
snp_matrix <- snp_matrix[keep,]

###########
# Write
###########

setwd(output_directory)
write.table(snp_matrix, "snp_matrix_X.txt", sep="\t", quote=F)

# Standardize it and create the kinship matrix
X <- snp_matrix
for(i in 1:ncol(X))
X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
K <- (t(X) %*% X)/nrow(X)
# K is a symmetric square matrix. I observed values 0.11 - 0.45 between unrelated individuals (and the diagonal should only contain 1s). Please keep in mind these are just correlations between SNP profiles in a set of SNPs that was detected in these samples. They cannot be used to determine, whether the individuals are e.g. 1st degree relatives.

# The trace norm should be 1. Check it:
a <- eigen(K, only.values = TRUE)$values
sum(abs(a))/ncol(K)

write.table(K, "kinship_SNP_profile_correlations.txt", sep="\t", quote=F)















