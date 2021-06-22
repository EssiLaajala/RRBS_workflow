# Developed under R version 3.6.3.

###################################################################
# Choose these before you start
###################################################################

## Load the SNP-filtered count matrices. The dimensions are number of CpG sites x number of samples
total <- # read.table("total_prefiltered_noSNP.txt", sep="\t", header=T)
meth <- # read.table("meth_prefiltered_noSNP.txt", sep="\t", header=T)

# Limit for the max coverage. If you insert 0.999, it means CpG sites with coverages above the 99.9th percentile will be removed (at each sample)
percentile <- # 0.999 # Choose this number!

## Choose a coverage threshold:
coverage_treshold <- # 10 # Choose this number!
## Choose the maximum number of samples that can have coverage below the coverage_treshold:
sample_treshold <- # 115

## Output file names (paths)
total_filename <- # "total_noSNP_filtered.txt"
meth_filename <- # "meth_noSNP_filtered.txt"

#####################################################################################
# Remove sites with a too great coverage (percentile 99.9, each sample separately)
#####################################################################################

# Currently, NA means either zero coverage or a SNP. Replace NA with 0.
total[is.na(total)] <- 0
meth[is.na(meth)] <- 0

# Run this only once!
limits <- vector()
for(i in 1:ncol(total)){
	sorted <- sort(total[,i], decreasing=T)
	sorted <- sorted[sorted>0]
	limit <- sorted[length(sorted)*(1-percentile)]
	meth[,i][total[,i]>limit] <- 0
	total[,i][total[,i]>limit] <- 0
	limits <- c(limits, limit)
}
names(limits) <- colnames(total)

#####################################################
# Filter
#####################################################

# Filter all rows where the total read count is low (either SNP or cov <10 or cov above 99.9th percentile) in more than sample_treshold number of samples

# Nbr of samples where each site (cytosine) was below coverage treshold
a <- rowSums(total < coverage_treshold)

# The cytosine will be removed if it was undetected in too many samples
remove_rows <- names(a)[a > sample_treshold]

total_filt <- total[setdiff(rownames(total),remove_rows),]
meth_filt <- meth[setdiff(rownames(meth),remove_rows),]

print(identical(rownames(total_filt), rownames(meth_filt))) # This should be true

write.table(total_filt, total_filename, sep="\t")
write.table(meth_filt, meth_filename, sep="\t")
