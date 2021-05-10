# Developed under R version 3.6.3.

library(calibrate)

#################################################################
# Choose these before you start:
#################################################################

## Coverage filtered (and SNP-filtered) count matrices
## The dimensions(rows, columns) should be the number of SNPs (after coverage filtering), the number of samples 
#total_filt <- as.matrix(read.table("total_noSNP_filtered.txt", sep="\t", header=T, check.names=F))
#meth_filt <- as.matrix(read.table("meth_noSNP_filtered.txt", sep="\t", header=T, check.names=F))

## For further coverage filtering
#cov_limit = 10 # Choose this number!
## The minimum number of samples that should have minimum coverage cov_limit
#sample_nbr_limit = 87 Choose this number!

## Choose an output file name or path:
#output_filename <- "pca.txt"

#################################################################
# Some more coverage filtering (not to have to impute too much)
#################################################################

a = rowSums(total_filt >= cov_limit) 
keep = names(a[a >= sample_nbr_limit])

#####################################################
# Convert counts to methylation proportions
#####################################################

# If the total read count is 0, the value is actually missing
# (Also, we avoid inf values by removing zeros)
total_filt[total_filt==0] <- NA
prop_table <- meth_filt/total_filt

# Coverage filtering:
prop_table <- prop_table[keep,]

#####################################################
# Median imputation
#####################################################

# Impute NAs with the median of the cytosine

for(i in 1:nrow(prop_table)){
	z <- prop_table[i,]
	z[is.na(z)] <- median(z, na.rm=T)
	prop_table[i,] <-z
}

#####################################################
# Remove X and Y
#####################################################

remove <- grep("chrX",rownames(prop_table))
remove <- c(remove, grep("chrY",rownames(prop_table)))

if(length(remove)>0) prop_table <- prop_table[-remove,]

########
# PCA
########

pca <- prcomp(t(prop_table))
write.table(pca$x, output_filename, sep="\t", quote=F)


