# This converts from the input format of PQLseq to the input format of RADMeth.

###########################################################################
# Design matrix for RADMeth: Continuous covariates transformed to binary
###########################################################################

## Choose these:

# Loaded the design matrix that was used for PQLseq (where continuous covariates are continuous):
covs <- # read.table("sex_and_other_covariates_173_samples.txt", sep="\t", header=T)

covs_continuous <- # covs[,c("PC1", "PC2","birth.weight","gestational.weight.gain..mother", "year", "transformed_month","BMI..mother","height..mother","age..mother")]
output_filename_design <- # "../Results/design_173_samples.txt"

####################################################################################################

covs_binary <- covs[,setdiff(colnames(covs), colnames(covs_continuous))]

new_continuous <- vector()
for(i in colnames(covs_continuous)){
        low_lim <- summary(covs_continuous[,i])["1st Qu."]
        high_lim <- summary(covs_continuous[,i])["3rd Qu."]
        small <- cbind(as.numeric(covs_continuous[,i] < low_lim), as.numeric(covs_continuous[,i] > high_lim))
        colnames(small) <- c(paste0("low_", i), paste0("high_", i))
        new_continuous <- cbind(new_continuous, small)
}

covs_rad <- cbind(covs_binary, new_continuous)

# Make sure it's still full rank
covs_rad = as.matrix(covs_rad)
qr(covs_rad)$rank
dim(covs_rad)

write.table(covs_rad, output_filename_design, sep="\t", quote=F)

###############################################################################################
# Proportion table chunks
###############################################################################################

# PQLseq works with two count matrices: methylated and total reads, whereas the input of RADMeth is one matrix, where columns total and methylated are next to each other, as described in the MethPipe manual section 3.3.2: http://smithlabresearch.org/downloads/methpipe-manual.pdf
# To parallelize RADMeth, I divided the input to 56 chunks of 50000 CpGs

## Choose these:

total <- # read.table("Step2_coverage_filtering/Results/total_noSNP_filtered.txt", sep="\t", header=T) # The total read counts for each CpG in each sample (dimensions: 2.7 million CpG sites, 173 samples)
meth <- # read.table("Step2_coverage_filtering/Results/meth_noSNP_filtered.txt", sep="\t", header=T) # The methylated read counts for each CpG in each sample (dimensions: 2.7 million CpG sites, 173 samples)
# To parallelize RADMeth, I divided the input to 56 chunks of 50000 CpGs
nbr_cpgs <- # 50000 
number_of_files <- # 56

output_filename <- # "../Results/proportion_chunks/prop_table_" # The final filename of each chunk will be paste0(filename,i,".txt"), e.g. prop_table_1.txt

###############################################################################################

total = total[,rownames(covs_rad)] # They are already in the same order but just to make sure
meth = meth[,rownames(covs_rad)] # They are already in the same order but just to make sure

# Missing values should remain missing but other values get transformed to avoid exactly zeros and exactly 1s
# If you skip this pseudo-count transformation, RADMeth's p-value inflation will be much worse
meth[total == 0] = NA
total[total == 0] = NA
meth <- meth + 1
total <- total + 2
meth[is.na(total)] = 0
total[is.na(total)] = 0
# Checked that methylation proportions are now between 0.001897533 and 0.997830803

prop_matrix = vector()
for(i in rownames(covs_rad)) prop_matrix = cbind(prop_matrix, total[,i], meth[,i])
rownames(prop_matrix) = paste0(rownames(total), ":+:CpG")

# Divide it into chunks and write:
from <- 1
to <- nbr_cpgs
for(i in 1:number_of_files){
write.table(prop_matrix[from:to,], paste0(output_filename,i,".txt"), sep="\t", quote=F, col.names=F)
from <- from+nbr_cpgs
to <- to+nbr_cpgs
if(to > nrow(prop_matrix)) to = nrow(prop_matrix)
}

# Added the column names later, e.g. cat colnames173.txt prop_table_1.txt > prop_table_1





