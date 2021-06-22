# Developed under R version 3.6.3.

source("PQLseq_EL.R") # The path to the modified pqlseq function
library(PQLseq)
library(Matrix) # for some reason this is not required by PQLseq but function nearPD is needed if the kinship matrix is singular (in this case it isn't so this package is not necessary)
# Versions PQLseq_1.1 and Matrix_1.2-18

#####################################
# Choose these before you start
#####################################

## This was written assuming there's a "covariates" matrix in addition to the "treatment" vector. I'm also assuming the covariates matrix includes some binary covariates. If it doesn't, skip the lines that check category-specific coverage and insert colnames(covariates) into covs_kept. Again, please modify this interactively. A couple gigabytes of memory is enough, and PQLseq runs within a minute or two if you set the rows_per_chunk to e.g. 100

## In this example, I'm reading in 30000 rows per job. This was run as an array job (92 independent jobs in this case, because my total number of rows is less than 92*30000 but more than 91*30000). The rows to read were identified based on the slurm task ID and rows_per_chunk (see the definition of start_idx below)
rows_per_chunk <- # 30000 # Choose this!

## At each CpG a binary covariate is excluded, unless there's a minimum number of samples with the minimum coverage in each category
min_coverage <- # 10 # Choose this!
min_nbr_samples <- # 3 # Choose this!

## Continuous covariates: character vector of column names corresponding to continuous covariates in the covariates matrix
continuous_covs <- # c("X.Intercept.", "PC1","PC2","birth.weight","gestational.weight.gain..mother", "year", "transformed_month", "BMI..mother", "height..mother", "age..mother") # Choose these!
## If you don't have any continuous covariates, insert the name of the intercept column: continuous_covs <- c("X.Intercept.")

# You might need to modify the "read the other inputs" section between rows 85 and 106 if your inputs have a different format (for example with respect to the presence/absence of column names)
# My "treatment" vector was saved as a text file in this format:
#Subject1 0
#Subject2 1
#Subject3 1
# etc.
# My "covariates" and "kinship" matrices were tab delimited text files that included column and row names

#####################################
# Read arguments 
#####################################

# I've provided an example shell script used together with this

args <- commandArgs(TRUE)
# The input arguments should be
# 1: path to meth matrix (nbr of methylated reads of each sample, each site)
# 2: path to total matrix (nbr of total reads of each site, each sample)
# 3: path to treatment vector
# 4: path to kinship matrix
# 5: output name
# 6: path to covariates matrix (covariates other than the "treatment")
# 7: slurm task ID # 1-92 in this case

#####################################
# Read the input count matrices
#####################################

start_idx = as.numeric(args[7])*rows_per_chunk-(rows_per_chunk-1)
nlines = length(readLines(args[1]))
if(start_idx + rows_per_chunk > nlines) rows_per_chunk = nlines - start_idx

cols_meth = as.vector(as.matrix(read.table(args[1], nrows=1)))
meth <- as.matrix(read.table(args[1], header=F,  skip=start_idx, nrows=rows_per_chunk, row.names=1))
colnames(meth) = cols_meth

cols_total = as.vector(as.matrix(read.table(args[2], nrows=1)))
total <- as.matrix(read.table(args[2], header=F,  skip=start_idx, nrows=rows_per_chunk, row.names=1))
colnames(total) = cols_total

# Remove CpGs with zero variation. They wouldn't throw an error but just a warning and then they'd end up as rows of NAs without rownames in results
total2 = total
total2[total2 == 0] = NA
prop = meth/total2
sds = apply(prop, 1, function(z) sd(z, na.rm=T))
keep = names(sds[sds != 0])
total = total[keep,]
meth = meth[keep,]

# Pseudo count transformation: Missing values should remain missing (PQLseq will leave them out) but other values get transformed to avoid methylation proportions that are exactly zero or exactly 1
# Run this part only once!
meth[total == 0] = NA
total[total == 0] = NA
meth <- meth + 1
total <- total + 2
meth[is.na(total)] = 0
total[is.na(total)] = 0
min_coverage <- min_coverage + 2
# My methylation proportions are now between 0.001897533 and 0.997830803

#####################################
# Read the other inputs
#####################################

treatment <- as.matrix(read.table(args[3], header=F))
#     V1        V2
#[1,] "Subject1" "0"
#[2,] "Subject2" "1"
#[3,] "Subject3" "1"
# etc.
a = as.numeric(treatment[,2]) # Careful! Hard-coded!
names(a) = treatment[,1] # Careful! Hard-coded!
treatment = a # Careful! Hard-coded!

kinship <- read.table(args[4], sep="\t", header=T, check.names=F)
covariates <- read.table(args[6], sep="\t", header=T)

meth_matches_total <- identical(rownames(meth), rownames(total))
covariates_match_total <- identical(rownames(covariates), colnames(total))
kinship_matches_meth <- identical(rownames(kinship), colnames(meth))
treatment_matches_meth <- identical(names(treatment), colnames(meth))

#####################################
# Run PQLseq
#####################################

# I used tryCatch, because without it PQLseq always ran into some error at some point (some exception among 2.7 million CpGs). If you don't need it and if you don't want to check the category-specific coverage for binary covariates, you don't need a loop at all. Just run:
# results <- pqlseq_EL(CountData=meth, LibSize=total, Phenotypes=treatment, Covariates=covariates, RelatednessMatrix=kinship, fit.maxiter=5000)
# (In that case you don't need a list structure either. The output of PQLseq is a matrix.)

cfit <- list() # List instead of matrix because some covariates won't be present for some CpG sites (different output dimensions)
if(meth_matches_total && covariates_match_total && kinship_matches_meth && treatment_matches_meth){
        for(i in 1:nrow(total)){
                write(rownames(total)[i],stderr())
                print(rownames(total)[i])

                ## Check which covariates can be taken into account for this CpG site (enough coverage in many enough samples)
                a <- total[i,] >= min_coverage 
				# a is a logical vector, whose length is the number of samples
                covs_kept <- continuous_covs
                binary_covs <- setdiff(colnames(covariates), covs_kept)
                for(bin_cov in binary_covs){
                        smallest <- min(table(a, covariates[,bin_cov])["TRUE",])
                        if(smallest >= min_nbr_samples) covs_kept = c(covs_kept, bin_cov)
                }

                ## Run PQLseq (just for covariates chosen above)
                tryCatch({
                cfit[[length(cfit)+1]] = pqlseq_EL(CountData=meth[i,,drop=F], LibSize=total[i,,drop=F], Phenotypes=treatment, Covariates=covariates[,covs_kept], RelatednessMatrix=kinship, fit.maxiter=5000)
				# If you use the original pqlseq instead of pqlseq_EL you need to specify that you are using a binomial mixed effect model (fit.model="BMM")!
                }, error = function(e) {
                        cat("ERROR :",conditionMessage(e), "\n")
                })
        }
}
names(cfit) = sapply(cfit, rownames)

###########################################################
# Create an output matrix (from the list created above)
###########################################################

nbr = unlist(sapply(cfit, ncol))
one_site_with_all_columns = names(nbr[nbr==max(nbr)])[1]
all_cols = colnames(cfit[[one_site_with_all_columns]])
results = matrix(nrow=length(cfit), ncol=length(all_cols))
rownames(results) = names(cfit)
colnames(results) = all_cols
for(i in 1:length(cfit)) results[rownames(cfit[[i]]), colnames(cfit[[i]])] = as.vector(as.matrix(cfit[[i]]))

###########
# Write
###########

write.table(results, args[5], sep="\t", quote=F)








