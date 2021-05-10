
# This needed almost 30 GB of memory, since I have quite a large number of columns in the PQLseq output (and 2.7 million rows)

########################################################################
# Choose these before you start
########################################################################

# Metadata
cov <- as.matrix(read.table("../Data/covariates_173_samples.txt", sep="\t", header=T)) # Choose this path!
sex <- read.table("../Data/treatment_sex_173.txt", header=F) # Choose this path! Choose also a correct name for the variable (sex in this example), which was given as input vector "Phenotypes" for PQLseq
#     V1        V2
#[1,] "Subject1" "0"
#[2,] "Subject2" "1"
#[3,] "Subject3" "1"
# etc.
a = sex[,2] # Careful! Hard-coded!
names(a) = sex[,1] # Careful! Hard-coded!
sex = a
if(identical(names(sex), rownames(cov))) cov = cbind(cov, sex)

# Choose the location of PQLseq output files. This assumes the working directory only contains those files and only from one analysis.
setwd("/scratch/project_2002156/essilaaj/DIPP_cord_blood_RRBS/Vaihe4_PQLseq_173/Vaihe2_diff_meth/Tulokset/results_permuted_induced_labor_trick_3")

###################################################################
# Combine outputs, pick converged
###################################################################

# Combine PQLseq output files into one data frame of approx 2.6 million sites (this takes a few minutes)

files <- list.files()
output <- vector()
for(i in files){
        out1 <- read.table(i, header=F, sep="\t", skip=1)
        output <- rbind(output, out1)
}
output <- as.matrix(output)

output <- output[output[,1]!="NULL",]
rownames(output) <- output[,1]
output <- output[,-1]

names_out <- as.vector(as.matrix(read.table(i, header=F, sep="\t", nrows=1)))
colnames(output) <- names_out

output_converged <- output[output[,"converged"] == " 1",]

# Modify this!
colnames(output_converged) <- gsub("Phenotypes", "sex", colnames(output_converged)) # Careful! Hard-coded! Choose this!

##########################################################################
# Load count matrices and pick CpGs for which PQLseq converged
##########################################################################

# Here I also make sure the columns correspond to the rows of the design matrix (although that has been double checked in earlier steps)

meth <- as.matrix(read.table("../../../../Vaihe2_filter/Tulokset/meth_noSNP_filt_173.txt", sep="\t", header=T, check.names=F))
meth <- meth[rownames(output_converged),rownames(cov)]
total <- as.matrix(read.table("../../../../Vaihe2_filter/Tulokset/total_noSNP_filt_173.txt", sep="\t", header=T, check.names=F))
total <- total[rownames(output_converged),rownames(cov)]

###########################################################################
# Sort total and meth
###########################################################################

chr <- sapply(rownames(total), function(z) strsplit(z, split=":")[[1]][1])
coord <- sapply(rownames(total), function(z) strsplit(z, split=":")[[1]][2])

a <- is.na(as.numeric(sapply(chr, function(z) strsplit(z, split="r")[[1]][2]))) # Chromosomes such as X,Y produce NAs so they will be TRUE
# It will warn that NAs are introduced by coercion (that's supposed to happen)
nbr <- rownames(total)[!a] # Locations in chromosomes 1-22
letter <- rownames(total)[a] # locations in other chromosomes (X,Y etc.)
chr_nbr <- chr[!a]
chr_letter <- chr[a]

### Separately sort chromosomes with numbers (1-22) and the others with letters ###

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

# for_sorting <- rbind(for_sorting_nbr, for_sorting_letter)
# I only kept chromosomes 1 - 22
total <- total[rownames(for_sorting_nbr),] 
meth <- meth[rownames(total),]

############################################################
############################################################
#### No need to re-run the lines above for each covariate
############################################################
############################################################

#####################################################
# Choose these before you start!
#####################################################

# The covariate you want to build the bed file for (also coverage filtering will be done for that one)
covariate <- "induced.labor" # Choose this!
output_path <- "sorted_PQLseq_induced_labor_3.bed" # Choose this! 

coverage_limit <- 10 # Choose this number!
continuous_limit <- 0.33 # Choose this number! Minimum proportion of samples with non-missing values
per_level_limit <- 0.33 # # Choose this number! Minimum proportion of samples with non-missing values
per_level_limit_second <- 5 # # Minimum number of samples with non-missing values for binary covariates (because some binary variables might have only e.g. 8 samples representing one category, 33% is not enough)

# Choose these!
continuous_covs = c("birth.weight", "gestational.weight.gain..mother", "year", "transformed_month", "BMI..mother", "height..mother", "age..mother")
binary_covs = c("class", "epidural.anaesthetic", "apgar_low", "C.section","smoking.during.pregnancy", "nbr.of.earlier.miscarriages", "insulin.treated.diabetes..mother", "induced.labor", "sex")

#####################################################
# Coverage-filtering for the covariate of interest
#####################################################

output_converged <- output_converged[rownames(total),] # This way round because meth and total have been sorted

if(covariate %in% continuous_covs){
        lim <- round(continuous_limit*ncol(total), digits=1)
        keep <- rownames(total[rowSums(total >= coverage_limit) > lim,])
}

if(covariate %in% binary_covs){
        a <- total[,cov[,covariate] == 1]
        b <- total[,cov[,covariate] == 0]
        lim_a <- round(per_level_limit*ncol(a), digits=0)
        lim_b <- round(per_level_limit*ncol(b), digits=0)
        keep1 <- rownames(a[rowSums(a >= coverage_limit) > lim_a,])
        keep2 <- rownames(b[rowSums(b >= coverage_limit) > lim_b,])
        keep3 <- rownames(a[rowSums(a >= coverage_limit) > per_level_limit_second,])
        keep4 <- rownames(b[rowSums(b >= coverage_limit) > per_level_limit_second,])
        keep <- intersect(keep1, keep2)
        keep <- intersect(keep, keep3)
        keep <- intersect(keep, keep4)
}
print(identical(intersect(keep1, keep2), keep)) # This should be true for most covariates (if lim_a and lim_b are at least 5)

##############################################################
# BED format (needed for the spatial adjustment of P values)
##############################################################

# The bed file should contain 9 columns: chromosome, location, strand, context, P value, total_group1, meth_group1, total_group0, meth_group0
#chr1    10497   +       CpG     0.192358        9553    7629    9841    7740
#chr1    10525   +       CpG     0.903816        15835   14900   14041   13230
#chr1    10542   +       CpG     0.768604        15854   15103   14053   13311
#chr1    10563   +       CpG     0.210359        15853   14265   14057   12705

bed <- output_converged[keep,]

chr <- sapply(rownames(bed), function(z) strsplit(z, split=":")[[1]][1])
coord <- sapply(rownames(bed), function(z) strsplit(z, split=":")[[1]][2])

pvals <- bed[,paste0("pvalue_", covariate)]
bed <- data.frame(chr=chr, coord=coord, strand=rep("+", nrow(bed)), context=rep("CpG",nrow(bed)), pValue=pvals)

bed <- bed[keep,]

# put bed in the same format with RADMeth output
# radmeth adjust didn't work without these columns although it doesn't use them at all (I tested putting 0 to all these columns but the adjusted p-values were the same)
if(covariate %in% binary_covs){
        tot_1 <- rowSums(total[rownames(bed), cov[,covariate] == 1], na.rm=T)
        meth_1 <- rowSums(meth[rownames(bed), cov[,covariate] == 1])
        tot_0 <- rowSums(total[rownames(bed), cov[,covariate] == 0], na.rm=T)
        meth_0 <- rowSums(meth[rownames(bed), cov[,covariate] == 0])
        bed <- cbind(bed, tot_1, meth_1, tot_0, meth_0)
}


if(covariate %in% continuous_covs){
        low_lim <- as.numeric(summary(cov[,covariate])["1st Qu."])
        high_lim <- as.numeric(summary(cov[,covariate])["3rd Qu."])
        tot_1 <- rowSums(total[rownames(bed), cov[,covariate] >= high_lim], na.rm=T)
        meth_1 <- rowSums(meth[rownames(bed), cov[,covariate] >= high_lim])
        tot_0 <- rowSums(total[rownames(bed), cov[,covariate] <= low_lim], na.rm=T)
        meth_0 <- rowSums(meth[rownames(bed), cov[,covariate] <= low_lim])
        bed <- cbind(bed, tot_1, meth_1, tot_0, meth_0)
}

#######################
# write
#######################

write.table(bed, output_path, sep="\t", row.names=F, col.names=F, quote=F) # Check

