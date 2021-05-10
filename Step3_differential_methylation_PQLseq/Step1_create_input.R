# This is a documentation of the input files in my project. They are not publicly available, because they contain sensitive data, but the format is described here. Two input files for PQLseq are created: The "treatment vector" (design vector for one covariate) and the "covariates" matrix (design matrix for all other covariates except the treatment vector). The other three input files were already created in the coverage filtering step (methylated and total read count matrices) and the BS-SNPer step (kinship matrix)

# It doesn't matter which covariate is in the treatment vector if you use my modification of PQLseq, which outputs all coefficients, standard errors and P values. Just, you might want to choose some covariate with enough samples in each category if it's a binary covariate. If you only have a handful of samples in a category, PQLseq will fail at CpG sites, where they all have zero coverage. I double checked that PQLseq doesn't treat the "treatment" any differently from the other covariates. In contrast, the earlier version of PQLseq, which is called MACAU, puts a different prior on the treatment. (I wouldn't recommend using MACAU anyway, because it doesn't give you MCMC convergence diagnostics and at least in my data it was very unstable).

# Read in a couple lines from one of the count matrices just to make sure the design matrix will be in the same order
total_173 = read.table("total_noSNP_filtered.txt", sep="\t", header=T, nrows=2)
clinical_data = clinical_data[colnames(total_173),]

# Before this transformation, the month of birth was encoded as numbers 1-12
transformed_month <- cos((2*pi)/12*as.numeric(clinical_data[,"month"]))
clinical_data <- cbind(clinical_data, transformed_month)

model = c("class", "birth weight", "epidural anaesthetic", "sex", "apgar_low", "C-section", "gestational weight gain, mother", "smoking during pregnancy", "year", "transformed_month", "BMI, mother", "height, mother", "nbr of earlier miscarriages", "age, mother", "insulin-treated diabetes, mother", "induced labor")

################################################
# treatment vector
################################################

treatment = clinical_data[,"sex"]
model = setdiff(model, "sex")

################################################
# Covariates matrix
################################################

pca = as.matrix(read.table("pca.txt", sep="\t", header=T))
pca = pca[rownames(clinical_data),]

batch <- as.factor(clinical_data[,"Batch"])
batch <- model.matrix(~batch, batch)

covariates <- cbind(batch, pca[,1:2], clinical_data[,model])

##########################################################
# Impute missing values with medians of each covariate
##########################################################

# Check that these medians make sense for(i in 1:ncol(covariates)) {print(colnames(covariates)[i]); print(median(as.numeric(covariates[,i]), na.rm=T))}
# Most covariates don't have missing values anyway (max 4)

for(i in 1:ncol(covariates)) covariates[is.na(covariates[,i]),i] <- median(as.numeric(covariates[,i]), na.rm=T)

######################################
# z-transform continuous covariates
######################################

# PQLseq p-values are invariant to z-transformation (or any other linear scaling) but just to improve interpretability and convergence:

continuous = c("PC1", "PC2", "birth weight", "gestational weight gain, mother", "year", "transformed_month", "BMI, mother", "height, mother", "age, mother")

for(i in continuous){ 
	a = as.numeric(covariates[,i])
	covariates[,i] = (a-mean(a))/sd(a)
}

##########################################################
# Write tables
##########################################################

# Check that it's full-rank
mode(covariates) = "numeric"
qr(covariates)$rank == ncol(covariates)

write.table(covariates, "/Results/covariates_173_samples.txt", sep="\t", quote=F)
write.table(treatment, "/Results/treatment_173_sex.txt", quote=F, col.names=F)



