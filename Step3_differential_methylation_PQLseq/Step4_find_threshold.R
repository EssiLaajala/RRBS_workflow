#############################################
# RADMeth spatial adjustment of P values
#############################################

# First you need to install Methpipe and the GNU scientific library. Then run this for each real and permuted covariate you performed the PQLseq analysis for:
# ~/apps/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 sorted_permuted_induced_labor.bed > sorted_adjusted_permuted_induced_labor.bed # Choose this path and file names!
# ~/apps/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 sorted_induced_labor.bed > sorted_adjusted_induced_labor.bed # Choose this path and file names!

# radmeth adjust runs within 30 seconds and requires less than 16G memory for a bed file of approx 2.5 million rows.

###################################
###################################
### Significance
###################################
###################################

# The code below was built under R version 4.0.4. 
# I recommend running this interactively, 16G memory is enough, runs in a few minutes.

###################################
# Choose these before you start
###################################

adjusted_bed_permuted <- # read.table("sorted_adjusted_permuted_induced_labor.bed", sep="\t", header=F) # Path to the adjusted bed file (output of radmeth adjust, produced above). This should contain P values and spatially adjusted P values for differential methylation associated with a permuted covariate (ideally, no differential methylation should be observed).
colnames(adjusted_bed_permuted) <- # c("chr","start","strand","context","pval","adj.pval", "FDR", "total.1","meth.1","total.0","meth.0") # This is the order of columns in the current output format of radmeth adjust (Methpipe version 3.4.3). Check the Methpipe manual to make sure that's the case.

adjusted_bed_real <- # read.table("sorted_adjusted_induced_labor.bed", sep="\t", header=F) # Path to the adjusted bed file (output of radmeth adjust, produced above) This should contain P values and spatially adjusted P values for differential methylation associated with a covariate of interest
colnames(adjusted_bed_real) <- # c("chr","start","strand","context","pval","adj.pval", "FDR", "total.1","meth.1","total.0","meth.0")

FDR_threshold <- # 0.05

DMC_path <- # "lonely_dmcs_induced_labor.bed" # Choose this output path/filename for the DMCs that have a Benjamini-Hochberg corrected P value (before spatial adjustment) under P_threshold. The output will be a bed file (some rows from adjusted_bed_real)

CpGs_within_candidate_DMRs_path <- # "CpGs_FDR5_induced_labor.bed" # Choose this output path/filename for the CpGs that are differentially methylated based on empirically FDR controlled spatially adjusted P values

#################################################
# Small modifications to adjusted_bed_permuted
#################################################

# Calculate the coverage corrected mean methylation difference ("dif")
dif <- (adjusted_bed_permuted[,"meth.1"]/adjusted_bed_permuted[,"total.1"])-(adjusted_bed_permuted[,"meth.0"]/adjusted_bed_permuted[,"total.0"])
dif <- round(dif, digits=3)
adjusted_bed_permuted <- cbind(adjusted_bed_permuted, dif)

## Check how many CpG sites would be differentially methylated based on default thresholds:
# false_discoveries <- adjusted_bed_permuted[adjusted_bed_permuted[,"FDR"]<0.05,]
# nrow(false_discoveries)
## How many would be significantly different based on BH-corrected raw P value (PQLseq P value, before spatial adjustment)?
# fdr_false = p.adjust(as.numeric(adjusted_bed_permuted[,"pval"]), method="BH")
# sum(fdr_false < 0.05)
## In my experience, PQLseq P values were fine, sum(fdr_false < 0.05) was almost always zero (sometimes 1 or 2). In contrast nrow(false_discoveries) was between 200 and 2500, indicating that the spatially adjusted P values were inflated.

# Correct format
adjusted_bed_permuted <- as.matrix(adjusted_bed_permuted)
adjusted_bed_permuted <- gsub(" ", "", adjusted_bed_permuted)
rownames(adjusted_bed_permuted) <- paste(adjusted_bed_permuted[,"chr"], adjusted_bed_permuted[,"start"], sep=":")

#################################################
# Small modifications to adjusted_bed_real
#################################################

# Calculate the coverage corrected mean methylation difference ("dif")
dif <- (adjusted_bed_real[,"meth.1"]/adjusted_bed_real[,"total.1"])-(adjusted_bed_real[,"meth.0"]/adjusted_bed_real[,"total.0"])
dif <- round(dif, digits=3)
adjusted_bed_real <- cbind(adjusted_bed_real, dif)

# Correct format
adjusted_bed_real <- as.matrix(adjusted_bed_real)
adjusted_bed_real <- gsub(" ", "", adjusted_bed_real)
rownames(adjusted_bed_real) <- paste(adjusted_bed_real[,"chr"], adjusted_bed_real[,"start"], sep=":")

adjusted_bed_real_original <- adjusted_bed_real

#####################################################
# Benjamini-Hochberg for raw P values
#####################################################

fdr_before_spatial_adjustment <- p.adjust(as.numeric(adjusted_bed_real[,"pval"]), method="BH")
if(sum(fdr_before_spatial_adjustment < FDR_threshold) > 0){
	lonely_dmcs <- adjusted_bed_real[fdr_before_spatial_adjustment < FDR_threshold,,drop=F]
	write.table(lonely_dmcs, DMC_path, sep="\t", quote=F)
}

#############################################################
# Empirical FDR control for the spatially adjusted P values
#############################################################

# Take only the rows for which PQLseq converged for both (for a fair comparison)
adjusted_bed_real <- adjusted_bed_real[intersect(rownames(adjusted_bed_real), rownames(adjusted_bed_permuted)),]
adjusted_bed_permuted <- adjusted_bed_permuted[intersect(rownames(adjusted_bed_real), rownames(adjusted_bed_permuted)),]
identical(rownames(adjusted_bed_real), rownames(adjusted_bed_permuted))
# TRUE

adjusted_bed_real <- cbind(adjusted_bed_real, rep("real", nrow(adjusted_bed_real)))
adjusted_bed_permuted <- cbind(adjusted_bed_permuted, rep("permuted", nrow(adjusted_bed_permuted)))

both <- rbind(adjusted_bed_permuted[,c(6,13)], adjusted_bed_real[,c(6,13)]) # Watch out! Hard-coded! Column 6 contains the spatially adjusted P values and column 13 contains categories "real" and "permuted".
both <- both[order(as.numeric(both[,"adj.pval"])),]
rownames(both) <- paste(rownames(both), both[,2], sep="_") # Watch out! Hard-coded!
# both is a matrix with the following format. It's sorted by the spatially adjusted P value and column 2 indicates whether the spatially adjusted P value was observed in the real or the permuted analysis (permuted meaning the covariate was permuted when PQLseq was run)
#                      adj.pval
#chr5:1594595_permuted "4.36207e-13" "permuted"
#chr5:1594605_permuted "4.36207e-13" "permuted"
#chr5:1594607_permuted "4.84612e-13" "permuted"

# Unless you expect more than 1 % of the CpGs to be differentially methylated, you can stop seeking for the threshold when you've run through top 1 %
stop_at <- round( nrow(both) * 0.01) # Watch out! Hard-coded! You can replace stop_at with any large enough value, e.g. 20000 (unless you expect to have tens of thousands of differentially methylated CpGs). You can also comment this out if you like, but in that case the following for-loop will take 100 minutes instead of one minute. 
both <- both[1:stop_at,]
fdr <- vector()
for(i in 1:nrow(both)) # This takes approx 1 minute if stop_at is approx 50000.
	fdr <- c(fdr, sum(both[1:i,2] == "permuted")/sum(both[1:i,2] == "real")) # Watch out! Hard-coded! The categories permuted and real need to be in column 2
idx <- tail(which(fdr < FDR_threshold), 1)
# The numbers of discoveries you'd get for the real and the permuted covariate with this threshold (both[idx,1]):
table(both[1:idx,2])

thres <- both[idx,1]
a <- adjusted_bed_real_original[as.numeric(adjusted_bed_real_original[,"adj.pval"]) <= thres,]
write.table(a, CpGs_within_candidate_DMRs_path, sep="\t", quote=F)

# I ran this script three times for each covariate of interest (three different permutations of the covariate). Chose the median empirically determined threshold for the spatially adjusted P values.

