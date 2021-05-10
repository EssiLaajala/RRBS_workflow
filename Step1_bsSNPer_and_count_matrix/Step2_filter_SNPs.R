
# Set read counts to NA at detected SNPs, specifically for each individual
# This doesn't take long (approx 20 seconds per subject) and doesn't require too much memory (e.g. 16G was enough) so I recommend doing this interactively

##############################################################
# Choose these before you start
##############################################################

## You need count matrices (to filter the SNPs from) and the output vcf files of BS-SNPer

## Count matrices: Numbers of methylated reads and total numbers of reads at each CpG site at each sample. Dimensions (rows, columns) are number of CpGs, number of samples (e.g. 3 million rows and 173 columns). The row names are such as chr1:12345
#meth_filt <- read.table("path_to_file/meth_prefiltered.txt", header=T, sep="\t")
#total_filt <- read.table("path_to_file/total_prefiltered.txt", header=T, sep="\t")

## The vcf files:
## Set working directory to the output directory you used for BS-SNPer. It contains SNP out files named such as "SNP_Subject1.out" + temporary files. Check that all SNP out files end with "SNP finished", indicating that BS-SNPer wasn't interrupted.
## SNP out files are vcf files that contain 10 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, NA00001. Please check that this is the format you have.
colnames_vcf <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NA00001")
#files <- list.files(pattern="SNP")
#a <- sapply(files, function(z) strsplit(z, split="SNP_")[[1]][2])
#names(files) <- sapply(a, function(z) strsplit(z, split="\\.")[[1]][1])
## files is a character vector containing the names of SNP files ("SNP_Subject1.out" etc.) and names(files) should correspond to the column names of meth_filt and total_filt (sample names such as "Subject1" in this example).

## File names (paths) to write to:
#filename_total <- "/Results/total_prefiltered_noSNP.txt"
#filename_meth <- "/Results/meth_prefiltered_noSNP.txt"
#filename_SNPs <- "/Results/meth_prefiltered_noSNP.txt"

##################################################################
# Create a list of CpGs with a SNP marked "PASS" (in each subject)
# Remove SNPs from meth_filt and total_filt
##################################################################

# The working directory should still be the BS-SNPer output folder

snp_list <- vector("list", length(files))
for(sample in 1:length(files)){

	# Read a SNP.out file (for one sample)
	snp_out <- as.matrix(read.table(files[sample], header=F, sep="\t"))
	colnames(snp_out) <- colnames_vcf

	# Rownames of the snp_out matrix need to be of the format chr1:12345 (same as in total_filt and meth_filt)
	rownames(snp_out) <- paste(snp_out[,"CHROM"], snp_out[,"POS"], sep=":")
	rownames(snp_out) <- gsub(" ", "", rownames(snp_out)) # just to make sure

	# Pick our CpG sites with snps marked as "PASS"
	our_snps <- intersect(rownames(snp_out), rownames(total_filt))
	snp_out_our <- snp_out[our_snps,]
	snp_pass <- rownames(snp_out_our[snp_out_our[,"FILTER"]=="PASS",])

	snp_list[[sample]] <- snp_pass
	names(snp_list)[sample] <- names(files)[sample]

	total_filt[snp_pass,names(files)[sample]] <- NA
	meth_filt[snp_pass,names(files)[sample]] <- NA		
}

write.table(total_filt, filename_total, sep="\t", quote=F)
write.table(meth_filt, filename_meth, sep="\t", quote=F)
write.table(unlist(snp_list), filename_SNPs)


