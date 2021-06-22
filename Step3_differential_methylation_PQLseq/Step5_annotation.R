library(genomation)
library(GenomicRanges)
# These can be installed from Bioconductor

# Developed under R version 4.0.4, genomation 1.22.0, GenomicRanges 1.42.0

# I did a bit of manual work in the middle (with UCSC genome browser, marked STOP)

#############################################################
# Choose these before you start
#############################################################

sigf <- # read.table("CpGs_FDR5_sex.bed", sep="\t", header=T) # Path to CpGs that were part of candidate differentially methylated regions based on empirically FDR controlled spatially adjusted P values. This is a tab-delimited data frame (bed file) with rownames such as chr1:1099982 and the following columns:
# chr	start	strand	context	pval	adj.pval	FDR	total.1	meth.1	total.0	meth.0	dif
lonely_dmcs <- # read.table("lonely_dmcs_sex.bed", sep="\t", header=T) # Path to DMCs (based on Benjamini-Hochberg corrected P values before spatial adjustment

path_to_annotation_file <- # "../hg19_knownGenes.bed" # Choose this!

output_path <- # "../Results/sigf_sex_annotated.txt" # Path to the output file. 

#############################################################
# Add the type of significance to the last column
#############################################################

sigf <- cbind(sigf, rep("significant_after_spatial_adjustment", nrow(sigf)))
colnames(sigf)[ncol(sigf)] <- "type"
sigf[intersect(rownames(lonely_dmcs), rownames(sigf)), "type"] <- "significant_with_both_criteria"

lonely_dmcs <- lonely_dmcs[setdiff(rownames(lonely_dmcs), rownames(sigf)),]
lonely_dmcs <- cbind(lonely_dmcs, rep("significant_before_spatial_adjustment", nrow(lonely_dmcs)))
colnames(lonely_dmcs)[ncol(lonely_dmcs)] = "type"
sigf = rbind(sigf, lonely_dmcs)

#############################################################
# Turn it into a GRanges object
#############################################################

end <- sigf[,"start"]
sigf <- cbind(sigf, end)
sigf.gr <- as(sigf, "GRanges")

##########################################
# Annotation to nearest gene by TSS
##########################################

gene.obj <- readTranscriptFeatures(path_to_annotation_file)
sigf.ann <- annotateWithGeneParts(sigf.gr, gene.obj)
ann <- getAssociationWithTSS(sigf.ann)
ann <- as.matrix(ann)
## ann is an annotation matrix with 4 columns:
#             target.row dist.to.feature feature.name feature.strand
#chr1:1099982 "   1"     "  -2502"       "uc001acw.2" "+"           
#chr1:1099986 "   2"     "  -2498"       "uc001acw.2" "+"      

###########################
###########################
# STOP
###########################
###########################

# Translate uc codes to gene names

## Used UCSC table browser (in Tools) to get gene symbols for the uc-codes in unique(ann[,"feature.name"]). Saved them to a file called uc_to_gene_names_sex.txt.
# Group: Genes and gene predictions
# Assembly: hg19
# Identifiers: paste list
# Track: UCSC genes
# Table: KnownGene
# Output: Selected fields from primary and related tables (then select name and geneSymbol)
ann_genes <- as.matrix(read.table("uc_to_gene_names_sex.txt", sep="\t", header=T)) # Use the correct filename!

#     hg19.knownGene.name hg19.kgXref.geneSymbol
#[1,] "uc001acw.2"        "MIR200B"             
#[2,] "uc001aiz.2"        "BC018779"            
#[3,] "uc001alm.1"        "AJAP1"          

###########################
###########################
# CONTINUE
###########################
###########################

# Add gene symbols to ann
rownames(ann) <- ann[,"feature.name"]
rownames(ann_genes) <- ann_genes[,1] # Watch out! Hard-coded!
ann <- cbind(ann, ann_genes[rownames(ann), 2]) # Watch out! Hard-coded!

# Add gene symbols to sigf
# This is done in a bit complicated way since all sites don't have a nearest gene (ann is missing some rows).
gene.symbol <- vector(length=nrow(sigf))
gene.symbol[as.numeric(ann[,"target.row"])] <- ann[,5] # Watch out! Hard-coded!
sigf <- cbind(sigf, gene.symbol)

# Similarily, add the distance to feature
dist.to.feature <- vector(length=nrow(sigf))
dist.to.feature[as.numeric(ann[,"target.row"])] <- ann[,"dist.to.feature"]
sigf <- cbind(sigf, dist.to.feature)

# Similarily, add the feature strand
feature.strand <- vector(length=nrow(sigf))
feature.strand[as.numeric(ann[,"target.row"])] <- ann[,"feature.strand"]
sigf <- cbind(sigf, feature.strand)

# rownames(sigf) <- paste(sigf[,"chr"], sigf[,"start"], sep=":")

##################################################
# Add genomic parts (exon, intron, promoter)
##################################################

parts <- sigf.ann@members
rownames(parts) <- paste(sigf[,"chr"], sigf[,"start"], sep=":")

part <- rep("intergenic", nrow(parts))
parts <- cbind(parts, part)
parts[parts[,"exon"]=="1","part"] <- "exon"
parts[parts[,"intron"]=="1","part"] <- "intron"
parts[parts[,"prom"]=="1","part"] <- "promoter"

sigf <- cbind(sigf, parts[rownames(sigf),"part"])
colnames(sigf)[ncol(sigf)] <- "genomic part"

######################################################
# Write
######################################################

sigf <- sigf[,setdiff(colnames(sigf), "end"),]

write.table(sigf, output_path, sep="\t", quote=F)

# The output contains these columns:
# chr	start	strand	context	pval	adj.pval	FDR	total.1	meth.1	total.0	meth.0	dif	type	gene.symbol	dist.to.feature	feature.strand	genomic part





