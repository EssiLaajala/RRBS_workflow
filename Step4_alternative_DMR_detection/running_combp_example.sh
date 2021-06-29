#!/bin/bash -l
#SBATCH -J combp_permuted_epidural
#SBATCH -t 00:55:00
#SBATCH -n 1
#SBATCH -p small
#SBATCH --mem=1G
#

# Cloned this repository: https://github.com/brentp/combined-pvalues. There's an excellent Readme. This is just a documentation of how I used it:

# The P-values from PQLseq (or RADMeth beta-binomial regression or any other differential methylation analysis) were used as input in the form of a bed-file
# Please notice that for combp the bed-file needs to be sorted alphabetically (chr10 should come before chr2 etc.)

# In this example, the raw p-values were in column 5, window size was 200, step size 1. I also tried window size 500, step size 10 but the results were very similar.

# Estimate the autocorrelation
python ../combined-pvalues/cpv/acf.py -d 1:200:1 -c 5 ../Data/sorted_permuted_epidural.bed > ../Results/acf_permuted_epidural.bed

# Perform Stouffer-Lipták-Kechris correction:
python ../combined-pvalues/cpv/slk.py --acf ../Results/acf_permuted_epidural.bed -c 5 ../Data/sorted_permuted_epidural.bed > ../Results/pvals_slk_permuted_epidural.bed

# Detect peaks (candidate regions with lots of small spatially adjusted p-values):
python ../combined-pvalues/cpv/peaks.py --dist 200 --seed 0.1 ../Results/pvals_slk_permuted_epidural.bed > ../Results/pvals_regions_permuted_epidural.bed

# Get slk-adjusted Sidák corrected p-values for differential methylation at each candidate region:
python ../combined-pvalues/cpv/region_p.py -p ../Data/sorted_permuted_epidural.bed -r ../Results/pvals_regions_permuted_epidural.bed -s 1 -c 5 > ../Results/regions_sig_permuted_epidural.bed

