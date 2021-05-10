#!/bin/bash -l
#SBATCH -J PQLseq_cord_blood_RRBS
#SBATCH -o Results/output_PQLseq_%A_%a.txt
#SBATCH -e Results/errors_PQLseq_%A_%a.txt
#SBATCH -t 36:00:00
#SBATCH -n 1
#SBATCH -p small
#SBATCH --mem=3G
#SBATCH --array=1-92
#

module load r-env-singularity/3.6.3

# The input arguments should be
# 1: path to meth matrix (nbr of methylated reads of each sample, each site)
# 2: path to total matrix (nbr of total reads of each site, each sample)
# 3: path to treatment vector
# 4: path to kinship matrix
# 5: output name
# 6: path to confounding covariates matrix
# 7: Slurm task ID

# Location of the coverage filtered count matrices:
FOLDER="/path_to_project/Step2_coverage_filtering/Results/"
# Location of the other inputs (treatment, coverage and kinship)
FOLDER_METADATA="/path_to_project/Step3_PQLseq/Data/"
OUTPUT_FOLDER="/path_to_project/Step3_PQLseq/Results/"

METH="${FOLDER}meth_noSNP_filtered.txt"
TOTAL="${FOLDER}total_noSNP_filtered.txt"
TREATMENT="${FOLDER_METADATA}treatment_173_sex.txt"
KINSHIP="${FOLDER_METADATA}kinship_SNP_profile_correlations.txt"
OUTPUT_NAME="${OUTPUT_FOLDER}results_$SLURM_ARRAY_TASK_ID"
COVARIATES="${FOLDER_METADATA}covariates_173_samples.txt"

srun singularity_wrapper exec Rscript --no-save Step2_run_PQLseq_with_pseudo_count.R $METH $TOTAL $TREATMENT $KINSHIP $OUTPUT_NAME $COVARIATES $SLURM_ARRAY_TASK_ID


