
# Before this, the methylation extractor needs to be run once without trimming (it outputs M-bias files that can be used to decide how much to trim).
# I trimmed everything that differed from the methylation levels of middle positions by more than 3 sd (Read 1 3', Read 2 5' and Read 2 3')
# This is an example of reading M bias files and printing the names of those samples that need further trimming
# This is only intended as a documentation of what was done and is not a nice script for others to use. 
# If you want to use it, you need to look at the M bias files first to know how many lines to skip and how many to read (nrows)
# This doesn't work unless all M-bias files in the working directory are for the same read length

## Choose these:
# Working directory: It should contain the M-bias files.
# files <- list.files(pattern="txt") # files should now contain M-bias file names, one file for each sample (in this example the working directory contains no other txt files)
sd_limit <- # How many standard deviations away from the mean is considered extreme? I've used value 3.
middle <- # I used 10:90 (middle 80 % of read length 100)
nbr_rows_to_skip <- # For example, if you want to look at read 2 CpG M-bias in the example file, skip 318 rows
read_length <- # For read 2 it's 99 in the example file 
start_positions <- # If middle positions are 10:90, try 1:9 here first
end_positions <- # If middle positions are 10:90, try 91:read_length here first
column <- 4 # Methylation percentages are in column 4 (in the current Bismark methylation extractor M-bias file format)

###########
# Start
###########

to_be_trimmed <- vector()
for(i in files){
	mbias <- read.table(i, skip=nbr_rows_to_skip, nrows=read_length)
	# TRUE/FALSE-vectors (one value for each start position):
	above_normal_variation <- mbias[start_positions,column] > (mean(mbias[middle,column]) + sd_limit * sd(mbias[middle,column]))
	below_normal_variation <- mbias[start_positions,column] < (mean(mbias[middle,column]) - sd_limit * sd(mbias[middle,column]))
	# If any base at the start positions is above or below normal variation, the sample needs further trimming
	if(sum(above_normal_variation, na.rm=T) + sum(below_normal_variation, na.rm=T) > 0) to_be_trimmed <- c(to_be_trimmed,i)
}
to_be_trimmed
# This just printed M-bias file names that indicate more trimming is needed. How much? Try different start_positions. In this example, if you start from the 4th base, this prints an empty vector indicating that bases 4:9 are within the normal variation. It's enough to trim away the first base
# For an individual sample, the vectors above_normal_variation and below_normal_variation will also tell you where to trim. 

# This example sample had a value outside normal variation at position 3. Other samples had them at positions 1 or 2 (of read 2). For simplicity, first 3 bases were trimmed from the 5' end of read 2 of all samples. The other option is to do different trimmings for different samples.

###########
# End
###########

to_be_trimmed <- vector()
for(i in files){
	mbias <- read.table(i, skip=nbr_rows_to_skip, nrows=read_length)
	# TRUE/FALSE-vectors (one value for each end position):
	above_normal_variation <- mbias[end_positions,column] > (mean(mbias[middle,column]) + sd_limit * sd(mbias[middle,column]))
	below_normal_variation <- mbias[end_positions,column] < (mean(mbias[middle,column]) - sd_limit * sd(mbias[middle,column]))
	# If any base at the end positions is above or below normal variation, the sample needs further trimming
	if(sum(above_normal_variation, na.rm=T) + sum(below_normal_variation, na.rm=T) > 0) to_be_trimmed <- c(to_be_trimmed,i)
}
to_be_trimmed
# In this example it's enough to trim away one base. If you try end_positions = 91:98, it prints an empty vector

# All samples had the same pattern in this project: One outlier at the last base of read 2.


