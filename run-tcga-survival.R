#
# example code to run tcga-survival-functions
#

setwd("~/Dropbox/tcga-survival")
source("./tcga-survival-functions.R")

# the inputs file needs to have the 
# following format:
# col1 study_id
# col2 drugs_file_name
# col3 followup_file_name
# col4 nte_file_name
data_to_process <- read.table(
	file="inputs.txt",
	header=TRUE,
	row.names=1,
	sep="\t",
	stringsAsFactors=FALSE
	)

# find all drug names used and report
observed_drugs <- get_observed_drugs(
	data_to_process
	)


write.table(
	observed_drugs,
	"observed_drugs.txt",
	col.names=FALSE,
	row.names=FALSE,
	quote=FALSE,
	sep="\t"
	)

# we need to somehow flag samples treated
# with one or more treatments of interest
# idealy, retaining some info about other
# classes.
#
# Try creating a matrix of all samples *
# all treatments. Then dynamically group
# cell lines with certain column names

samples_and_drugs <- get_samples_and_drugs(
	data_to_process,
	drugs=observed_drugs
	)

write.table(
	samples_and_drugs,
	file="samples_and_drugs.txt",
	col.names=TRUE,
	row.names=TRUE,
	sep="\t",
	quote=FALSE
	)


# at this point we need to hand-check
# the returned list of drugs observed.
# There are typos and variations on drug
# names to deal with.

drugs <- read.table(
	file="./drugs.txt",
	header=FALSE,
	sep="\t",
	stringsAsFactors=FALSE
	)

platinum_treated <- rep(0, times=nrow(samples_and_drugs))
names(platinum_treated) <- rownames(samples_and_drugs)
for(sample in names(platinum_treated)){
	platinum_treated[sample] <- max(
		samples_and_drugs[sample,which(as.matrix(drugs) %in% colnames(samples_and_drugs))],
		na.rm=TRUE
		)
}


# get survival estimates for each sample
# and study
overall_survival <- get_overall_survival(
	data_to_process
	)
write.table(
	overall_survival,
	"overall_survival.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t"
	)

progression_free_survival <- get_progression_free_survival(
	data_to_process
	)
write.table(
	progression_free_survival,
	"progression_free_survival.txt",
	col.names=TRUE,
	row.names=FALSE,
	sep="\t"
	)


