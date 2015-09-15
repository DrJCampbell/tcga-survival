#
# R functions for the tcga-sruvival project
#

get_observed_drugs <- function(x){
	drugs <- NULL
	for(i in nrow(x)){
		drug_data <- read.table(
			file=x[i,"drugs"],
			sep="\t",
			header=TRUE,
			comment.char="",
			stringsAsFactors=FALSE
			)
		drugs <- unique(
			c(
				drugs,
				drug_data[,"pharmaceutical_therapy_drug_name"]
				)
			)
	}
	return(drugs)
}

get_overall_survival <- function(x){
	overall_survival <- NULL
	for(i in 1:nrow(x)){
		study_id <- rownames(x)[i]
		followup <- read.table(
			file=x[i,"followup"],
			sep="\t",
			header=TRUE,
			comment.char="",
			stringsAsFactors=FALSE
			)
			# drop the first two rows (after the header)
			# they are supplementary headers
			followup <- followup[-c(1,2),]
		for(j in 1:nrow(followup)){
			# last_contact_days_to
			# death_days_to
			patient <- followup[j,"bcr_patient_barcode"]
			last_contact <- followup[j,"last_contact_days_to"]
			days_to_death <- followup[j,"death_days_to"]
			patient_overall_survival <- NA
			outcome <- NA
			if(followup[j,"death_days_to"] != "[Not Applicable]"){
				patient_overall_survival <- followup[j,"death_days_to"]
			}else{
				if(followup[j,"last_contact_days_to"] != "[Not Available]"){
					patient_overall_survival <- followup[j,"last_contact_days_to"]
				}	
			}
			vital_status <- followup[j,"vital_status"]
			tumor_status <- followup[j,"tumor_status"]
			if(vital_status == "Alive"){
				outcome <- vital_status
			}else{
				# if someone is dead and tumor free
				# they must not have died because of
				# the tumor?
				if(vital_status == "Dead" & tumor_status != "TUMOR FREE"){
					outcome <- vital_status
				}
			}
			overall_survival <- rbind(
				overall_survival,
				c(
					study_id,
					patient,
					patient_overall_survival,
					outcome
					)
				)
		}
	}
	overall_survival <- data.frame(
		study_id=as.character(overall_survival[,1]),
		patient=as.character(overall_survival[,2]),
		overall_survival=as.numeric(overall_survival[,3]),
		outcome=as.character(overall_survival[,4])
		)

	# at this stage each patient can have multiple entries.
	# we need to get the list of patients and for each one
	# select the largest value (i.e. not yet dead) from the
	# set of multiple follow-up vists. This needs to be
	# done within each study
	overall_survival_filtered <- NULL
	studies <- unique(overall_survival$study_id)
	for(study in studies){
		patients <- unique(overall_survival[which(
			overall_survival$study_id == study
			),"patient"])
		for(patient in patients){
			patient_data <- overall_survival[which(overall_survival$study_id == study & overall_survival$patient == patient),]
			latest_date_for_patient_row <- which(patient_data$overall_survival == max(patient_data$overall_survival))
			latest_patient_data <- patient_data[latest_date_for_patient_row[1],]
			overall_survival_filtered <- rbind(
				overall_survival_filtered,
				latest_patient_data
				)
		}
	}
	return(overall_survival_filtered)
}


get_progression_free_survival <- function(x){
	progression_free_survival <- NULL
	tumor_free_patient_data <- NULL
	for(i in 1:nrow(x)){
		study_id <- rownames(x)[i]
		nte <- read.table(
			file=x[i,"nte"],
			sep="\t",
			header=TRUE,
			comment.char="",
			stringsAsFactors=FALSE
			)
		# also read in the followup file so we can find tumour free
		# patients. Need to add these to the data set at the end
		followup <- read.table(
			file=x[i,"followup"],
			sep="\t",
			header=TRUE,
			comment.char="",
			stringsAsFactors=FALSE
			)
		# drop the first two rows (after the header)
		# they are supplementary headers
		nte <- nte[-c(1,2),]
		followup <- followup[-c(1,2),]
		for(j in 1:nrow(nte)){
			patient <- nte[j,"bcr_patient_barcode"]
			new_tumor_event_dx_days_to <- nte[j,"new_tumor_event_dx_days_to"]
			new_neoplasm_event_type <- nte[j,"new_neoplasm_event_type"]
			progression_free_survival <- rbind(
				progression_free_survival,
				c(
					study_id,
					patient,
					new_tumor_event_dx_days_to,
					new_neoplasm_event_type
					)
				)
		}
		j <- NULL
		for(j in 1:nrow(followup)){
			patient <- followup[j,"bcr_patient_barcode"]
			last_contact <- followup[j,"last_contact_days_to"]
			tumor_status <- followup[j,"tumor_status"]
			if(tumor_status == "TUMOR FREE"){
				tumor_free_patient_data <- rbind(
					tumor_free_patient_data,
					c(
						study_id,
						patient,
						last_contact,
						tumor_status
						)
					)
			}
		}
	}
	progression_free_survival <- data.frame(
		study_id=as.character(progression_free_survival[,1]),
		patient=as.character(progression_free_survival[,2]),
		progression_free_survival=as.numeric(progression_free_survival[,3]),
		outcome=as.character(progression_free_survival[,4])
		)

	# filter the progression free survival data to find the earliest
	# visit where a new tumour event was recorded. The add to this
	# patients from the followup file that are tumour free - using the
	# date of the latest visit to record the non-event
	progression_free_survival_filtered <- NULL
	studies <- unique(progression_free_survival$study_id)
	for(study in studies){
		patients <- unique(progression_free_survival[which(
			progression_free_survival$study_id == study
			),"patient"])
		for(patient in patients){
			patient_data <- progression_free_survival[which(
				progression_free_survival$study_id == study &
				progression_free_survival$patient == patient
				),]
			earliest_date_for_patient_row <- which(
				patient_data$progression_free_survival == min(patient_data$progression_free_survival, na.rm=TRUE)
				)
			earliest_patient_data <- patient_data[earliest_date_for_patient_row[1],]
			progression_free_survival_filtered <- rbind(
				progression_free_survival_filtered,
				earliest_patient_data
				)
		}
	}
	# now get the data for tumour free patients
	tumor_free_patients <- setdiff(
		tumor_free_patient_data[,2],
		progression_free_survival[,2]
		)
		
	tumor_free_patient_data <- data.frame(
		study_id=as.character(tumor_free_patient_data[,1]),
		patient=as.character(tumor_free_patient_data[,2]),
		progression_free_survival=as.numeric(tumor_free_patient_data[,3]),
		outcome=as.character(tumor_free_patient_data[,4])
		)
	
	# get the latest visit date for each tumour-free
	# patient
	tumor_free_patient_data_filtered <- NULL
	for(patient in unique(tumor_free_patient_data[,"patient"])){
		patient_data <- tumor_free_patient_data[which(
			tumor_free_patient_data[,"patient"] == patient
			),]
		latest_patient_visit <- which(
			patient_data[,3] == max(patient_data[,3], na.rm=TRUE)
			)[1]
		tumor_free_patient_data_filtered <- rbind(
			tumor_free_patient_data_filtered,
			patient_data[latest_patient_visit,]
			)
	}
	
	progression_free_survival_filtered <- rbind(
		progression_free_survival_filtered,
		tumor_free_patient_data_filtered[(tumor_free_patient_data_filtered[,2] %in% tumor_free_patients),]
		)
	return(progression_free_survival_filtered)
}









