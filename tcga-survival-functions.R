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
#				print(overall_survival[,1])
#				print(overall_survival[,2])
#				print(overall_survival[,3])
#				print(overall_survival[,4])
		}
	}
	overall_survival <- data.frame(
		study_id=as.character(overall_survival[,1]),
		patient=as.character(overall_survival[,2]),
		overall_survival=as.numeric(overall_survival[,3]),
		outcome=as.character(overall_survival[,4])
		)
	return(overall_survival)
}


get_progression_free_survival <- function(x){
	progression_free_survival <- NULL
	for(i in 1:nrow(x)){
		study_id <- rownames(x)[i]
		nte <- read.table(
			file=x[i,"nte"],
			sep="\t",
			header=TRUE,
			comment.char="",
			stringsAsFactors=FALSE
			)
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
	}
	progression_free_survival <- data.frame(
		study_id=as.character(progression_free_survival[,1]),
		patient=as.character(progression_free_survival[,2]),
		progression_free_survival=as.numeric(progression_free_survival[,3]),
		outcome=as.character(progression_free_survival[,4])
		)
	return(progression_free_survival)
}









