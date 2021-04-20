#!/usr/bin/env Rscript

## Script to filter and prioritize MT variant identified with MToolBox.
## Based on Maddie Couse Variant Prioritization : https://team.bcchr.ca/display/TGA/MToolbox+Mitochondrial+Analysis
## Developped by Solenne Correard on March 29, 2019
##Last update: SC, April 3rd, 2019

library(plyr)

#Open pedfile
ResultsDirectory=getwd()
setwd(ResultsDirectory)

#Read ped file
pedfiles=list.files(pattern=".ped")
if (length(pedfiles)>1) {
	print ("ERROR, several ped files")
} else {
	#Create the filtered csv for affected members of the family
	ped=read.csv(pedfiles, sep="\t", header=TRUE, na.strings=c("","NA"))
	Affected=subset(ped, ped[,6]=="2")
	Number_of_affected=nrow(Affected)
	
	for (i in (1: Number_of_affected)){
		Proband_ID_i<- gsub("_GRCh38", "", Affected[i,2])
		
		##Open the csv files
		MToolBox_files=list.files(pattern="MToolBox_")
		if (length(MToolBox_files)<1) {
			print ("ERROR, no MToolBox files")
		} else {
			proband_i_directory=paste0("MToolBox_", Proband_ID_i,"/OUT_", Proband_ID_i)
			setwd(proband_i_directory)
			annotation_CSV_proband_i=list.files(pattern=".annotation.csv")
			csv_proband_i=read.csv(annotation_CSV_proband_i, sep="\t", header=TRUE, na.strings=c("","NA"))
			setwd(ResultsDirectory)

			##Filter the proband file to keep variant with (HF>0.18 or empty score) & Nt.Variability>0.0026
			csv_proband_i_filtered= subset(csv_proband_i, ((csv_proband_i$HF>0.18 | is.na(csv_proband_i$HF)) & csv_proband_i$Nt.Variability <0.0026))
			nrow(csv_proband_i)
			nrow(csv_proband_i_filtered)

			#Rename the columns with the proband number
			colnames(csv_proband_i_filtered)=c("Sample_p", ".Variant.Allele", paste0(Proband_ID_i, "_HF_p"), paste0(Proband_ID_i, "_CI_lower.CI_upper_p"), "var_RSRS", "var_MHCS", "var_rCRS", "var_Haplogroup", "var_Other.Haplogroups", "var_Locus", "var_Nt.Variability", "var_Codon.Position", "var_Aa.Change", "var_Aa.Variability", "var_tRNA.Annotation", "var_Disease.Score", "var_RNA.predictions", "var_MutPred.pred", "var_MutPred.prob", "var_PolyPhen.2.HumDiv.pred", "var_PolyPhen.2.HumDiv.prob", "var_PolyPhen.2.HumVar.pred", "var_PolyPhen.2.HumVar.prob", "var_PANTHER.pred", "var_PANTHER.prob", "var_PhD.SNP.pred", "var_PhD.SNP.prob", "var_SNPs.GO.pred", "var_SNPs.GO.prob", "var_Mitomap.Associated.Disease.s.", "var_Mitomap.Homoplasmy", "var_Mitomap.Heteroplasmy", "var_Somatic.Mutations", "var_SM.Homoplasmy", "var_SM.Heteroplasmy", "var_ClinVar", "var_OMIM.link", "var_dbSNP.ID", "var_Mamit.tRNA.link", "var_PhastCons20Way", "var_PhyloP20Way", "var_AC.AN.1000.Genomes", "var_X1000.Genomes.Homoplasmy", "var_X1000.Genomes.Heteroplasmy")

			write.table(csv_proband_i_filtered, file=paste0("MToolBox_annotated_p", Proband_ID_i), sep="\t", row.names = FALSE, quote=FALSE)
			}
		}
	#Create the filtered csv for mother of the family
	Mother_ID <- gsub("_GRCh38", "", Affected[1,4])
	if (Mother_ID=="-9") {
		print ("No Mother in the ped file")
	} else {
		setwd(paste0("MToolBox_", Mother_ID,"/OUT_", Mother_ID, "/"))
		annotation_CSV_mother=list.files(pattern=".annotation.csv")
		csv_mother=read.csv(annotation_CSV_mother, sep="\t", header=TRUE, na.strings=c("","NA"))
		setwd(ResultsDirectory)
			
		#Rename the columns of the mother file with the mother number
		colnames(csv_mother)=c("Sample_m", ".Variant.Allele", paste0(Mother_ID, "_HF_m"), paste0(Mother_ID, "_CI_lower.CI_upper_m"), "var_RSRS", "var_MHCS", "var_rCRS", "var_Haplogroup", "var_Other.Haplogroups", "var_Locus", "var_Nt.Variability", "var_Codon.Position", "var_Aa.Change", "var_Aa.Variability", "var_tRNA.Annotation", "var_Disease.Score", "var_RNA.predictions", "var_MutPred.pred", "var_MutPred.prob", "var_PolyPhen.2.HumDiv.pred", "var_PolyPhen.2.HumDiv.prob", "var_PolyPhen.2.HumVar.pred", "var_PolyPhen.2.HumVar.prob", "var_PANTHER.pred", "var_PANTHER.prob", "var_PhD.SNP.pred", "var_PhD.SNP.prob", "var_SNPs.GO.pred", "var_SNPs.GO.prob", "var_Mitomap.Associated.Disease.s.", "var_Mitomap.Homoplasmy", "var_Mitomap.Heteroplasmy", "var_Somatic.Mutations", "var_SM.Homoplasmy", "var_SM.Heteroplasmy", "var_ClinVar", "var_OMIM.link", "var_dbSNP.ID", "var_Mamit.tRNA.link", "var_PhastCons20Way", "var_PhyloP20Way", "var_AC.AN.1000.Genomes", "var_X1000.Genomes.Homoplasmy", "var_X1000.Genomes.Heteroplasmy")
		
		write.table(csv_mother, file=paste0("MToolBox_annotated_m", Mother_ID), sep="\t", row.names = FALSE, quote=FALSE)
		}
		
	list_filtered_csv_p=list.files(pattern="MToolBox_annotated_p")
	list_filtered_csv_m=list.files(pattern="MToolBox_annotated_m")
	if (length(list_filtered_csv_p)==1 & length(list_filtered_csv_m)==0){
		#If there is only one proband file (Singleton)
		proband_csv = read.csv(list_filtered_csv_p, sep="\t", header=TRUE, na.strings=c("","NA"))
		proband_csv_ordered= proband_csv[order(proband_csv$var_Disease.Score), ]
		proband_csv_ordered_ordered= proband_csv_ordered[c(2, 1 , 3, 4, 10, 37, 11, 13, 14, 36, 38, 30, 31, 32, 33, 34, 35, 42, 43, 44, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 40, 41, 5, 6, 7, 8, 9, 12, 15, 39)]
		write.table(proband_csv_ordered_ordered, file="MToolBox_annotated.txt", sep="\t", row.names = FALSE, quote=FALSE)

	} else if (length(list_filtered_csv_p)>1 & length(list_filtered_csv_m)==0){
		#If there is no mother and several affected individuals	
		csv_filtered_merged = read.csv(list_filtered_csv_p[1], sep="\t", header=TRUE, na.strings=c("","NA"))
		for (i in 2:length(list_filtered_csv_p)){
	    	currentFile = read.csv(list_filtered_csv_p[i], sep="\t", header=TRUE, na.strings=c("","NA"))
	  		csv_filtered_merged = merge(csv_filtered_merged, currentFile, by=".Variant.Allele", all=TRUE)
	  		csv_filtered_merged =rename(csv_filtered_merged, c("Sample_p.x"="Sample_p"))
	  		}
	} else if (length(list_filtered_csv_p)==1 & length(list_filtered_csv_m)==1){
	  	#If there is a mother and 1 affected individual
	  	proband_csv = read.csv(list_filtered_csv_p, sep="\t", header=TRUE, na.strings=c("","NA"))
		mother_csv = read.csv(list_filtered_csv_m, sep="\t", header=TRUE, na.strings=c("","NA"))
	  	csv_filtered_merged = merge(proband_csv, mother_csv, by=".Variant.Allele", all.x=TRUE)
	} else if (length(list_filtered_csv_p)>1 & length(list_filtered_csv_m)==1){
	  	#If there is a mother and several affected individuals
	  	mother_csv = read.csv(list_filtered_csv_m, sep="\t", header=TRUE, na.strings=c("","NA"))
	  	csv_filtered_merged_p = read.csv(list_filtered_csv_p[1], sep="\t", header=TRUE, na.strings=c("","NA"))
		for (i in 2:length(list_filtered_csv_p)){
	    	currentFile = read.csv(list_filtered_csv_p[i], sep="\t", header=TRUE, na.strings=c("","NA"))
	  		csv_filtered_merged_p = merge(csv_filtered_merged_p, currentFile, by=".Variant.Allele", all=TRUE)
	  		csv_filtered_merged_p =rename(csv_filtered_merged_p, c("Sample_p.x"="Sample_p"))
	  		}
	  	csv_filtered_merged = merge(csv_filtered_merged_p, mother_csv, by=".Variant.Allele", all.x=TRUE)	  	
	} else {
		print ("No solution with the files present")	
	  	}
	
	##Check if there is a final file already, if not, create it.
	list_final_files=list.files(pattern="MToolBox_annotated.txt")  	
	if (length(list_final_files)==0) {
	  	#Order the variants in the file in :
			#Decreasing Disease.Score (High disease score top of the table - threshold for significance disease score > 0.4311)
		csv_filtered_merged_ordered= csv_filtered_merged[order(csv_filtered_merged$var_Disease.Score.x), ]

	##Merge Samples number in one column, separated by a ","
		csv_filtered_merged_ordered_sample= csv_filtered_merged_ordered[grep("Sample_",colnames(csv_filtered_merged_ordered))]
		Samples=t(t(apply(csv_filtered_merged_ordered_sample, 1, paste, collapse=",")))
		
		##Merge HF values in one column, separated by a ","
		csv_filtered_merged_ordered_HF= csv_filtered_merged_ordered[grep("_HF_",colnames(csv_filtered_merged_ordered))]
		HF=t(t(apply(csv_filtered_merged_ordered_HF, 1, paste, collapse=",")))

		##Merge CI values in one column, separated by a ","
		csv_filtered_merged_ordered_CI= csv_filtered_merged_ordered[grep("_CI_lower.CI_upper_",colnames(csv_filtered_merged_ordered))]
		CI_lower_CI_upper=t(t(apply(csv_filtered_merged_ordered_CI, 1, paste, collapse=",")))

		##Merge the new data frame with the previous one
		csv_filtered_merged_ordered_merged=cbind(csv_filtered_merged_ordered, Samples, HF, CI_lower_CI_upper)

		##Select the columns of interest
		csv_filtered_merged_ordered_head=head(csv_filtered_merged_ordered_merged[,c(".Variant.Allele", "Samples", "HF", "CI_lower_CI_upper",colnames(csv_filtered_merged_ordered_merged)[grep(".x",colnames(csv_filtered_merged_ordered_merged))])])
		
		##Reorder the columns to fit with the excel template
		csv_filtered_merged_ordered_head_ordered <- csv_filtered_merged_ordered_head[c(1, 2, 3, 4, 10, 37, 11, 13, 14, 36, 38, 30, 31, 32, 33, 34, 35, 42, 43, 44, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 40, 41, 5, 6, 7, 8, 9, 12, 15, 39)]

		write.table(csv_filtered_merged_ordered_head_ordered, file="MToolBox_annotated.txt", sep="\t", row.names = FALSE, quote=FALSE)

}}

list_file_to_remove=list.files(pattern="MToolBox_annotated_") 	
file.remove(list_file_to_remove)
		













