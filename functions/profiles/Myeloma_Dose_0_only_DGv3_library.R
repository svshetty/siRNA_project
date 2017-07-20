#### Experiment specific parameters
load_profile_Myeloma_Dose_0_only_DGv3_library <- function(profile){

	profile$rootfolder = profile$workdir
	profile$datafile=paste(profile$rootfolder,'data',profile$ExperimentName,"Myeloma_Dose_0_only_DGv3_library.csv",sep=profile$pathSep)
	profile$cellHTS_file_folder = paste(profile$rootfolder,'data',profile$ExperimentName,sep=profile$pathSep)
	profile$cellHTS_reportfolder= paste(profile$rootfolder,'reports',profile$ExperimentName,sep=profile$pathSep)

	#in raw data
 #input the column number where the luminescence value is stored in csv file
 
  profile$celllines_inRawdata = list(SP=c(7)) #use with original data
  
  
 
	profile$negCtrl = c("siRNA Buffer","AS-Non-Sil","Non-Sil","GFP")
	profile$posCtrl = c("UBBs1","ACDC")
	profile$other=c("Media 1", "Media 2")


	#select "all" cell lines
	profile$celllines = names(profile$celllines_inRawdata)

	#in augmented data
	# required field for the list 
	profile$augDataFields = c("plate","wells","gene","siRNA",profile$celllines)
	profile$augDataNULLtemplate = sapply(profile$augDataFields,function(x) NULL)
	profile$plate_prefix="SP"
	#profile$plate_well_sep_inRawFile="-"
	profile$plate_inRaw_Col=2
	profile$siRNA_inRaw_Col=4
	profile$well_inRaw_Col=3
	profile$gene_inRaw_Col=6
	
	return(profile)
}

#### Create augmented dataset from raw
raw2augment_Myeloma <- function(data_raw,profile){
	data_aug_Myeloma = profile$augDataNULLtemplate
	#"plate"
	plates_unpurified = stringr::str_extract(data_raw[,profile$plate_inRaw_Col],paste(profile$plate_prefix,"([0-9]+)",sep=""))
	data_aug_Myeloma$plate=as.numeric(stringr::str_extract(plates_unpurified,"([0-9]+)"))
	#str(data_aug_Myeloma$plate)
	#"wellRow and Col"

	data_aug_Myeloma$wells= data_raw[,profile$well_inRaw_Col]
	#str(data_aug_Myeloma$wells)
	#convert from A1 to A01
	data_aug_Myeloma$wells=stringr::str_replace(data_aug_Myeloma$wells,"([A-Pa-p])([1-9]$)","\\10\\2")

	#"gene"
	data_aug_Myeloma$gene = data_raw[,profile$gene_inRaw_Col]
	
	#str(data_aug_Myeloma$gene)
	#"siRNA"
	data_aug_Myeloma$siRNA = data_raw[,profile$siRNA_inRaw_Col]
	#str(data_aug_Myeloma$siRNA)
	
 #cell lines
	for (c in profile$celllines){
		data_aug_Myeloma[[c]] = as.matrix(data_raw[,profile$celllines_inRawdata[[c]]])
	} 
	#str(data_aug_Myeloma[[c]])

	return(data_aug_Myeloma)
}

### create Plate Confirguration file according to cellHTS format
create_cellHTS_Plateconfig_Myeloma<-function(data_augmented,profile){     
	for (cellline in profile$celllines){
		PlateConfig_flname=paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_PlateConf",".txt", sep="")
		file.create(PlateConfig_flname,showWarnings=TRUE)
		cat("Wells:\t", profile$totalWells,"\n","Plates:\t", length(unique(data_augmented$plate)), "\n", append=TRUE, file=PlateConfig_flname,sep="" )
		cat(paste("Plate", "Well","Content",sep="\t"), "\n",file=PlateConfig_flname, append=TRUE,sep="")
		cat(paste("*", "*", "sample",sep="\t"),"\n",file=PlateConfig_flname, append=TRUE,sep="")


		for (i in 1:length(data_augmented$siRNA)){
			if(any(data_augmented$siRNA[i]==profile$posCtrl)){cat(paste(data_augmented$plate[i],data_augmented$wells[i], "pos", sep="\t"),"\n", file=PlateConfig_flname,append=TRUE,sep="")}
			else if(any(data_augmented$siRNA[i]==profile$other)){cat(paste(data_augmented$plate[i],data_augmented$wells[i], "other", sep="\t"), "\n",file=PlateConfig_flname,append=TRUE,sep="")}
			else if (any(data_augmented$siRNA[i]==profile$negCtrl)){cat(paste(data_augmented$plate[i],data_augmented$wells[i], "neg", sep="\t"),"\n", file=PlateConfig_flname,append=TRUE,sep="")}
		}

		ScreenLog_flname= paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_Screenlog",".txt", sep="")
		file.create(ScreenLog_flname,showWarnings=TRUE)
		cat(paste("Plate", "Sample","Well", "Flag", "Comment",sep="\t"),"\n",file=ScreenLog_flname,"\n", append=TRUE,sep="")
		#this file has to be filled in manually; rememeber to press 'Tab' after filling in data for each column and 'Enter' at the end of the line
	}   
}
