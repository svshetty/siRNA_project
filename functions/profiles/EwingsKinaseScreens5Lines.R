####Experiment specific parameters

load_profile_EwingsKinaseScreens5Lines <- function(profile){
    
	profile$rootfolder = profile$workdir
	profile$datafile=paste(profile$rootfolder,'data',profile$ExperimentName,"Ewing's Kinase Screens 5 Lines Data.csv",sep=profile$pathSep)
	profile$cellHTS_file_folder = paste(profile$rootfolder,'data',profile$ExperimentName,sep=profile$pathSep)
	profile$cellHTS_reportfolder= paste(profile$rootfolder,'reports',profile$ExperimentName,sep=profile$pathSep)

	#in raw data 
	profile$celllines_inRawdata = list(TC32080328=c(3,4),
										TC71080324=c(5,6),
										SKES080403=c(7,8),
										RDES080417=c(9,10),
										GM080565=c(11,12))

	profile$negCtrl = c("siRNA Buffer","siRNA Buffer_","Non-Sil","Non-Sil_","Scrambled","Scrambled_","GFP","GFP_")
	profile$posCtrl = c("UBB_s1","UBBs1_")
	#profile$other= Reference (first 2 columns for all plates); use regular expression 
	
	#select "all" cell lines
# 	profile$celllines = names(profile$celllines_inRawdata)
	# OR selected cell lines
	profile$celllines = c("TC71080324")

	#in augmented data
	# required field for the list 
	profile$augDataFields = c("plate","wells","gene","siRNA",profile$celllines)
	profile$augDataNULLtemplate = sapply(profile$augDataFields,function(x) NULL)
	profile$plate_prefix="Plt"
	profile$plate_well_sep_inRawFile="-"
	profile$plate_well_inRaw_Col=1
	profile$gene_siRNA_inRaw_Col=2

	return(profile)
}

#### Create augmented dataset from raw

raw2augment_Ewings <- function(data_raw,profile){
  data_aug_Ewings = profile$augDataNULLtemplate

  #"plate"
  plates_unpurified = stringr::str_extract(data_raw[,profile$plate_well_inRaw_Col],paste(profile$plate_prefix,"[0-9]+",profile$plate_well_sep_inRawFile,sep=""))
  data_aug_Ewings$plate=as.numeric(stringr::str_extract(plates_unpurified,"[0-9]+"))
  
  #"wellRow and Col"
  wells_unpurified = stringr::str_extract(data_raw[,profile$plate_well_inRaw_Col],paste(profile$plate_well_sep_inRawFile,"[A-P][0-9]{1,2}$",sep=""))
  data_aug_Ewings$wells = stringr::str_replace(wells_unpurified,profile$plate_well_sep_inRawFile,"")
  
  #"gene"
  data_aug_Ewings$gene = stringr::str_extract(data_raw[,profile$gene_siRNA_inRaw_Col],"([A-Za-z0-9]+)")
  
  #"siRNA"
  data_aug_Ewings$siRNA = data_raw[,profile$gene_siRNA_inRaw_Col]
  #cell lines
  for (c in profile$celllines){
    data_aug_Ewings[[c]] = as.matrix(data_raw[,profile$celllines_inRawdata[[c]]])
  }  
  return(data_aug_Ewings)
}

### create Plate Confirguration file according to cellHTS format
  
create_cellHTS_Plateconfig_Ewings<-function(data_augmented,profile){     
  for (cellline in profile$celllines){
    PlateConfig_flname=paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_PlateConf",".txt", sep="")
    file.create(PlateConfig_flname,showWarnings=TRUE)
    cat("Wells:\t", profile$totalWells,"\n","Plates:\t", length(unique(data_augmented$plate)), "\n", append=TRUE, file=PlateConfig_flname,sep="" )
    cat("Plate", "Well","Content", file=PlateConfig_flname,"\n", append=TRUE,sep="\t")
    cat("*", "*", "sample",file=PlateConfig_flname,"\n", sep="\t", append=TRUE)
    cat("*", "[A-P]0[1-2]","other", file=PlateConfig_flname, "\n", sep="\t",append=TRUE)
    
    for (i in 1:length(data_augmented$siRNA)){
      if(any(data_augmented$siRNA[i]==profile$posCtrl)){cat(data_augmented$plate[i],data_augmented$wells[i], "pos",file=PlateConfig_flname,"\n", sep="\t", append=TRUE)}
      else if (any(data_augmented$siRNA[i]==profile$negCtrl)){cat(data_augmented$plate[i],data_augmented$wells[i], "neg",file=PlateConfig_flname,"\n", sep="\t", append=TRUE)}
    }
     
    ScreenLog_flname= paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_Screenlog",".txt", sep="")
    file.create(ScreenLog_flname,showWarnings=TRUE)
    cat("Plate\t", "Sample\t","Well\t", "Flag\t", "Comment", file=ScreenLog_flname,"\n", append=TRUE,sep="")
    #this file has to be filled in manually; rememeber to press 'Tab' after filling in data for each column and 'Enter' at the end of the line
  }   
}
  

  
  