### function to create raw, augmeted and cellHTS data sets
read_raw_data <- function(profile){
  data_raw = read.csv(profile$datafile, sep=profile$dataDelimiter,header=T,check.names=T,na.strings = profile$naString,stringsAsFactors = FALSE)
  data_augmented = raw2augment(data_raw,profile)
  data_augmented_check(data_augmented,profile)
  
  if (!profile$cellHTS_FilesReady){
    augment2cellHTSwrite(data_augmented,profile)}
  data_cellHTS = cellHTS_dataread(profile)
  
# 	data_cellHTS = NULL
      
  return(list(raw=data_raw,augmented=data_augmented,cellHTS=data_cellHTS))
}

### function to convert raw data to augmented
raw2augment <- function(data_raw,profile){
	switch(profile$ExperimentName, 
		EwingsKinaseScreens5Lines={data_aug=raw2augment_Ewings(data_raw,profile);},
		Myeloma_Dose_0_only_DGv3_library={data_aug=raw2augment_Myeloma(data_raw,profile);}
	)

	return(data_aug)
}

create_cellHTS_Plateconfig<-function(data_augmented, profile){
	switch(profile$ExperimentName, 
		EwingsKinaseScreens5Lines={create_cellHTS_Plateconfig_Ewings(data_augmented,profile);},
		Myeloma_Dose_0_only_DGv3_library={create_cellHTS_Plateconfig_Myeloma(data_augmented,profile);}
	)
}

data_augmented_check <- function(data_augmented,profile){
  for (fld in profile$augDataFields) 
    if (is.null(unlist(data_augmented[fld])))
      stop("data_augmented_check failed")
}

### function to create cellHTS data from augmented data
augment2cellHTSwrite <- function(data_augmented,profile){
  plate_nums = unlist(unique(data_augmented$plate))
  
  for (cellline in profile$celllines){
    dir_name<-paste(profile$cellHTS_file_folder,cellline, sep=profile$pathSep)
    if (!file.exists(dir_name))
      dir.create(dir_name)
  
    dir_reportflname=paste(profile$cellHTS_reportfolder,cellline,sep=profile$pathSep)
    if (!file.exists(dir_reportflname))
      dir.create(dir_reportflname,recursive = TRUE)
    
    GeneID_flname=paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_Gene_ID",".txt", sep="")
    file.create(GeneID_flname,showWarnings=TRUE)
    write.table(cbind(data_augmented$plate,data_augmented$wells,data_augmented$gene, data_augmented$siRNA),
    file=GeneID_flname, row.names = FALSE,col.names=c("Plate", "Well","HFAID" ,"GeneID"), append = FALSE, ,quote=FALSE,sep = "\t")
    
    Description_flpath= paste(profile$cellHTS_file_folder,profile$pathSep,cellline, sep="")
    templateDescriptionFile(filename=paste(cellline,"_Description.txt",sep=""), path=Description_flpath, force=TRUE)
    #this file has to be filled in manually;   
  }
 
  ####creating data files within the function "augment2cellHTSwrite"
  
  for (cellline in profile$celllines){   
    numReplicates_thisCellline = length(unlist(profile$celllines_inRawdata[cellline]))
    Platelist_flname=paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,cellline,"_Platelist",".txt", sep="")
    file.create(Platelist_flname,showWarnings=TRUE)
    cat("Filename", "Plate", "Replicate", file=Platelist_flname, sep="\t","\n") 
    
    for (plate in plate_nums){
      thisPlateIndex = (data_augmented$plate == plate)
      
      for (r in 1:numReplicates_thisCellline){
        data_subset_lumvalues= data_augmented[[cellline]][thisPlateIndex,r]
        data_subset_wellID = data_augmented$well[thisPlateIndex]
        data_subset=data.frame(data_subset_wellID,data_subset_lumvalues)
        
        fl_name=paste(cellline, "P",as.character(plate), "R", as.character(r),".txt", sep="")
        fl_Fullname=paste(profile$cellHTS_file_folder,profile$pathSep,cellline,profile$pathSep,fl_name, sep="")
        
        write.table(cbind(fl_name,data_subset),file= fl_Fullname,append=FALSE, quote=FALSE,sep="\t", row.names = FALSE,col.names = FALSE)
        cat(fl_name, plate, r, file=Platelist_flname, sep="\t","\n",append=TRUE) 
        
        if (length(data_subset_lumvalues) != profile$totalWells){
          cat("number of wells = ",length(data_subset), "in plate ",plate,"cellline ",cellline,"replicate ",r,"\n")
          stop("ERROR")}
               
      }
    }
  }
  ##creating plate confirguration files"augment2cellHTSwrite"
  create_cellHTS_Plateconfig(data_augmented,profile)

}
      
### function for reading the raw datafiles created in cellHTS format and converting to cellHTS objects
cellHTS_dataread<-function(profile){
  data_cellHTS =  sapply(profile$cellline,function(x) NULL)
   
  for (cellline in profile$celllines){
    dataPath=paste(profile$cellHTS_file_folder,cellline, sep=profile$pathSep)
    experimentName=cellline
    x=readPlateList(paste(cellline,"Platelist.txt",sep="_"),name=experimentName,path=dataPath,verbose=profile$cellHTS_verbose)
	if (profile$cellHTS_screenlog_empty){
		x=configure(x,path=dataPath,
			descripFile=paste(cellline,"Description.txt",sep="_"), 
			confFile=paste(cellline,"PlateConf.txt",sep="_"))  
	}
	else {
		x=configure(x,path=dataPath,
			descripFile=paste(cellline,"Description.txt",sep="_"), 
			confFile=paste(cellline,"PlateConf.txt",sep="_"), 
			logFile=paste(cellline,"Screenlog.txt",sep="_"))
	}                                  
    data_cellHTS[[cellline]] = x
    
  }           
  return(data_cellHTS)
}

### function for normalizing the data from cellHTS objects created
normalize_data<- function(raw_data,profile){
	normalizedData = sapply(profile$cellline,function(x) NULL)             

	for (cellline in profile$celllines){
		x=raw_data$cellHTS[[cellline]]
		xn= sapply(profile$normalization_mtds,function(x) NULL)
		for (normMethods in profile$normalization_mtds){
			#do parameter check here
			newparameter=parametercheck(profile,normMethods)
			if (newparameter$methodPrefix=="cellHTS"){
				xn[[normMethods]]=normalizePlates(x,
					verbose=profile$cellHTS_verbose,
					scale=newparameter$profile$normalization_scale,
					log=newparameter$profile$normalization_log, 
					method=newparameter$method, 
					varianceAdjust= newparameter$profile$normalization_VarianceAdjust) 
          Zfactor=getZfactor(xn[[normMethods]],robust=TRUE)
    print(Zfactor)
			}
		}  
		normalizedData[[cellline]]=xn
  
	}
 return(normalizedData)              
}

score_data<- function(normalized_data,profile){
	scoredData = sapply(profile$cellline,function(x) NULL)             

	for (cellline in profile$celllines){
		xsc= sapply(profile$normalization_mtds,function(x) NULL)   
    dataPath=paste(profile$cellHTS_file_folder,cellline, sep=profile$pathSep)
		for (normMethods in profile$normalization_mtds){  
			xn=normalized_data[[cellline]][[normMethods]]
			xsc[[normMethods]]=scoreReplicates(xn,sign=profile$score_sign,method=profile$score_mtds)
			xsc[[normMethods]]=summarizeReplicates(xsc[[normMethods]],summary=profile$summarize_mtds)
      xsc[[normMethods]]=annotate(xsc[[normMethods]],geneIDFile=paste(cellline,"Gene_ID.txt",sep="_"), path=dataPath)
		} 
		scoredData[[cellline]]=xsc 
	}      
	return(scoredData)
}

