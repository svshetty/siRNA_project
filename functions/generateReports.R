### functions to create reports in cellHTS 
generateReports <- function(reportType,profile,raw_data,normalized_data=NULL,scored_data=NULL){
    switch(reportType,
       raw={generateReports_raw(raw_data,profile);},
       normalized={generateReports_normalized(raw_data,normalized_data,profile);},
       score={generateReports_score(raw_data,normalized_data,scored_data,profile);}
	)
}

generateReports_raw <- function(raw_data,profile){
	for (cellline in profile$celllines){
		x=raw_data$cellHTS[[cellline]]
		reportflname=paste(profile$cellHTS_reportfolder,profile$pathSep,cellline,profile$pathSep,"report-raw",sep="")
	  setSettings(list(plateList=list(reproducibility=list(include=TRUE),
                                   intensities=list(include=TRUE), average=list(include=TRUE,map=TRUE))))
	  writeReport(raw=x,outdir=reportflname, force=TRUE) 
    
  
		#to produce scatterplots of intensity values across all plates
		QCfolderFull = paste(reportflname,"QCplots",sep=profile$pathSep)
		if (!file.exists(QCfolderFull))
		dir.create(QCfolderFull)
		numReplicates_thisCellline = length(unlist(profile$celllines_inRawdata[cellline]))

		b=Data(x)
		t=x@featureData
  
    pdf(file=paste(QCfolderFull, profile$pathSep, "Boxplot of Intensity by Probe.pdf",sep=""))
    boxplot(b ~ wellAnno(x),col="pink", outline=FALSE)
    dev.off()
		
    controlst=get("controlStatus",t@data) 
		plotColor=vector(mode="character",length=length(controlst)) 
		nmap=names(profile$ColorMap)

		for (i in 1:length(profile$ColorMap)){plotColor[controlst==nmap[i]]=profile$ColorMap[i]}

		for (i in 1:numReplicates_thisCellline){
			pdf(file=paste(QCfolderFull, profile$pathSep, "rawIntensityDotPlot","Replicate", i, ".pdf",sep=""))
			plot(b[,i,1],col=unlist(plotColor),
            			 main=paste("Raw intensities for replicate ",i,sep=""), 
            			 xlab="Well Index", ylab="Raw Intensity")   
			dev.off()
      
		}

		if (numReplicates_thisCellline>1){
			for (i in unlist(unique(t@data$plate))){
				thisPlateIndex = t@data$plate==i
				for (j in 1:(numReplicates_thisCellline-1)){
					for (k in (j+1):numReplicates_thisCellline){

						pdf(file=paste(QCfolderFull, profile$pathSep,"ScatterplotofReplicates","_plate_",i,"Replicate",j, "and Replicate",k, ".pdf", sep=""))
						plot(b[thisPlateIndex,j,1],b[thisPlateIndex,k,1],col=unlist(plotColor),
						main="Correlation between replicates:Raw Intensities", 
						sub= paste("plate",i,sep=""), 
						xlab=paste("Replicate",j,sep=""), ylab= paste("Replicate",k,sep=""))
						dev.off()

					}
				}
			}
		}

	}
}

generateReports_normalized <- function(raw_data,normalized_data,profile){
  if (is.null(normalized_data)){
    stop("no normalized data")
  }
  else {
   for (cellline in profile$celllines){
    for (normMethods in profile$normalization_mtds){
      x=raw_data$cellHTS[[cellline]]
      xn=normalized_data[[cellline]][[normMethods]] 
      reportflname=paste(profile$cellHTS_reportfolder,cellline,"report-normalized",normMethods,sep=profile$pathSep)
      writeReport(raw=x,normalized=xn,outdir=reportflname, force=TRUE)
    }
   }
  }
}

generateReports_score <- function(raw_data,normalized_data,scored_data,profile){
  if (is.null(normalized_data)){
    stop("no normalized data")}
  else {
   for (cellline in profile$celllines){
    for (normMethods in profile$normalization_mtds){
      x=raw_data$cellHTS[[cellline]]
      xn=normalized_data[[cellline]][[normMethods]]
      xsc=scored_data[[cellline]][[normMethods]] 
  
      setSettings(list(plateList=list(reproducibility=list(include=TRUE,map=TRUE),
                                      intensities=list(include=TRUE,map=TRUE)),
                                      screenSummary=list(scores=list(range=c(-4,8),map=TRUE)),
                                      controls=list(col=profile$cellHTScolorMap)))
      
      reportflname=paste(profile$cellHTS_reportfolder,cellline,"report-scored",normMethods,sep=profile$pathSep)
      writeReport(raw=x,normalized=xn,scored=xsc,outdir=reportflname, force=TRUE)            
      b=Data(xsc)
      pdf(file=paste(reportflname, profile$pathSep, "Boxplot of Z-score by Probe.pdf",sep=""))
      boxplot(b ~ wellAnno(x),col="lightblue", outline=FALSE)
      dev.off()
    }
   }
  }
}