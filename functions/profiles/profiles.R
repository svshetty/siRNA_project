#### General parameters

load_exp_profile <- function(ExpName,workdir,pathSep,optsShow,cellHTS_FilesReady=TRUE){

	#default values
	profile <- list(
		ExperimentName = ExpName,
		workdir=workdir,
		pathSep=pathSep,
		cellHTS_FilesReady=cellHTS_FilesReady,
		datafile="",
		Rows = c(LETTERS)[1:16],
		Cols = c(1:24),
		dataDelimiter = ',',
		naString = " ",
		totalWells = 384,
		cellHTS_verbose = FALSE,
		cellHTS_screenlog_empty = TRUE,
		ColorMap=list(pos="red", neg="blue", sample="black", other="green",empty="yellow"),
    cellHTScolorMap=c("sample"="black","neg"="blue", "controls"="#4DAF4A", "other"= "green","empty"= "#FF7F00","flagged"= "#A65628","act"= "#E41A1C","inh"= "yellow" ,"pos"= "red"),
		normalization_mtds=c("cellHTS_robustZscore"),
		#normalization_mtds=c("cellHTS_POC", "cellHTS_negatives", "cellHTS_NPI","cellHTS_robustZscore","cellHTS_median", "cellHTS_mean", "cellHTS_shorth", "cellHTS_Bscore", "cellHTS_locfit"),
		normalization_scale="additive", #(options available: "additive", "multiplicative")
		normalization_VarianceAdjust="none",#(options available: "byBatch", "byExperiment", "none")
		normalization_log=FALSE, #(options available:TRUE,FALSE)
		score_mtds="zscore", #options available:"none","zscore", "NPI"
		score_sign="-", #options available:"+", "-"
		summarize_mtds="max" #options available: "mean", "median", "max", "min", "closestToZero", "furthestFromZero", "rms"
	)


	# load exp specific parameters
	exp_set = FALSE
	switch(
		ExpName,
		EwingsKinaseScreens5Lines=
		{profile=load_profile_EwingsKinaseScreens5Lines(profile);
			exp_set=TRUE;},
		Myeloma_Dose_0_only_DGv3_library=
		{profile=load_profile_Myeloma_Dose_0_only_DGv3_library(profile);
			exp_set=TRUE;}

	)
	if (!exp_set) stop("unregonanized experiment name")


	#options
	options(warn=1) #0 -- (the default) warnings are stored until the top–level function returns
	#1 -- warnings are printed as they occur


	#show
	if (optsShow == "show"){
		allnames = names(profile)
		for (i in 1:length(profile)){
			cat(allnames[i]," -- ",as.character(profile[allnames[i]]),"\n");
		}
	}


	return(profile)
}


