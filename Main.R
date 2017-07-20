##################
# HTSuite ver1.0 #
# 7/12/11
##################

### Init
rm(list = ls(all = TRUE))
pathSep = ifelse(.Platform$OS.type == "windows","\\","/")
workdir = getwd()
#load local functions
functionPath = paste(workdir,"functions",sep=pathSep)
listFuns = dir(path = functionPath, pattern = "*.R", all.files = F, full.names = T, recursive = T, ignore.case = T)
  for (rsource in listFuns) source(rsource)

### Experiment
ExpName = 'Myeloma_Dose_0_only_DGv3_library'
	# EwingsKinaseScreens5Lines
	# Myeloma_Dose_0_only_DGv3_library
profile = load_exp_profile(ExpName,workdir,pathSep,'noshow',FALSE)
# if cellHTS files are ready use the command below
#profile = load_exp_profile(ExpName,workdir,pathSep,'noshow')


### Data
require(cellHTS2)
require(stringr)
raw_data = read_raw_data(profile)
  #creates the following datasets
  #data$data_raw
  #data$augmented
  #data$cellHTS
   
# QC: first pass
# creates raw data QC plots
generateReports('raw',profile,raw_data)

### normalization
#function prototype depends on function calls to cellHTS
#normalizes the data using cellHTS methods
normalized_data=normalize_data(raw_data,profile)

### QC2: second pass
#creates normalized data plots
generateReports('normalized',profile,raw_data,normalized_data)


### hit identification
scored_data=score_data(normalized_data,profile)
generateReports('score',profile,raw_data,normalized_data,scored_data)

## (optional) RSV: multifple siRNAs for a gene