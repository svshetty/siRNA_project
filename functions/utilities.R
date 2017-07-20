### function for checking parameters in the normalization methods in cellHTS
parametercheck<- function(profile,normMethods){
  splitmethod=strsplit(normMethods,"_")
  methodPrefix=splitmethod[[1]][[1]]
  method=splitmethod[[1]][[2]]                                                   
  if (methodPrefix=="cellHTS"){
    # checking for compatibility of scale and log transformation of data
    if (profile$normalization_scale=="additive" && profile$normalization_log){
      stop("If scale is additive the log cannot be TRUE")};
    # assigning specific methods with the parameters as defined in cellHTS documentation
    if(method=="robustZscore"){method="median" 
                               profile$varianceAdjust="byPlate"}; 
    if(method=="Bscore" || method=="locfit"){profile$varianceAdjust="byPlate"}
  }
  
  return(list(profile=profile,methodPrefix=methodPrefix,method=method))                                                             
}
                                                              
                                                              
