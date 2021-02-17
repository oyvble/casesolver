#' @title getModelSettings
#' @description Obtain model settings from environment
#' @param nnTK an environment object from stored CaseSolver object
#' @export
#load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/TutorialDataCaseSolver/cs_error.Rdata")
getModelSettings = function(nnTK) {
  popFreq = get("popFreq",envir=nnTK) #get population freqs
  setupModel = get("setupModel",envir=nnTK) #get model options
  setupAdvanced = get("setupAdvanced",envir=nnTK)
  setupRare = get("setupRare",envir=nnTK) #get options of rare alleles
  
  kit = NULL #no kit specified by default
  xiBW <- xiFW <- 0 #no stutters assumed by default
  if( setupModel$degrad==1 ) kit = getEnvirKit(nnTK) #if degrad model
  if( setupModel$stuttBW==1 ) xiBW = NULL #if stutter model
  if( setupModel$stuttFW==1 ) xiFW = NULL #if stutter model
  
  #obtain global settings:
  pC=setupModel$dropinC
  lambda=setupModel$dropinL
  threshT=setupModel$threshT
  fst=setupModel$fst
  
  #Checking for marker based settings (update params:
  setupMarkers = get("setupMarkers",envir=nnTK)
  if(!is.null(setupMarkers)) {# && length(setupMarkers)==5) {
    if( !all(setupMarkers[[1]]%in%names(popFreq)) ) {
      stop("Wrong markers considered in marker specific settings. Please re-configure!")
      #return(NULL) #return from function
    }
    threshT = setNames(setupMarkers[[2]],setupMarkers[[1]])
    pC = setNames(setupMarkers[[3]],setupMarkers[[1]])
    lambda = setNames(setupMarkers[[4]],setupMarkers[[1]])
    fst = setNames(setupMarkers[[5]],setupMarkers[[1]])
  }
  nDone=setupAdvanced$nDone
  
  #Rare allele settings
  normalize = as.logical(setupRare$normalize)
  minFreq = NULL #default is NULL (using minfreq as min observed)
  if(!is.na(setupRare$minFreq) && nchar(setupRare$minFreq)>0) minFreq = as.numeric(setupRare$minFreq) #ensure numbers
    
  return(list(popFreq=popFreq,kit=kit,xiBW=xiBW,xiFW=xiFW,
         pC=pC,lambda=lambda,threshT=threshT, fst=fst,nDone=nDone,
         normalize=normalize,minFreq=minFreq))
}