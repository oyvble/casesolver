#' @title calcQuanMLE
#' @description A function to calculate under a specific hypotheses: a given evidence and given refData
#' @param evidData A list with evidence data evidData[[sample]][[locus]]
#' @param refData A list with reference data refData[[refname]][[locname]]$adata
#' @param condOrder Conditional vector (hypothesis)
#' @param NOC Assumed number of contributors
#' @param nnTK an environment object from stored CaseSolver object
#' @param isWOE Indicating whether fst should be used or not
#' @param verbose Whether progress should be printed
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @export

calcQuanMLE = function(evidData,refData,condOrder,NOC,nnTK,isWOE=FALSE,verbose=FALSE,knownRef=NULL) { 

  mod <- casesolver::getModelSettings(nnTK) #Obtain model settings
  if(!isWOE) mod$fst = 0*mod$fst #set as zero if not WOE

  #Prepare which loci to use
  locs <-  intersect( names(mod$popFreq), names(evidData[[1]])) #loci to consider
  refD <- list() #reference list must have other structure than considered here...  
  evidData <- lapply(evidData,function(x) x[locs]) #loci to consider
  if(is.null(refData)) {
    refD = NULL
  } else {
    for(loc in locs) refD[[loc]] <- lapply(refData,function(x) x[[loc]]$adata) #extract reference
  }
  
  #prosess data (regarding frequencies and rare alleles etc)
  data <- euroformix::Qassignate(evidData, mod$popFreq[locs],refD,incS=FALSE,incR=FALSE,normalize=mod$normalize,min=mod$minFreq) #popFreq must be given with correct order?
  
  MLEfit <- euroformix::contLikMLE(NOC,data$samples,data$popFreq,data$refData,condOrder=condOrder,xi=mod$xiBW,xiFW=mod$xiFW,prC=mod$pC,lambda=mod$lambda,nDone=mod$nDone,threshT=mod$threshT,kit=mod$kit, fst=mod$fst, verbose=verbose, knownRef=knownRef)
  return(MLEfit) #fit
}# end  