#' @title calcHp
#' @description A function to calculate under a specific hypotheses: a given evidence and given refData
#' @export

calcHp = function(evidData,refData,nC,popFreq,kit,pC,lambda,threshT,nDone=4) { 
  uselocs <- names(popFreq)#[!names(popFreq)%in%"AMEL"] #EXCLUDE AMEL
  getAlleles = function(x) {
 #   if(length(x)==1) x <- rep(x,2)
   if(length(x)==1) return(as.character()) #dont consider if only 1 allele
    return(x)
  }
  refD <- list() #reference list must have other structure than considered here...  
  evidData <- lapply(evidData,function(x) x[uselocs]) #loci to consider
  for(loc in uselocs) refD[[loc]] <- lapply(refData,function(x) getAlleles(x[[loc]]$adata)) #extract reference
  if(is.null(refData)) refD = NULL
  data <- Qassignate(evidData, popFreq[uselocs],refD,incS=FALSE,incR=FALSE) #popFreq must be given with correct order?

  #Calc Lik(E|Hp)
  condOrder <- NULL
  if(!is.null(refData)) condOrder <- 1:length(refData) #conditional order

  if(nC>2) fitHp <- contLikMLEpara(nC,evidData,data$popFreq,data$refData,condOrder=condOrder,xi=0,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
  if(nC<=2) fitHp <- contLikMLE(nC,evidData,data$popFreq,data$refData,condOrder=condOrder,xi=0,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)

  return(fitHp) #return only matchlist
}# end  