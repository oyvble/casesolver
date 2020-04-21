#' @title calcHp
#' @description A function to calculate under a specific hypotheses: a given evidence and given refData
#' @param evidData A list with evidence data
#' @param refData A list with reference data
#' @param nC Assumed number of contributors
#' @param popFreq A list with population frequences popFreq[[locusname]
#' @param kit Selected kit
#' @param pC Dropin probability
#' @param lambda Dropin PH hyperparameter 
#' @param threshT Detection threshold used
#' @param xi Stutter proportion parameter
#' @param nDone Number of required optimisations
#' @export

calcHp = function(evidData,refData,nC,popFreq,kit,pC,lambda,threshT,xi,nDone=2) { 
  uselocs <- names(popFreq)[toupper(names(popFreq))%in%toupper(names(evidData[[1]]))] #use loci given in data
  efmversion = packageVersion("euroformix") #obtain efm version
  if(is.null(xi) && efmversion < "3.0.0") uselocs <- setdiff(uselocs,"AMEL") #EXCLUDE AMEL IF STUTTER CONSIDERD (QUICK SOLUTION) 

  getAlleles = function(x) {
  if(length(x)==1) return(as.character()) #dont consider if only 1 allele
    return(x)
  }
  refD <- list() #reference list must have other structure than considered here...  
  evidData <- lapply(evidData,function(x) x[uselocs]) #loci to consider
  for(loc in uselocs) refD[[loc]] <- lapply(refData,function(x) getAlleles(x[[loc]]$adata)) #extract reference
  if(is.null(refData)) refD = NULL
  data <- euroformix::Qassignate(evidData, popFreq[uselocs],refD,incS=FALSE,incR=FALSE) #popFreq must be given with correct order?

  #Calc Lik(E|Hp)
  condOrder <- NULL
  if(!is.null(refData)) condOrder <- 1:length(refData) #conditional order

#samples=evidData;popFreq=data$popFreq;refData=data$refData
  if(nC>2) fitHp <- euroformix::contLikMLEpara(nC,evidData,data$popFreq,data$refData,condOrder=condOrder,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=TRUE)
  if(nC<=2) fitHp <- euroformix::contLikMLE(nC,evidData,data$popFreq,data$refData,condOrder=condOrder,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)

  return(fitHp) #return only matchlist
}# end  