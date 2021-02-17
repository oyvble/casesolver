#' @title calcQuanLRcomparison
#' @description  Function calculating LR using EuroForMix for remaining matches (from Step 1)

#' @param DBmix A list with evidence information [[sample]][[loc]] = list(adata,hdata)
#' @param DBref A (nR x nL) matrix with reference information. LIST not efficient 
#' @param matchlist A matrix with overview of what references are matches each of samples (SampleName,Reference,MAC)
#' @param popFreq A list of allele frequencies for a given population.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: Returned by euroformix::getKit()
#' @param xiBW Model parameter for Backward stutter. NULL indicates stutter model.
#' @param xiFW Model parameter for Forward stutter. NULL indicates stutter model.
#' @param pC A numeric for allele drop-in probability. 
#' @param lambda Parameter in modeled peak height shifted exponential model. 
#' @param threshT The detection threshold given. 
#' @param nDone Maximum number of random evaluations nlm-optimizing routing. Default is 4.
#' @param maxC Maximum number of contributors possible to assume. Default is 3.
#' @param nContr A vector of number of contributors for each of the row in the matchlist. If NULL (default), the max(nA)/2 rule is applied.
#' @param normalize A boolean of whether normalization should be applied or not.
#' @param minFreq The freq value included for new alleles. Default NULL is using min.observed in popFreq.
#' @export

calcQuanLRcomparison = function(DBmix,DBref,matchlist,popFreq,kit,xiBW=0,xiFW=0,pC=0.05,lambda=0.01,threshT=50,nDone=4,maxC=3,nContr=NULL,normalize=FALSE,minFreq=NULL) { 
  #Output 1=matchlist: Unique Match matrix (normalized number of allele match counting) for all combinations
  #Output 2=storedFitHp: Stored fit under Hp
  #Aim: Calcualate LR for all situations in matchlist
 
  calcMLE =  function(cond=NULL) {
  	 return( euroformix::contLikMLE(nC,sample,data$popFreq,data$refData,condOrder=cond,xi=xiBW,xiFW=xiFW,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE) )
  }
  
  log10LR <- numContr <- rep(NA,nrow(matchlist)) #vector to store LR values and assumed number of contributors
  storedFitHp<- vector(mode="list",nrow(matchlist)) #list-object for each row
  locs <-  intersect( names(popFreq), names(DBmix[[1]])) #loci to consider (overlap)
  
  #Obtain list of reference data: notice that single-alleles are removed from calculation
  refLIST <- casesolver::tabToListRef(tab=DBref,setEmpty=TRUE) #FORCING 1-alleles to be empty
  
  print(paste0("Calculating QUAN based LR for ",nrow(matchlist)," comparisons... this may take a while (hours)"))
  systime <- system.time( {
  unEvid <- unique(matchlist[,1]) #get unique evidence 
  for(ss in unEvid) { #for each unique stain we estimate number of contr.
    # ss = unEvid[1]
    #Notice: Empty markers important because of information about allele dropouts.
  
    sample <- lapply(DBmix[ss],function(x) x[locs]) #extract sample
    data <- euroformix::Qassignate(sample, popFreq[locs],incS=FALSE,incR=FALSE,normalize=normalize,minF=minFreq) #don't include stutters, use all loci
    if(!is.null(nContr)) {
      if(length(nContr)!=nrow(matchlist)) stop("The length of the argument nContr must be qual the number of rows in matchlist")
      nClow <- unique(as.integer(nContr[matchlist[,1]==ss])) #estimated number of contributors from qualLR
    } else {
      nClow <- ceiling(max(sapply(sample[[1]],function(x) length(x$adata)))/2) #get lower boundary of #contr
    }
    nC <- min(nClow,maxC) #assumed number of contributors: Restrict to 3 unknowns by default can be changed in GUI
    print(paste0("Assumed number of contributors for sample ",ss,": ",nC))
    
    #Calc Lik(E|Hd):   
    fitHd <- calcMLE()  #euroformix::contLikMLE(nC,sample,data$popFreq,xi=xi,xiFW=xiFW,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
    
   #NOTICE: These values can be stored for later for deconvolution
    loghd <- fitHd$fit$loglik #used under all combination for this stain.
   
    #Calc. Hp to get LR:
    whatR <- matchlist[which(ss == matchlist[,1]),2] #what refs to consider for particular sample
    for(rr in whatR ) { #for each references
     #rr=whatR[1]
     if(!rr%in%rownames(DBref)) { #if reference not found in dataList, we will create list:
      next #skip if not found
     }
     refD <- list() #reference list must have other structure than considered here...
     for(loc in locs) refD[[loc]] <- lapply(refLIST[rr],function(x) x[[loc]]$adata) #extract reference
     data <- euroformix::Qassignate(sample, popFreq[locs],refD,incS=FALSE,incR=FALSE,normalize=normalize,minF=minFreq) #popFreq must be given with correct order?
  
     #Calc Lik(E|Hp)   
     fitHp <- calcMLE(cond=1) #euroformix::contLikMLE(nC,sample,data$popFreq,data$refData,condOrder=1,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
     insind <- matchlist[,1]==ss &  matchlist[,2]==rr #index to insert
     log10LR[insind]  <-  (fitHp$fit$loglik  - loghd)/log(10)  #insert LR on log10 scale
     numContr[insind] <- nC
     storedFitHp[insind] <- list(fitHp) #store object
    } #end for each references given unique stain
    ii <- which(ss==unEvid)  
    print(paste0(round(ii/length(unEvid)* 100), "% LR quan calculation complete...")) 
 } #end for each stains
 })[3]
 print(paste0("FINISHED: Calculating Quan LR took ",ceiling(systime), " seconds"))

 #Add score to match list:
 matchlist <- cbind(matchlist,log10LR,numContr)

 return(list(MatchList=matchlist,storedFitHp=storedFitHp)) #return only matchlist
}# end  calcLRcomparison 
