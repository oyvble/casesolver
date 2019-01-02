#' @title calcQuanLRcomparison
#' @description  Function calculating LR using EuroForMix for remaining matches (from Step 1)

#' @param DBmix A list with evidence information [[sample]][[loc]] = list(adata,hdata)
#' @param DBref A (nR x nL) matrix with reference information. LIST not efficient 
#' @param matchlist A matrix with overview of what references are matches each of samples (SampleName,Reference)
#' @param popFreq A list of allele frequencies for a given population.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}. 
#' @param pC A numeric for allele drop-in probability. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param nDone Maximum number of random evaluations nlm-optimizing routing. Default is 4.
#' @param maxC Maximum number of contributors possible to assume. Default is 3.
#' @export

calcQuanLRcomparison = function(DBmix,DBref,matchlist,popFreq,kit,xi,pC,lambda,threshT,nDone=4,maxC=3) { 
 require(euroformix) 
#Output 1=MatchMatrix: Unique Match matrix (normalized number of allele match counting) for all combinations
#Output 2=MatchList: A match list giving all combinations having score (equal or) greater than threshMAC
  #Aim: Calcualate LR for all situations in matchlist
  #Question: USE AMEL OR NOT?

  getAlleles = function(x) {
#    if(length(x)==1) x <- rep(x,2)
    if(length(x)==1) return(as.character()) #dont consider if only 1 allele
    return(x)
  }
  
  log10LR <- numContr <- rep(NA,nrow(matchlist)) #vector to store LR values and assumed number of contributors
  storedFitHp<- vector(mode="list",nrow(matchlist)) #list-object for each row
  uselocs <-  names(popFreq)[toupper(names(popFreq))%in%toupper(names(DBmix[[1]]))] #loci to consider
  if(is.null(xi)) uselocs <- setdiff(uselocs,"AMEL") #EXCLUDE AMEL IF STUTTER CONSIDERD (QUICK SOLUTION) 

  print(paste0("Calculating QUAN based LR for ",nrow(matchlist)," comparisons... this may take a while (hours)"))
 systime <- system.time( {
  unEvid <- unique(matchlist[,1]) #get unique evidence 
  for(ss in unEvid) { #for each unique stain we estimate number of contr.
  # ss = unEvid[1]
  #Notice: Empty markers important because of information about allele dropouts.

  sample <- lapply(DBmix[ss],function(x) x[uselocs]) #extract sample
  data <- euroformix::Qassignate(sample, popFreq[uselocs],incS=FALSE,incR=FALSE) #don't include stutters, use all loci
  nClow <- ceiling(max(sapply(sample[[1]],function(x) length(x$adata)))/2) #get lower boundary of #contr

  if(!is.null(matchlist) && ncol(matchlist)>=5) { #Only considered if QualLR performed (used to estimate #contr)
   nClow <- as.integer(matchlist[matchlist[,1]==ss,5]) #estimated number of contributors from qualLR
  }
  nC <- min(nClow,maxC) #assumed number of contributors: Restrict to 3 unknowns by default can be changed in GUI
  print(paste0("Assumed number of contributors for sample ",ss,": ",nC))
  
  #Calc Lik(E|Hd): 
  if(nC>2) fitHd <- euroformix::contLikMLEpara(nC,sample,data$popFreq,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
  if(nC<=2) fitHd <- euroformix::contLikMLE(nC,sample,data$popFreq,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
 #NOTICE: These values can be stored for later for deconvolution
  loghd <- fitHd$fit$loglik #used under all combination for this stain.
 
  #Calc. Hp to get LR:
  whatR <- matchlist[which(ss == matchlist[,1]),2] #what refs to consider for particular sample
  for(rr in whatR ) { #for each references
   #rr=whatR[1]
   if(!rr%in%rownames(DBref)) { #if reference not found in dataList, we will create list:
    next #skip if not found
   } else {
    refD <- strsplit(DBref[rownames(DBref)==rr,colnames(DBref)%in%uselocs],"/") #extract reference
   } 
   refD2 <- list() #reference list must have other structure than considered here...
   for(loc in uselocs) {
     refD2[[loc]] <- list(getAlleles(refD[[loc]]))
     names(refD2[[loc]]) <- rr #insert name here. Important for showing TopEPG
   }
   data2 <- euroformix::Qassignate(sample, popFreq[uselocs],refD2,incS=FALSE,incR=FALSE) #popFreq must be given with correct order?

   #Calc Lik(E|Hp)
   if(nC>2) fitHp <- contLikMLEpara(nC,sample,data2$popFreq,data2$refData,condOrder=1,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
   if(nC<=2) fitHp <- contLikMLE(nC,sample,data2$popFreq,data2$refData,condOrder=1,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=FALSE)
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
