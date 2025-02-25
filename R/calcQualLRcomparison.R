#' @title calcQualLRcomparison
#' @description  Function calculating LR using EuroForMix for remaining matches (from Step 1)
#' @param DBmix A list with evidence information [[sample]][[loc]] = list(adata,hdata)
#' @param DBref A (nR x nL) matrix with reference information. LIST not efficient 
#' @param matchlist A matrix with overview of what references are matches each of samples (SampleName,Reference)
#' @param popFreq A list of allele frequencies for a given population.
#' @param pC A numeric for allele drop-in probability. Can be a vector.
#' @param maxC Maximum number of contributors possible to assume. Default is 6.
#' @param useMinK1 A boolean of whether minimum number of contributors to test should be one.
#' @param normalize A boolean of whether normalization should be applied or not.
#' @param minFreq The freq value included for new alleles. Default NULL is using min.observed in popFreq.
#' @param maxIter Maximum number of iterations for nlm
#' @param useMAC Whether to use maximum allele counter (MAC) as the estimated NOC
#' @export

#library(casesolver);gui()
calcQualLRcomparison = function(DBmix,DBref,matchlist,popFreq,pC=0.05,maxC=6,useMinK1=TRUE,normalize=FALSE,minFreq=NULL, maxIter=5, useMAC=FALSE) { 
  #DBmix<<-DBmix;DBref<<-DBref;matchlist<<-matchlist;popFreq<<-popFreq;
  #pC=0.05;maxC=5;useMinK1=TRUE;normalize=TRUE;minFreq=NULL; maxIter=5
  Qallele = "99"
  
  log10LR <- numContr <- rep(NA,nrow(matchlist)) #vector to store LR values and assumed number of contributors
  locs <-  intersect( names(popFreq), names(DBmix[[1]])) #loci to consider (overlap)
  print(paste0("Calculating QUAL based LR for ",nrow(matchlist)," comparisons..."))
   
  #Update drop-in parameter to be per-marker specific
  if(length(pC)==1) pC = setNames(rep(pC,length(locs)),locs)
  
  #Obtain list of reference data: notice that single-alleles are removed from calculation
  refLIST <- casesolver::tabToListRef(tab=DBref,setEmpty=TRUE) #FORCING 1-alleles to be empty
  
systime <- system.time( {
  unEvid <- unique(matchlist[,1]) #get unique evidence 
  for(ss in unEvid) { #for each unique stain we estimate number of contr.
    # ss = unEvid[1]
    #Notice: Empty markers important because of information about allele dropouts.
    sample <- lapply(DBmix[ss],function(x) x[locs]) #extract sample
    Qset <- euroformix::Qassignate(sample, popFreq[locs],incS=FALSE,normalize=normalize,minF=minFreq)
    maxCsample = maxC #store upper limit of NOC to traverse (may depend on sample)
    
    #traverse number of contr under Hd.
    if(useMinK1) { #TRUE if K=1 contributors should be lower number of contributors
      nClow <- 1 
    } else {
      nClow <- ceiling(max(sapply(sample[[1]],function(x) length(x$adata)))/2) #get lower boundary of #contr
      if(useMAC) maxCsample = nClow #dont iterate higher NOC than this if method to use
    }
    bestfoo <- NULL
    for(nC in nClow:maxCsample) {
      foohd <- euroformix::calcQualMLE(nC,Qset$samples,Qset$popFreq, prC=pC,fst=0, maxIter = maxIter)
      if(is.null(bestfoo)) {
        bestfoo <- foohd
      } else { #if not first
        if( (foohd$min+1) < bestfoo$min ) { #if new max was more than 1 better (AIC criterion)
          bestfoo <- foohd
        } else {
          break; #stop if not better
        }
      }
      bestfoo$nC = nC 
    }
    loghd <- -bestfoo$min #maximum
    nC <- bestfoo$nC  #number of contr to use
    #1/(1+exp(-bestfoo$est)) #estimated dropout
  
    #Calc. Hp to get LR:
    whatR <- matchlist[which(ss == matchlist[,1]),2] #what refs to consider for particular sample
    for(rr in whatR ) { #for each references
      if(!rr%in%rownames(DBref))  next #skip if not found
      refD <- list() #reference list must have other structure than considered here...
      for(loc in locs) refD[[loc]] <- lapply(refLIST[rr],function(x) x[[loc]]$adata) #extract reference
      
      #NEED TO TAKE INTO ACCOUNT RARE ALLELES (Good idea to use Qassignate function)
      QsetHp <- euroformix::Qassignate(sample, popFreq[locs],refD,incS=FALSE,incR=FALSE,normalize=normalize,minF=minFreq) #popFreq must be given with correct order?
      foohp <- euroformix::calcQualMLE(nC,QsetHp$samples,QsetHp$popFreq, QsetHp$refData, condOrder = 1, prC=pC,fst=0, maxIter = maxIter)
      lr = (-foohp$min - loghd)/log(10)  #obtain LR on  log10 scale
      insind <- matchlist[,1]==ss &  matchlist[,2]==rr #index to insert
      log10LR[insind]  <- lr  #insert LR 
      numContr[insind] <- nC #insert number of contributors
    } #end for each references given unique stain
    ii <- which(ss==unEvid)  
    print(paste0(round(ii/length(unEvid)* 100), "% LR qual calculation complete...")) 
  } #end for each stains
})[3]
  print(paste0("FINISHED: Calculating qual LR took ",ceiling(systime), " seconds"))
  
  #Add score to match list:
  matchlist <- cbind(matchlist,log10LR,numContr)
  
  return(list(MatchList=matchlist)) #return only matchlist
}# end  calcLRcomparison 
