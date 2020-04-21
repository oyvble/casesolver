#' @title calcRMP
#' @description Random match probability calculations
#' @details A function for calculating Random match probabilities/RMNE for each profiles (Evidence/Refs)
#' @param nnTK an environment object from stored CaseSolver object
#' @param sig Number of significant values to use
#' @return A list with elements evid and ref
#' @export

calcRMP = function(nnTK,sig=3) {  
  refL <- get("refDataTABLE",envir=nnTK) 
  mixL <- get("mixDataLIST",envir=nnTK) #import mixtures

  popL <- get("popFreq",envir=nnTK) #import popFreq
  if(is.null(popL)) return() #Return if no data
  sn <- names(mixL) #get evid names
  rn <- rownames(refL) #get ref names
  locs <- names(popL) #loci to consider
  locs  <- setdiff(locs,"AMEL") #remove AMEL
  print( paste0("Considered loci from freq file: ",paste(locs,collapse="/")) )
  minF = min(unlist(popL)) #use minimum observed # #0.001 #minimium freq inserted when missing alleles
  
  calcRMNE = function(X) { #import a list X[[loc]] with allele vector and compare incl alleles in popFreqs
    #notice newly observed alleles are assumed to have allele freq=0
    rmne = 1 #default rmne value
    for(loc in locs) {
      if( is.null(X[[loc]]) || length(X[[loc]])==0) next; #skip if empty or not existing
      av <- names(popL[[loc]]) #get allele names
      freqs <- popL[[loc]][match( X[[loc]],av )] #[av%in%X[[loc]]] #get allele freqs
      freqs[is.na(freqs)] <- minF #.001 #insert minimum allele equal 0.001 if prev. not seen, alternative=min(unlist(popL))
      rmne <- rmne*sum(freqs)^2 #apply formula
    }
    return(rmne) #return calculated
  }
  
  
  #CALCULATE RMNE FOR EVIDENCE
  evidList <- matrix("",ncol=3,nrow=length(sn))
  for(ss in sn) { #for each evidence
    ind <- which(sn==ss)
    rmne <- calcRMNE(X=lapply(mixL[[ss]][locs],function(x) x$adata)) #return only adata: DONT CONSIDER THE DETECTION THRESHOLD HERE
    evidList[ind,] <- c(paste0("#",ind),ss,format(rmne,digits=sig))
  }
  
  #CALCULATE RMP FOR REFERENCES
  popL = lapply(popL,function(x) x/sum(x)) #require sum to 1
  Glist = euroformix::getGlist(popL)
  rmp = rep(1,nrow(refL)) #default rmp value
  for(loc in locs) {
    sub = refL[,toupper(colnames(refL))==toupper(loc)] 
    G = paste0(Glist[[loc]]$G[,1],"/",Glist[[loc]]$G[,2])
    indUse = match(sub,G)
    isna = is.na(indUse) #get those not found
    rmp[!isna] = rmp[!isna]*Glist[[loc]]$Gprob[indUse[!isna]]
    rmp[isna] = rmp[isna]*minF #IN SITUATION OF MISSING 
  }
  refList = cbind(paste0("#",1:length(rn)),rn,format(rmp,digits=sig))
  return( list(evid=evidList,ref=refList) ) #return results
} #end function
