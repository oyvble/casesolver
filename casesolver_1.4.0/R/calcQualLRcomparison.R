#' @title calcQualLRcomparison
#' @description  Function calculating LR using EuroForMix for remaining matches (from Step 1)
#' @param DBmix A list with evidence information [[sample]][[loc]] = list(adata,hdata)
#' @param DBref A (nR x nL) matrix with reference information. LIST not efficient 
#' @param matchlist A matrix with overview of what references are matches each of samples (SampleName,Reference)
#' @param popFreq A list of allele frequencies for a given population.
#' @param pC A numeric for allele drop-in probability. Default is 0.
#' @param maxC Maximum number of contributors possible to assume. Default is 6.
#' @param useMinK1 A boolean of whether minimum number of contributors to test should be one.
#' @export

calcQualLRcomparison = function(DBmix,DBref,matchlist,popFreq,pC=0.05,maxC=6,useMinK1=TRUE) { 
  require(forensim) 
  
  log10LR <- numContr <- rep(NA,nrow(matchlist)) #vector to store LR values and assumed number of contributors
  locs <- names(popFreq)[toupper(names(popFreq))%in%toupper(names(DBmix[[1]]))] #find loci to consider (only those in evid)
  print(paste0("Calculating QUAL based LR for ",nrow(matchlist)," comparisons..."))
  pDv = c(0.1,0.35,0.7) #Necessary to look for potential better start point than (0.1) to make optimizer more robust
  iterlim0=5 #max number of iterations in optimization
  
  #Model calculation:
  neglikhd <- function(pD) {
       pDhd <-  rep(1/(1+exp(-pD)),nC)
       hdvec <- rep(1,length(locs))
       for(loc in locs) hdvec[which(loc==locs)] <- likEvid( Evidlist[[loc]],T=NULL,V=NULL,x=nC,theta=0, prDHet=pDhd, prDHom=pDhd^2, prC=pC, freq=popFreqQ[[loc]])
       return( -sum(log(hdvec)) )
  }
  neglikhp <- function(pD) {
     pDhp <- rep(1/(1+exp(-pD)),nC)
     hpvec <- rep(1,length(locs))
     for(loc in locs) hpvec[which(loc==locs)] <- likEvid( Evidlist[[loc]],T=Reflist[[loc]],V=NULL,x=nCList[[loc]],theta=0, prDHet=pDhp, prDHom=pDhp^2, prC=pC, freq=popFreqQ[[loc]])
     return( -sum(log(hpvec)) )
  }

 systime <- system.time( {
  unEvid <- unique(matchlist[,1]) #get unique evidence 
  for(ss in unEvid) { #for each unique stain we estimate number of contr.
   # ss = unEvid[39]
   #Notice: Empty markers important because of information about allele dropouts.
   sample <- lapply(DBmix[ss],function(x) x[locs]) #extract sample
   Qset <- euroformix::Qassignate(sample, popFreq[locs],incS=FALSE)

   #DATA IS MADE READY FOR OPTIMIZATION
   Evidlist <- popFreqQ <- list()
   for(loc in locs) {
   #loc=locs[1]
    adata <- sample[[1]][[loc]]$adata
    popFreqQ[[loc]]  <- Qset$popFreq[[loc]]
    if(length(adata)==0) {
      adata=0 #is empty
    } else {
      adata <- 1:length(adata) #use length
      names(popFreqQ[[loc]]) <- 1:length(popFreqQ[[loc]])
    }
    Evidlist[[loc]] <-adata #evidence to store
   }
 
   #traverse number of contr.
   if(useMinK1) { #TRUE if K=1 contributors should be lower number of contributors
     nClow <- 1 
   } else {
     nClow <- ceiling(max(sapply(sample[[1]],function(x) length(x$adata)))/2) #get lower boundary of #contr
   }
   bestfoo <- NULL
   for(nC in nClow:maxC) {

     #find suitable optim startpoint
  	tmp = Vectorize(neglikhd)( log(pDv/(1-pDv)) )
     pD0 = pDv[which.min(tmp)]
     foohd <- nlm(Vectorize(neglikhd),log(pD0/(1-pD0)), iterlim=iterlim0)

     if(is.null(bestfoo)) {
      bestfoo <- foohd
     } else { #if not first
      if( (foohd$min+1) < bestfoo$min ) { #if new max was more than 1 better (AIC criterion)
       bestfoo <- foohd
      }else {
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
    #rr=whatR[1]
    Reflist <- nCList <- list() #create the reflist with encoded alleles
    for(loc in locs) {
    #loc=locs[4]
     freq <- Qset$popFreq[[loc]]
     refA <- unlist(strsplit(DBref[rownames(DBref)==rr,colnames(DBref)==loc],"/"))  #get alleles
     if(length(refA)%in%c(0,1)) { #MUST HANDLE THAT REF ON LOCUS CAN BE MISSING OR 1 allele: ADD 1 UNKNOWN + SET TO ZERO
      refA = NULL #set as empty
      nCList[[loc]] <- nC  #number of unknowns = #contributors if missing
     } else {
      refA <- match( refA, names(freq) ) #get index in frequency
      refA[is.na(refA)] <- length(freq) #last index 
      nCList[[loc]] <- nC-1 #number of unknowns = #contributors - 1
     }
     Reflist[[loc]] <- refA
    } #end 
 
    #find suitable optim startpoint
    tmp = Vectorize(neglikhp)( log(pDv/(1-pDv)) )
    pD0 = pDv[which.min(tmp)] 
    foohp <- nlm(Vectorize(neglikhp),log(pD0/(1-pD0)), iterlim=iterlim0 )
    #1/(1+exp(-foohp$est)) #estimated dropout
  
    insind <- matchlist[,1]==ss &  matchlist[,2]==rr #index to insert
    log10LR[insind]  <-  (-foohp$min - loghd)/log(10)  #insert LR on log10 scale
    numContr[insind] <- nC
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
