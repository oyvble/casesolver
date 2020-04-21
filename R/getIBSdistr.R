#' @title getIBSdistr
#' @description  A function to get the distribution of IBS for random individuals
#' @param popFreq A list with population frequences popFreq[[locusname]
#' @param verbose Boolean of whether progress should be printed to console
#' @export

getIBSdistr = function(popFreq,verbose=TRUE) {
  #Theoretical distribution:
  popFreq = sapply(popFreq,function(x) x/sum(x)) #all alleles must sum to one (requirement of getGlist function)
  Glist <- euroformix::getGlist(popFreq)
  locs <- names(Glist)
  propagM <- numeric()
  for(loc in locs) {
   amat <- Glist[[loc]]$G
   pmat <- Glist[[loc]]$Gprob
   nR <- length(pmat) #span of uniques
   #get marginal of number of matches: 
   mac<- matrix(0,nR,nR) 
   for(i in 1:nrow(amat)) {
    t1 <- amat[i,1]==amat[,1] | amat[i,1]==amat[,2]
    t2 <- amat[i,2]==amat[,1] | amat[i,2]==amat[,2]
    mac[i,] <- as.numeric(t1)+as.numeric(t2)
    if( amat[i,1]==amat[i,2]) {
     ind <- amat[,1]!=amat[,2] & mac[i,]==2 #those not same but got MAC=2
     mac[i,ind] <- mac[i,ind] - 1 #subtrac 1
    }
   }
   pmac <- outer(pmat,pmat)
   #calc marginals of number of matches:
   agg <- aggregate(c(pmac),by=list(c(mac)),sum)

   if(which(locs==loc)>1) {
    newP <- outer(pmac,propagM[,2])
#    newM <- round(log(outer(exp(mac),exp(propagM[,1]))))
    newM <- outer(mac, propagM[,1],FUN="+") #adding instead of summing
    agg <- aggregate(c(newP),by=list(c(newM)),sum)   #calculate marginals for updated vector again:
   }
   propagM <- cbind(agg[,1],agg[,2])
   ii <- which(loc==locs)  
   if(verbose) print(paste0(round(ii/length(locs)* 100), "% completed")) 
  }
  rm(newP,newM) 
  gc() #clean memory
  return(propagM)
}  #End function