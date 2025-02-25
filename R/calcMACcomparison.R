#' @title calcMACcomparison
#' @description  A function calc number of alleles in ref-table which are within evidence-table
#' @param DBmix A (nS x nL) matrix with evidence (only alleles) information. LIST not efficient 
#' @param DBref A (nR x nL) matrix with reference information. LIST not efficient 
#' @param threshMAC Required MAC threshold for being a candidate match
#' @export

calcMACcomparison = function(DBmix,DBref,threshMAC) {
#Output 1=MatchMatrix: Unique Match matrix (normalized number of allele match counting) for all combinations
#Output 2=MatchList: A match list giving all combinations having score (equal or) greater than threshMAC
   locs <- colnames(DBmix) #get loci
   nL <- length(locs)

   #Ensure that these are matrix and not list
   DBmix = as.matrix(DBmix)
   DBref = as.matrix(DBref)
   
   #Find loci for refs having only 1 allele:
   isNA <- is.na(DBref) #get loci which is NA
   isHom <- matrix(sapply(strsplit( DBref, "/"),length),nrow=nrow(DBref ))==1  
   DBref[isHom & !isNA] <- paste0(DBref[isHom & !isNA],"/",DBref[isHom & !isNA]) #make sure it has two alleles per loci

   DBmixN <- rownames(DBmix) #get names of evidences
   DBrefN <- rownames(DBref) #get names of refs
   nS <- nrow(DBmix)
   if(nS==0) {
     print(paste0("No mixture evidence profiles detected. Nothing is done."))
     return(NULL) 
   }
   nR <- nrow(DBref)
   print(paste0("Calculating MAC for all ",nS*nR," comparisons: All refs against all stains"))
   bigMAC <- rep(0,nS*nR) #keep MAC in a vector (sample1-ref1,sample1-ref2,...,sample2-ref1 etc.)

  systime <- system.time( { #register timeusage of following comparison:
    for(ss in 1:nS) { #for each sample: Limited in how the samples are looking
    #ss=1
    bigInd <-  nR*(ss-1) + 1:nR  #index in bigMAC matrix
    macS <- nLocs <- rep(0,nR) #make vector for all references 
    for(loc in locs) { #for each locus: Vectorizing as much as possible!
     ll=which(loc==locs)
     if(is.na(DBmix[ss,ll]) || DBmix[ss,ll]=="")  next #skip if no data
     sttmp <- unique(unlist(strsplit(DBmix[ss,ll],"/"))) #get unique alleles only some alleles may have been given twice

     isna <- is.na(DBref[,ll]) | DBref[,ll]==""  #get which refs are non-zero. Not counted if missing!
     numWithin <- rep(0,sum(!isna)) #number of unique alleles of reference that are in stain
     if(length(numWithin)==0) next #skip if no references
     AvecList <- t(matrix(unlist(strsplit(DBref[!isna,ll],"/")),nrow=2)) #get alleles of refs (in case of 2)
     for(aa in sttmp) numWithin <- numWithin + rowSums(AvecList==aa) #sum up for each alleles in stain 
     macS[!isna] <- macS[!isna] + numWithin #add number of matching alleles
     nLocs[!isna] <- nLocs[!isna] + 1  #add locus  
#print(loc)
#print(macS)
#print(nLocs)
    } #end for each locus
    bigMAC[bigInd] <- macS/(2*nLocs)  #divide number of matching-alleles in refernce by maximum possible
    #hist(macS)
    #insert MAC to long vector
    if (ss%%100 == 0)  print(paste0(round(ss/nS* 100), "% MAC calculation complete")) 
   } #end number of samples
  })[3]
  print(paste0("Calculating MAC scores took ",ceiling(systime), " second(s)"))
  matchmat <- matrix(round(bigMAC,2),nrow=nR) #allele-match matrix
  colnames(matchmat) <- DBmixN 
  rownames(matchmat) <- DBrefN 
  #print(matchmat)

  #cbind(c(t(replicate(nR,rownames(DBmix)))),rep(rownames(DBref),nS),bigMAC)
  keepInd <- which(bigMAC>=threshMAC)
  #bigMAC[keepInd]
  Ctab <- cbind(  floor((keepInd-1)/nR) + 1, (keepInd-1)%%nR + 1 , bigMAC[keepInd]) #convert back indices
  colnames(Ctab) <- c("Evid","Refs","score") #note the order: stainID first
 
  #(1) Remove because it was the same stain:
  tab <- cbind( DBmixN[ Ctab[,1] ], DBrefN[ Ctab[,2] ], Ctab[,3] ) #get table of matches 
  colnames(tab) <- c("Evidence","Reference","MAC")

  #sort list-output 
  #tab=matchMACList$MatchList 
  #ord <- order( as.numeric(tab[,3]),tab[,2],decreasing=TRUE) #sort first by score, then by REF-name
  tab[,3] <- round(as.numeric(tab[,3]),2) #round score to 2 decimals
  return(list(MatchMatrix=matchmat,MatchList=tab))
} # end calcMACcomparison

