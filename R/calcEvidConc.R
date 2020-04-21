#' @title calcEvidConc
#' @description Function for calculating evidence concordance
#' @param nnTK an environment object from stored CaseSolver object (includes data)
#' @param nLarge A number indicating VERY large numeber of references
#' @return A list with canidate results
#' @export

calcEvidConc = function(nnTK,nLarge=10000) {  
  DBref <- get("mixDataTABLE",envir=nnTK)  #get table
  allrefsn <- rownames( DBref )
  nR <- length(allrefsn)
  if(nR<2) return() #return if less than 2
  locs <- colnames(DBref)
  
  popL = get("popFreq",envir=nnTK) #get popFreq
  if(!is.null(popL)) locs = locs[ toupper(locs)%in%toupper(names(popL)) ] #Use ONLY THOSE IN POPFREQ (AVOIDS Y-markers) 
  print(paste0("Calculating ",nR*(nR-1)/2," comparisons for ",length(locs)," loci."))
  
  allrefs <- list()
  for(loc in locs) {
    tmp <- DBref[,colnames(DBref)==loc] 
    isna <- is.na(tmp) | tmp=="" #missing
    if(all(isna)) next; #skip if all is missing
    tmp = strsplit(tmp,"/") #string split
    nAlleles = sapply(tmp,length) #get number of alleles
    unAlleles = setdiff(unique(nAlleles),0) #remove 0 as a possibility
    amat <- matrix("",nrow=nR,ncol=max(unAlleles)) #require 
    for(nA in unAlleles) { #for each number of alleles, we insert alleles into matrix
      indUse <- which(nAlleles==nA)
      amat[indUse,1:nA] <- t(matrix(unlist(tmp[indUse]),nrow=nA)) #gets error if none have two alleles
    }
    #swap = amat[,2]<amat[,1] #sorting already done (assumed only?)
    #amat[swap,] = amat[swap,2:1]
    allrefs[[loc]] <- amat
  }
  rm(amat,tmp);gc()
  
  nlocs <- mac <- matrix(0,nrow = nR, ncol = nR) #Ordinary matrix (could use bigmemory to optimize memory usage)
  locs <- names(allrefs) #get non-empty loci
  for(loc in locs) { #calculate the pairwise comparison
    unmat <- allrefs[[loc]] #consider all outcome
    #unmat = unique(amat) #get only unique rows (can be done in earlier stage?)
    nAlleles = rowSums(unmat!="") #get number of alleles 
    unA = sort(setdiff(unique(unmat[nAlleles>0]),"")) #get unique alleles
    
    #Preprosessing: Indicate which alleles that are in different refs
    MACmatrix = matrix(FALSE,nrow=nrow(unmat),ncol=length(unA)) #create matrix of 
    for(aa in unA) MACmatrix[ rowSums(unmat==aa)==1 , which(aa==unA) ] = TRUE  #for each alleles
    
    for(i in 1:nrow(unmat)) { #for each reference
      if(nAlleles[i]==0) next #skip if no alleles (assume not because of dropout)
      indUse = which(MACmatrix[i,])
      cc = rowSums(MACmatrix[,indUse,drop=FALSE]) #count number of matching alleles
      mac[i,] <- mac[i,] + cc/nAlleles[i] #insert score
      nlocs[i,] <- nlocs[i,] + as.integer(nAlleles>0) #insert number of loci
    }
    if(nrow(allrefs[[loc]])>nLarge) print(paste0("Locus ",loc," calculated"))
  }
  mac = mac/nlocs #normalize wrt number of loci
  #print(mac)
  
  #get candidates based on MAC threshold
  x0 <- get("setupThresh",envir=nnTK)$MACthresh 
  ext = 1:(nR-1) #columns/rows to extract
  tabind = numeric()
  for(r in nR:2) { #looping through for each candidate:
    indcand1 = which(mac[r,ext]>=x0) #Lower triangle considered
    indcand2 = which(mac[ext,r]>=x0) #consider symmetry
    
    isSame = intersect(indcand1,indcand2) #get candidates that are same
    score1 = mac[r,ext[isSame]] #get score 1
    score2 = mac[ext[isSame],r] #get score 2
    
    insInd = numeric() #temporary adding indices
    if(length(isSame)>0) {
      insInd = cbind(r,isSame) #insert that score1 is largest (by default)
      insInd[score2>score1] = insInd[,2:1] #swap if not
    }
    rem1 = setdiff(indcand1,isSame) #get reminder1 
    rem2 = setdiff(indcand2,isSame) #get reminder2
    if( length(rem1)>0 )  insInd = rbind(insInd, cbind(r,rem1) ) #insert reminder1
    if( length(rem2)>0 )  insInd = rbind(insInd, cbind(rem2,r) ) #insert reminder2
    if( length(insInd)>0 ) tabind = rbind(tabind, insInd ) #add to matrix
    ext = ext[-length(ext)] #iteratively remove last index in extraction
  }
  
  if(length(tabind)>0) {
#    tab <- cbind(signif(mac[tabind],3),paste0(allrefsn[tabind[,1]]," - ",allrefsn[tabind[,2]]))
    tab <- cbind(signif(mac[tabind],3), allrefsn[tabind[,1]], allrefsn[tabind[,2]] )
    ord <- order(as.numeric(tab[,1]),decreasing=TRUE)
    out <- tab[ord,,drop=FALSE]
    return(out)
  } else {
    return(NULL)
  }
}
