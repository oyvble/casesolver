#' @title calcIBS
#' @description Calculate identical by state statistics between all references
#' @param nnTK an environment object from stored CaseSolver object (includes data)
#' @param nLarge A number indicating VERY large numeber of references
#' @param mixtureName A string indicating "mixture" (specific language)
#' @return A list with canidate results
#' @export

calcIBS = function(nnTK,nLarge=10000,mixtureName="mixture") {  #Function to calculate IBS between references
  DBref <- get("refDataTABLE",envir=nnTK) #consider lists
  if(nrow(DBref)>=nLarge) {
    cat(paste0("The number of references were massive (>",nLarge,").\n Only references given in MatchStatus are considered!\n"))
    cand <- get("mixDataMATCHSTATUS",envir=nnTK) #get matchstatus lists
    cand = unique(cand[cand!=mixtureName]) #consider only refs which fits to single sources 
    DBref = DBref[rownames(DBref)%in%cand,]
  }
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
    isOne <- !(isna | grepl("/",tmp)) #index of one allele
    isTwo <- !(isna | isOne) #index of two alleles
    amat <- matrix("",nrow=nR,ncol=2)
    if(any(isTwo)) amat[isTwo,] <- t(matrix(unlist(strsplit( tmp[isTwo],"/")),nrow=2)) #gets error if none have two alleles
    amat[isOne,1] <- tmp[isOne] #insert only one allele
    #swap = amat[,2]<amat[,1] #sorting already done
    #amat[swap,] = amat[swap,2:1]
    allrefs[[loc]] <- amat
  }
  rm(amat,tmp);gc()
  
  if(require(bigmemory)){ #install.packages("bigmemory")
    options(bigmemory.typecast.warning=FALSE)
    mac = bigmemory::big.matrix(init=0,nrow = nR, ncol = nR, type = "char") # use short?
  } else { #if bigmemory not found
    tryCatch({
      mac = matrix(0,nrow = nR, ncol = nR) #Ordinary matrix
    },error = function(e)  
      print("Please install the R-package bigmemory to handle large memory and try again!")   
    )
  } 
  #  colnames(mac) <- rownames(mac) <- allrefsn
  locs <- names(allrefs) #get non-empty loci
  
  for(loc in locs) {
    amat <- allrefs[[loc]]
    unmat = unique(amat) #get only unique
    for(i in 1:nrow(unmat)) { #for each reference
      if(unmat[i,1]=="") next #skip if no alleles
      t1 <- unmat[i,1]==amat[,1] | unmat[i,1]==amat[,2]
      t2 <- unmat[i,2]==amat[,1] | unmat[i,2]==amat[,2]
      insval <- as.numeric(t1) #value to insert
      if(unmat[i,2]!="") insval <- insval + as.numeric(t2) #add value if two alleles
      if( unmat[i,1]==unmat[i,2] ) { #if reference had homozygout genotype
        ind <- amat[,1]!=amat[,2] & insval==2 #those not same but got MAC=2
        insval[ind] <- insval[ind] - 1 #subtract 1
      }
      #    mac[insind,] <- mac[insind,] + t(replicate(sum(insind),insval)) #insert
      insind2 = insval>0 #only indices which is positive
      insind =  unmat[i,1]==amat[,1] & unmat[i,2]==amat[,2] #sum(insind)
      mac[insind2,insind] <- mac[insind2,insind] + as.integer(insval[insind2]) #insert
    }
    if(nrow(amat)>nLarge) print(paste0("Locus ",loc," calculated"))
  }
  #print(mac)
  
  #get candidate
  x0 <- get("setupThresh",envir=nnTK)$minIBS #max(as.numeric(names(tab1)))
  ext = 1:(nR-1) #columns to extract
  tabind = numeric()
  for(r in nR:2) { #looping through for each candidate:
    indcand = which(mac[r,ext]>=x0) #extract those columns which had above threshold 
    ext = ext[-length(ext)] #remove last index in extraction
    if( length(indcand)>0) tabind = rbind(tabind, cbind(r,indcand) ) #add to matrix
    
  }
  if(length(tabind)>0) {
#    tab <- cbind(mac[tabind],paste0(allrefsn[tabind[,1]]," - ",allrefsn[tabind[,2]]))
    tab <- cbind(mac[tabind], allrefsn[tabind[,1]], allrefsn[tabind[,2]])
    ord <- order(as.numeric(tab[,1]),decreasing=TRUE)
    out <- tab[ord,,drop=FALSE]
    return(out) #return result table
    
  } else {
    return(NULL)
  }
}
