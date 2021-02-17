#' @title clusterOrder
#' @description Order samples wrt to clustered hierachical structure (similarity)
#' @details Create distance between all genotypes (over all markers) and perform hierachical clustering (order returned)
#' DISTANCE OUTCOMES:
#' All shared (dist=0): (A/B vs A/B), (A/A vs A/A), (A vs A/B), (A vs A/A), (A vs A)
#'None shared (dist=2): (A/B vs C/D), (A/A vs B/B), (A vs B/C), (A vs B/B), (A vs B)
#' One shared (dist=1): (A/B vs B/C), (A/A vs A/B)
#' @param table Table with genotypes
#' @param plot Whether the dendrogram should be plotted
#' @return order vector (ranked wrt hierarchical clustering)
#' @export

clusterOrder = function(table,plot=FALSE) {
  nRows = nrow(table) #number of samples
  nCols = ncol(table) #number of loci
  if(nRows==1 || nCols==0) return(NULL)
    
  DIST = matrix(0,nRows,nRows) #rep(0,nRows*(nRows-1)) #lower tringular distance
  for(col in 1:nCols) {
    genos = table[,col] #obtain all deconvolved genotyes (are above DCthresholds)
    
    unGenos = unique(na.omit(genos)) #get only unique genotypes (zero not considered)
    nUnique = length(unGenos) #number of unique genotypes 
    if(nUnique==0) next #skip if none
    backIndex = match(genos,unGenos) #obtain indices to insert back
    # unGenos[backIndex]==genos
    
    isTwo <- grepl("/",unGenos) #index of two alleles
    isOne <- !isTwo #index of one allele
    amat <- matrix("",nrow=nUnique,ncol=2) #allele matrix
    if(any(isTwo)) amat[isTwo,] <- t(matrix(unlist(strsplit( unGenos[isTwo],"/")),nrow=2)) #gets error if none have two alleles
    amat[isOne,1] <- unGenos[isOne] #insert only one allele
    isHet =  amat[,1]!=amat[,2]  #indicate het genos
    
    for(i in 1:nUnique) { #traverse each unique genotype
      #   print(unGenos[i])
      VAL1 <- as.numeric(amat[i,1]==amat[,1] | amat[i,1]==amat[,2])
      mac = VAL1
      if(isTwo[i]) {
        VAL2 <- as.numeric(amat[i,2]==amat[,1] | amat[i,2]==amat[,2])
        mac <- mac + VAL2 #add value if two alleles
        if( !isHet[i] ) { #if reference had homozygout genotype
          ind <- isHet & mac==2 #those not same but got MAC=2
          mac[ind] <- mac[ind] - 1 #subtract 1
        }
        mac[isOne] = 2*mac[isOne] #all matching ones to have score=2
      } else {
        mac = 2*mac #Ensure all matching to have score=2
      }
      #         names(mac) = unGenos
      
      insval = 2 - mac #obtain distance score = 2 - mac
      insvalAll = insval[backIndex] #obtain value for all samples
      insind1 = which( unGenos[i]==genos) #obtain index having the particular compared genotype
      insind2 = which(insvalAll>0) #only indices which is positive (ensures don't inserting for same genotype)
      DIST[insind2,insind1] <- DIST[insind2,insind1] + as.integer(insvalAll[insind2]) #insert to matrix
    }
    #       colnames(DIST) <- rownames(DIST) <- genos
  } #end for each marker     
  #       View(DIST);
  if(!isSymmetric(DIST)) {
    print("WARNING: DISTANCE CALCULATION WAS NOT SYMMETRIC!")
  }
  hc <- hclust(as.dist(DIST))# init distance matrix and perform  hierarchical cluster
  if(plot)  plot(hc)
  mergeInfo = t(hc$merge) #transpose to obtain pairs as rows
  ord = abs(mergeInfo[mergeInfo<0]) #obtain final order (turn negatives to positive index)
  
  return(ord) #returned ordering of samples
}# end if clustering