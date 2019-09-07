#' @title tabToListRef
#' @description Function to convert a matrix to list format. 
#' @details Used to convert reference data stored in matrix to list format required by functions.
#' @param tab A (nR x nL) matrix with rows giving locus specific genotypes (for each column) for references
#' @param ln Loci names of return list  
#' @param forceDI Whether genotypes are restricted to contain 2 alleles or not.
#' @export

tabToListRef  = function(tab,ln=NULL,forceDi=TRUE) {
   locs = colnames(tab)

   if(!is.null(ln)) {
     rmLoc = locs[!locs%in%ln] #locs to be removed
     adLoc = ln[!ln%in%locs] #locs to be added (as empty)
     if(length(rmLoc)>0) print(paste0("Loci to be removed: ", paste0(rmLoc,collapse="/")))
     if(length(adLoc)>0) print(paste0("Loci to be added (as empty): ", paste0(adLoc,collapse="/")))
     locs <- intersect(locs,ln)
   }
   
   refn = rownames(tab)
   refL = list()
   for(ref in refn) { #for each reference
     refL[[ref]] = list()
     indref = which(ref==refn)
     for(loc in locs) { #for each locus
       av =  unlist(strsplit(tab[indref, colnames(tab)==loc],"/"))
       av = av[!is.na(av)] #remove NAs
       if(forceDi && length(av)==1) av = rep(av,2)
       refL[[ref]][[loc]] = list(adata=av)
     }
   }
   return(refL)
 }
 