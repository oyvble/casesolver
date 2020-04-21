#' @title getStructuredData
#' @description A function for structuring data 
#' @details The function must return in a specific format.
#' @param data a list with a evidence table (mix) and a reference table (ref)
#' @param ln markers to consider (a given order)
#' @param minLoc is minimum number of overlapping markers to be identical (for partial profiles).
#' @export
getStructuredData = function(data,ln,minLoc=10) {
   #data is a list with a evidence table (mix) and a reference table (ref):
   #ln the markers to consider (a given order)
   #popFreq is a list with allele freqs popFreq[[locus]]. NB: Must contain same loci as in ln
   #minLoc is minimum number of overlapping markers to be identical (for partial profiles).
   #Evidences: colnames(data$mix)=("SampleName","Markers","Alleles","Heights"). Separated with "/" in text format.
   #References: colnames(data$ref)=("SampleName","Markers","Alleles"). Separated with "/" in text format.
   if(nrow(data$mix)>0) data$mix[,2] <- toupper(data$mix[,2]) #LOCUS-names are assigned as Upper-case! This is important to do!
   if(nrow(data$ref)>0) data$ref[,2] <- toupper(data$ref[,2]) #LOCUS-names are assigned as Upper-case! This is important to do!

   sortGenotypes = function(tmp) {
    suppressWarnings({
     v1 <- as.numeric(tmp[,1])
     v2 <- as.numeric(tmp[,2])
     })
     isNA <- is.na(v1) | is.na(v2)
     swap <- tmp[,2]<tmp[,1]  #sort alleles by text
     swap[!isNA] <-  v2[!isNA]<v1[!isNA]  #sort alleles by number
     tmp[swap,1:2] <- tmp[swap,2:1]
     return(tmp)
   }

   #Fix ref-alleles: Insert single alleles, order alleles in ref by increasing order
   if(nrow(data$ref)>0) {
    isempty <- data$ref[,3]==""
    impute <- !isempty & !grepl("/",data$ref[,3]) #have only one allele. Must get correct format: 2 alleles
    data$ref[impute,3] <-  paste0(data$ref[impute,3],"/",data$ref[impute,3]) #include twice 
    ord <- !isempty & !impute 
    tmp <- unlist( strsplit( data$ref[ord,3],"/") )
    if(length(tmp)!=2*sum(ord)) { #don't allow more than two alleles
      stop("The number of alleles for some of your references was greater than two! Please make sure not to exceed two.")
    }
    tmp <- t(matrix(tmp ,nrow=2 )) #get genotypes for each markers
    tmp <- sortGenotypes(tmp) #get sorted genotypes
    data$ref[ord,3] <- paste0(tmp[,1],"/",tmp[,2]) #make sure it get correct order
   } 

   #FIRST: DUPLICATED REF-PROFILES WILL BE REMOVED, GIVEN A NOTIFICATION AS MESSAGE
   #SECOND: REF-PROFILES EQUAL EVID-PROFILES (SINGLE) WILL BE DETECTED, WILL BE NOTIFIED BEH
   #THIRD: SINGLE EVID-PROFILES WILL BE ASSIGNED AS AN UNKNOWN (1,2,3 etc)
   refNames = evidNames = NULL
   if(nrow(data$ref)>0) refNames <- unique(data$ref[,1]) #vector of unique ref-samples
   if(nrow(data$mix)>0) evidNames <- unique(data$mix[,1]) #vector of unique evid-samples

   stringRef <- rep("",length(refNames))
   stringEvid <- rep("",length(evidNames))
   nClow <- rep(0,length(evidNames)) #used to get lower boundary of number of contributors in evidence profiles
   for(loc in ln) { #go through each loci and create a string per sample profiles: Used to detect similarity
#loc = ln[1]
    subRef <- data$ref[data$ref[,2]==loc,,drop=F]
    subRef <- subRef[match(refNames,subRef[,1]),,drop=F] #get correct order of data (NECESSARY IF IMPORT NOT CONSISTENT)
    subEvid <- data$mix[data$mix[,2]==loc,,drop=F]
    subEvid <-  subEvid[match(evidNames ,subEvid[,1]),,drop=F] #get correct order of data (NECESSARY IF IMPORT NOT CONSISTENT)
    subEvid[is.na(subEvid[,3]),3] <- "NA" #add string NA
    subEvid[subEvid[,3]=="",3] <- "NA" #add string NA
    subRef[is.na(subRef[,3]),3] <- "NA" #add string NA
    subRef[subRef[,3]=="",3] <- "NA" #add string NA

    #add loci separation sign before adding alleles to vector 
    if(which(loc==ln)>1) {
      if(length(refNames)>0) stringRef <- paste0(stringRef,"-") 
      if(length(evidNames)>0) stringEvid <- paste0(stringEvid,"-")  #add loci separation sign 
    }
    stringRef <- paste0(stringRef, subRef[,3]) #the alleles of ref is inserted (notice the order)  
    stringEvid <- paste0(stringEvid, subEvid[,3]) #the alleles of ref is inserted (notice the order)  

    if(length(evidNames)>0) { #require at least one evidence
     nCtmp <- ceiling(sapply(strsplit(subEvid[,3],"/"),length)/2)
     nCtmp[subEvid[,3]=="NA"] <- 0 #set to zero if na
     ind <- nCtmp>nClow #evids to update
     nClow[ind] <- nCtmp[ind]
    }
   } #end for each loci
   #Fix single source strings:
   evidNamesSS <- evidNames[nClow==1] #consider only single source profiles
   stringEvidSS <- stringEvid[nClow==1]  #consider only single source profiles

   #Fix stringEvid2: Impute homozygous alleles, and order the alleles:
   if(length(evidNamesSS)>0) { #require at least 1 sample
    stringEvidSS2 <- t(matrix( unlist(strsplit(stringEvidSS,"-")),ncol=length(stringEvidSS) )) #get genotypes for each markers  
    isempty <- stringEvidSS2=="NA" 
    impute <- !isempty & !grepl("/",stringEvidSS2) #have only one allele. Must get correct format
    stringEvidSS2[impute] <-  paste0(stringEvidSS2[impute],"/",stringEvidSS2[impute]) #include alleles twice 
    ord <- !isempty & !impute #markers with heteroyzogous variants
    tmp <- unlist( strsplit( stringEvidSS2[ord],"/") )
    tmp <- t(matrix(tmp ,nrow=2 )) #get genotypes for each markers
    tmp <- sortGenotypes(tmp) #get sorted genotypes
    stringEvidSS2[ord] <- paste0(tmp[,1],"/",tmp[,2]) #

    #update the evid-strings stringEvidSS and stringEvid again:
    stringEvidSS <- numeric() #restart this
    for(loc in ln) { #go through each loci and create a string per sample profiles: Used to detect similarity
     ind <- which(loc==ln) 
     if(ind>1) {
      stringEvidSS <- paste0(stringEvidSS,"-")  #add loci separation sign 
     }
     stringEvidSS<- paste0(stringEvidSS, stringEvidSS2[,ind]) #the alleles of ref is inserted (notice the order)  
    }
    stringEvid[nClow==1] <- stringEvidSS #update  evid-string
   } #end for SS samples

   #DETECT SIMILAR REFERENCES:
   #evidMatch is a vector with identical=match info ("mixture","refname","unknown X")
   evidMatch <- rep("mixture",length(evidNames)) 

   if(length(refNames)>0) {#IF ANY REFERENCES
    refU <- unclass(factor(stringRef)) #get similar refs (each number is a unique reference)
    isdup <- duplicated(refU) #which are duplicates?
    refDups <- refNames[isdup] #these refs are duplicates
    data$ref <- data$ref[!data$ref[,1]%in%refDups,] #update refs such that non-duplicates are keeped
    if( length(refDups)>0 ) {
     dups <- refU[isdup] #these are duplicate grps
     dupsUn <- unique(dups) #get unique duplicate grps
     for(gg in dupsUn) {
      print(paste0("Equal references detected: ",  paste0(refNames[refU==gg],collapse="="), ". Keeping first."))
     }
    }
    refNames <- na.omit(refNames[!isdup]) #unique refs: refnames
    stringRef <- stringRef[!isdup] #unique refs: refalleles
    stringRef2 <- t(matrix( unlist(strsplit(stringRef,"-")),ncol=length(stringRef) )) #get genotypes for each markers  
    nR = length(stringRef) #get references

    #DETECT SIMILAR EVIDENCE AND REF - PROFILES: 
    #SHOULD BE SIMILAR TO EFFICIENT MAC COMPARISON
    if(length(evidNamesSS)>0) { #require at least 1 sample
     unStringEvidSS <- unique(stringEvidSS) #get unique SS-string 
     stringEvidSS2 <- unique(stringEvidSS2) #get corresponding uniques
     nS = length(unStringEvidSS) #get unique SS samples 
     print(paste0("Comparing references against SingleSource samples: ",nS*nR," comparisons."))
     matchvec <- rep(NA,nS) #create a vector for finding matches
     for(ss in 1:nS) { #for each sample: We calculate pairwise comparison
       #ss=1
       macS <- nLocs <- rep(0,nR) #make vector for all references 
       for(loc in ln) { #for each locus: Vectorizing as much as possible!
        ll=which(loc==ln)
        if(is.na(stringEvidSS2[ss,ll]) || stringEvidSS2[ss,ll]=="NA")  next #skip if no data
        sttmp <- unlist(strsplit(stringEvidSS2[ss,ll],"/")) #get alleles
        isna <- is.na(stringRef2[,ll]) | stringRef2[,ll]=="NA"  #get which refs are non-zero. Not counted if missing!
        numWithin <- rep(0,sum(!isna)) #number of unique alleles of reference that are in stain
        if(length(numWithin)==0) next #skip if no references
        AvecList <- t(matrix(unlist(strsplit(stringRef2[!isna,ll],"/")),nrow=2)) #get alleles of refs (in case of 2)
        for(aa in unique(sttmp)) numWithin <- numWithin + rowSums(AvecList==aa) #sum up for each alleles in stain 
	  
 	   if( sttmp[1]==sttmp[2] ) { #if SS had homozygout genotype
     	ind <- AvecList[,1]!=AvecList[,2] & numWithin==2 #those not same but got MAC=2
     	numWithin[ind] <- numWithin[ind] - 1 #subtract 1 because these were not homozygous 
    	   }
        macS[!isna] <- macS[!isna] + numWithin #add number of matching alleles
        nLocs[!isna] <- nLocs[!isna] + 1  #add locus  
    	  } #end for each locus
       missmatch = 2*nLocs - macS #get number of missmatches
       foundind = which(missmatch==0 & nLocs>=minLoc ) #index with zero mismatches AND EVALUATED AT LEAST minLoc
       if(length(foundind)>0) { #if found a match 
	   evidMatch[unStringEvidSS[ss]==stringEvid] = refNames[ foundind[which.max(nLocs[foundind])] ]
        if(length(foundind)>1) {
         print("Several references matched (SIMILAR REFERENCES). Keeping the one with most loci!")
        }
       } #else if no match found
      } #end for each SS sample
      rm(AvecList ,numWithin,missmatch ) ; gc()
    } #end if having SS samples
   } #end IF ANY REFS

   #sort structures: Hierarchical differences 
   evidord <- order(nchar(stringEvid)) #sort evid with respect to number of characters
   stringEvid <- stringEvid[evidord] #update order
   evidNames <- evidNames[evidord] #update order
   evidMatch <- evidMatch[evidord]
   names(evidMatch) <- evidNames  #get vector of match status for all evidence profiles
   nClow <- nClow[evidord]  #update order
  
   #CREATE STRUCTURE: data-frame tables + lists
   dfmix <- matrix(ncol=length(ln),nrow=length(evidNames)) 
   dfref <- matrix(ncol=length(ln),nrow=length(refNames ))
   colnames(dfref) <- colnames(dfmix) <- ln
   rownames(dfmix) <- evidNames
   rownames(dfref) <- refNames 
   mixL <- list()  #refL is not considered anymore

   #Structure Evidences into list and matrix
   for(mix in evidNames) { #for each evid profile: CAN BE TIMECONSUMING FOR MANY EVIDS!
#mix = evidNames[1]
      mixL[[mix]] = list() 
      subEvid <- data$mix[ data$mix[,1]==mix,,drop=F] #get evid data
      subEvid <- subEvid[match(ln,subEvid[,2]),,drop=F] #get correct order of loci  (NECESSARY IF IMPORT NOT CONSISTENT)
      subEvid[is.na(subEvid[,3]),3:4] <- "" #add strings "" if NA
      tryCatch({ 
        dfmix[which(evidNames==mix),] = subEvid[,3] #extract import data
      }, error = function(e) print( paste0("There was an inconsistent number of loci observations for evidence ",mix))  ) 
      adata = strsplit(subEvid[,3],"/")
      hdata = strsplit(subEvid[,4],"/")
      for(loc in ln) { #loop for each loci
       locind = which(loc==ln)
       mixL[[mix]][[loc]] = list(adata=adata[[locind]],hdata=as.numeric(hdata[[locind]]))
      }
   }
 
   #helpfunction of whether a locus is a Y-marker
   isY = function(x) return( toupper(substr(x,0,2))=="DY" || toupper(substr(x,0,1))=="Y" )
   
   #Structure References into matrix
   for(loc in ln) { #
#loc=ln[26]
      subRef <- data$ref[data$ref[,2]==loc,,drop=F]
      subRef <- subRef[match(refNames,subRef[,1]),,drop=F] #get correct order of data (NECESSARY IF IMPORT NOT CONSISTENT)
      subRef[is.na(subRef[,3]),3] <- "" #add string "" if NA
	  #HANDLING Y-markers for references:
	 if( isY(loc) ) { #check if Y-marker
	        consInd = which(subRef[,3]!="") #indices to consider
		    if(length(consInd)>0) {
				tmp = matrix(unlist(strsplit(subRef[consInd,3],"/")),nrow=2) #get all genotypes
				changeInd = tmp[1,]==tmp[2,] #homozygous variants are changed
				if(any(consInd)) subRef[consInd[changeInd],3] <- tmp[1,changeInd] #insert only 1 allele variant
			}		   
     } #end if Y-marker

      tryCatch({
        dfref[,which(ln==loc)] = subRef[,3] #extract import data
      }, error = function(e) print( paste0("There was an inconsistent number of reference observations at locus ",loc))  )
   }
   evidMatch[evidMatch=="mixture" & nClow==0] <- "empty" #set empty if no data
   #evidMatch[nClow==1] <- "mixture" #reset to test next block

   #CONCERNING UNKNOWN PROFILES DEDUCED FROM SINGLE SOURCE PROFILES
   #Assign unknowns for the rest:
   isUnknown <- evidMatch=="mixture" & nClow==1 #potential unknowns
   if(!is.na(isUnknown[1]) && any(isUnknown)) { #continue only if any potential unknowns
    #Need to check if profiles may be same unknowns (because of NA):
    #Must compare between all samples.
    stringEvidSSun <- unique(stringEvid[isUnknown]) #get uniques strings
    stringEvidSSun2 <- matrix(unlist(strsplit(stringEvidSSun,"-")),ncol=length(stringEvidSSun)) #get uniques strings
    nS <- ncol(stringEvidSSun2) #number of unknowns
    if(nS>1) { #needs at least 2 unknowns
     matCC <- matrix(ncol=3,nrow=nS*(nS-1)/2) #a list of match connections
     cc <- 1 #counter
#FOLLOWING PROCEDURE CAN BE VERY TIME CONSUMING FOR LARGE NUMBER OF UNKNOWNS
     for(evid1 in 1:(nS-1)) { #for each unique evidence: 
      for(evid2 in (evid1+1):nS) { #for each unique evidence
       evid1val <- stringEvidSSun2[,evid1]
       evid2val <- stringEvidSSun2[,evid2]
       evid1val[evid2val=="NA"] <- "NA" #allow missing genotypes in evid1
       evid2val[evid1val=="NA"] <- "NA" #allow missing genotypes in evid2
       nLocs <- sum(evid1val!="NA") #number of non-empty loci
       match <- all(evid1val==evid2val) & nLocs>=minLoc #must match + have sufficiently number of loci
       matCC[cc,] <- c(evid1,evid2,match)
       cc <- cc + 1
      }
     }
    } 

    #Label unknown profiles and create a data matrix (reflist not considered anymore)
    unknownInd <- rep(NA,nS)
    refU <- matrix(nrow=nS,ncol=length(ln)) #unknown data matrix 
    colnames(refU) = ln 
    for(evid1 in 1:nS) { #traverse through each evidence
      matchevid <- numeric() 
      if(nS>1) {
       subind <- matCC[,1]==evid1 | matCC[,2]==evid1
       subCC <- matCC[subind,,drop=FALSE]
       matchevid <- subCC[subCC[,3]==1,2] #find potential other matching evidence
      }
      unk <- length(na.omit(unique(unknownInd)))+1 #unknown tag
      if(is.na(unknownInd[evid1])) {
        unknownInd[evid1] <- unk
        unknownInd[matchevid] <- unknownInd[evid1] #copy name to matching 

        #Create concensus of unknown profile: 
        for(loc in ln) {
         locind = which(loc==ln)
         ins <- stringEvidSSun2[locind,c(evid1,matchevid)]
         av <- ins[ins!="NA"][1] #unlist(strsplit(ins[ins!="NA"][1],"/")) #get alleles (should be able to edit this)
         if(!is.na(av[1]) && isY(loc) ) { #check if Y-marker
           av2 = unlist(strsplit(av,"/"))
           if(av2[1]==av2[2]) av = av2[1] #given as one allele if marker name starts if Y-marker
         }
         refU[unk,locind] <- av
        } #end for each marker
      } #end if new unknown
    } #end for each evidence
    nU = max(unknownInd)
    unknname <- paste0("Unknown ",1:nU)
    refU = refU[1:nU,,drop=FALSE] #keep only data table with inserted unknowns
    rownames(refU) = unknname 

    #NOTE: CAN SOME refs may be IDENTICAL because loci-information is filled in?
    if(any(duplicated(refU))) {
 	 print("NOTE: Some unknown references are redundant. This is not taken into account!")
    }

    #Add unknowns to evidMatch:
    evidUnUn <- unknname[unknownInd] #name of the unknowns belonging to SSevids
    for(str in stringEvidSSun) { #for each unique string
     ind <- which(stringEvid==str) #get index of equal strings
     evidMatch[ind] <- evidUnUn[which(str==stringEvidSSun)] #insert unknown name
    }

    #Add unknowns (in refList) to refDataTABLE
    dfref <- rbind(dfref,refU) #add to other references
   } #End if any unknowns

   return( list(mixDataTABLE=dfmix, refDataTABLE=dfref,mixDataLIST=mixL,mixDataMATCHSTATUS=evidMatch))
 } #end function getStructuredData 
