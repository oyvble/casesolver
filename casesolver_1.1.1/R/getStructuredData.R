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
   dataRef <- data$ref #create a copy (used to avoid problem with single alleles for refs)

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
    impute <- !isempty & !grepl("/",data$ref[,3]) #have only one allele. Must get correct format
    data$ref[impute,3] <-  paste0(data$ref[impute,3],"/",data$ref[impute,3]) #include twice 
    ord <- !isempty & !impute 
    tmp <- unlist( strsplit( data$ref[ord,3],"/") )
    if(length(tmp)!=2*sum(ord)) { #don't allow more than two alleles
      stop("The number of alleles for some of your references was greater than two! Please make sure not to exceed two.")
    }
    tmp <- t(matrix(tmp ,nrow=2 )) #get genotypes for each markers
    tmp <- sortGenotypes(tmp) #get sorted genotypes
    data$ref[ord,3] <- paste0(tmp[,1],"/",tmp[,2]) #
   } 

   #FIRST: DUPLICATED REF-PROFILES WILL BE REMOVED, GIVEN A NOTIFICATION AS MESSAGE
   #SECOND: REF-PROFILES EQUAL EVID-PROFILES (SINGLE) WILL BE DETECTED, WILL BE NOTIFIED BEH
   #THIRD: SINGLE EVID-PROFILES WILL BE ASSIGNED AS AN UNKNOWN (1,2,3 etc)
   refNames <- unique(data$ref[,1]) #vector of unique ref-samples
   evidNames <- unique(data$mix[,1]) #vector of unique evid-samples
   stringRef <- rep("",length(refNames))
   stringEvid <- rep("",length(evidNames))
   nClow <- rep(0,length(evidNames)) #used to get lower boundary of number of contributors in evidence profiles
   for(loc in ln) { #go through each loci and create a string per sample profiles: Used to detect similarity
    subRef <- data$ref[data$ref[,2]==loc,,drop=F]
    subEvid <- data$mix[data$mix[,2]==loc,,drop=F]
    subEvid[subEvid[,3]=="",3] <- "NA" #add string NA
    subRef[subRef[,3]=="",3] <- "NA" #add string NA

    #insert alleles:
    refNew <- rep("NA",length(refNames)) #insert ref alleles
    evidNew <- rep("NA",length(evidNames)) #insert evid alleles
    ordRef <- match(subRef[,1],refNames) #assure that orders are consistent for the two vectors
    ordEvid <- match(subEvid[,1],evidNames) #assure that orders are consistent for the two vectors
    if(nrow(subRef)>0) refNew[ordRef] <- subRef[,3] #reftmp2 
    evidNew[ordEvid] <- subEvid[,3]

    #add loci separation sign before adding alleles to vector 
    if(which(loc==ln)>1) {
      stringRef <- paste0(stringRef,"-") 
      stringEvid <- paste0(stringEvid,"-")  #add loci separation sign 
    }
    stringRef <- paste0(stringRef, refNew) #the alleles of ref is inserted (notice the order)  
    stringEvid <- paste0(stringEvid, evidNew ) #the alleles of ref is inserted (notice the order)  

    nCtmp <- rep(0,length(evidNames))
    nCtmp[ordEvid] <- ceiling(sapply(strsplit(subEvid[,3],"/"),length)/2)
    nCtmp[ordEvid][subEvid[,3]=="NA"] <- 0 #set to zero if na
    ind <- nCtmp>nClow #evids to update
    nClow[ind] <- nCtmp[ind]
   }
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
   }

   #DETECT SIMILAR REFERENCES:
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

   #DETECT SIMILAR EVIDENCE AND REF - PROFILES: 
   #evidMatch is a vector with identical=match info ("mixture","refname","unknown X")
   evidMatch <- rep("mixture",length(evidNames)) 
   if(length(evidNamesSS)>0) { #require at least 1 sample
    unStringEvidSS <- unique(stringEvidSS) #get unique SS-string 
    stringEvidSS2 <- unique(stringEvidSS2) #get corresponding uniques
    if(length(unStringEvidSS)>0) {
     for(ref in refNames) { #for each unique ref
      for(evid in 1:length(unStringEvidSS)) { #for each unique evidence
       stringRef2 <- unlist(strsplit(stringRef[ref==refNames],"-")) #get string of ref (temp var)
       stringEvid3 <- stringEvidSS2[evid,] #temp variable
       stringEvid3[stringRef2=="NA"] <- "NA" #allow missing genotypes in reference
       stringRef2[stringEvid3=="NA"] <- "NA" #allow missing genotypes in evidence
   
       nLocs <- sum(stringEvid3!="NA") #number of non-empty loci
       match <- all(stringRef2==stringEvid3) & nLocs>=minLoc #must match + have sufficiently number of loci
       if(match) evidMatch[unStringEvidSS[evid]==stringEvid] <- ref #insert name of ref to matching evid profiles  
      }
     }
    } #end if having unique single sources
   } #end if having at least 1 sample
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
   mixL <- refL <- list() #mix and ref-list

   for(mix in evidNames) { #for each evidence profile
     sub <- data$mix[data$mix[,1]==mix,,drop=F] #extract import data
     dfmix[which(mix==evidNames),which(ln%in%sub[,2])] <- sub[,3] #exctract only alleles
     #for each loci: insert adata and hdata
     mixL[[mix]] <- list()
     for(loc in ln) {
       dat <- strsplit(sub[loc==sub[,2],3:4],"/")
       mixL[[mix]][[loc]] <- list(adata=dat$Alleles,hdata=as.numeric(dat$Heights))
     }
   }
   for(ref in refNames ) { #for each ref profile
#     sub <- data$ref[data$ref[,1]==ref,,drop=F] #extract import data
     sub <- dataRef[dataRef[,1]==ref,,drop=F] #extract import data
     locind <- match(sub[,2],ln) #get position of loci
     if(any(duplicated(locind))) {
       print(sub) #show data
       stop(paste0("The reference had inconstistent genotypes for the same marker. Please fix data and try again!"))
     }
     dfref[which(ref==refNames),locind] <- sub[,3] #exctract only alleles

     #for each loci: insert adata
     refL[[ref]] <- list()
     for(loc in ln) {
       dat <- strsplit(sub[loc==sub[,2],3],"/")
       refL[[ref]][[loc]] <- list(adata=unlist(dat))
     }
   }
   evidMatch[evidMatch=="mixture" & nClow==0] <- "empty" #set empty

  #CONCERNING UNKNOWN PROFILES
   #Assign unknowns for the rest:
   isUnknown <- evidMatch=="mixture" & nClow==1 #potential unknowns
   if(any(isUnknown)) { #continue only if any potential unknowns
    #Need to check if profiles may be same unknowns (because of NA):
    #Must compare between all samples.
    stringEvidSSun <- unique(stringEvid[isUnknown]) #get uniques strings
    stringEvidSSun2 <- matrix(unlist(strsplit(stringEvidSSun,"-")),ncol=length(stringEvidSSun)) #get uniques strings
    nS <- ncol(stringEvidSSun2)
    if(nS>1) { #needs at least 2 unknowns
     matCC <- matrix(ncol=3,nrow=nS*(nS-1)/2)
     cc <- 1 #counter
     for(evid1 in 1:(nS-1)) { #for each unique evidence
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

    #Label unknown profiles and add to reflist also:
    unknownInd <- rep(NA,nS)
    refU <- list() #unknown list
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
        refU[[unk]] <- list()
        for(loc in ln) {
         ins <- stringEvidSSun2[which(loc==ln),c(evid1,matchevid)]
         av <- unlist(strsplit(ins[ins!="NA"][1],"/")) #get alleles (should be able to edit this)
         if( is.na(av[1]) ) av <- numeric() #empty marker
         if(length(av)==1) av = rep(av,2)
         refU[[unk]][[loc]] <- list(adata=av)
        } #end for each marker
      } #end if new unknown
    } #end for each evidence

    #NOTE: SOME refL may be IDENTICAL because loci-information is filled in.
    #Idea: Check all unknown reflists against each others, reveal identical profiles (i.e. it is same unknown). Must take into account missing alleles
    unknname <- paste0("Unknown ",1:max(unknownInd))
    names(refU) <- unknname
    nU <- length(refU) 
    if(nU>1) { #check only if more than 1
     for(u in 1:(nU-1)) {
      for(v in (u+1):nU) { 
 	   u1 <- unlist(refU[[unknname[u]]])
	   u2 <- unlist(refU[[unknname[v]]])
       match <- length(u1)==length(u2) && all(u1==u2) #must have same number of alleles and same alleles
       if(match) {
         refU[[unknname[v]]] <- NULL #Set second to NULL
         unknownInd[unknownInd==v] <- u
       }
      }
     }
    } #end for more than 1 unknown
    unknname <- paste0("Unknown ",1:length(refU)) #update names
    refL[unknname] <- refU #add unknown refs to ref-list

    #Add unknowns to evidMatch:
    evidUnUn <- unknname[unknownInd] #name the unknowns
    for(str in stringEvidSSun) { #for each unique string
     ind <- which(stringEvid==str) #get index of equal strings
     evidMatch[ind] <- evidUnUn[which(str==stringEvidSSun)] #insert unknown name
    }

    #Add unknowns (in refList) to refDataTABLE
    DBref2 <- numeric()
    for(ref in unknname) { #for each unknowns
      tmp  <- sapply(refU[[ref]], function(x) {
       av <- x$adata
       return(paste0(av,collapse="/"))
      })
      DBref2 <- rbind(DBref2 , tmp)
    } #end for each unknown
    rownames(DBref2) <- unknname
    dfref <- rbind(dfref,DBref2) #add to other references
   } #End if any unknowns

   return( list(mixDataTABLE=dfmix, refDataTABLE=dfref,mixDataLIST=mixL,refDataLIST=refL,mixDataMATCHSTATUS=evidMatch))
 } #end function getStructuredData 
