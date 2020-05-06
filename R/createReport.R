#' @title createReport
#' @description A function for creating report based on CaseSolver object
#' @details The function returns the names of the selected profiles (when window is closed)
#' @param nnTK an environment object from stored CaseSolver object
#' @export

createReport = function(nnTK) { #this function loads data and put them in a HTML script
  LUSsymbol <- "_"
  colonsymbol = ":" #use variable for colon
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  # L = getLanguage( "English")
  
  #helpfunction when extracting reference data
  getRefL = function(refs,forceDi=FALSE) { #return list with same order as for refs
    refL <- casesolver::tabToListRef(tab=refDataTABLE[ match(refs,rownames(refDataTABLE)),,drop=FALSE],forceDi=forceDi) #FORCING DUP alleles
    return(refL)
  }
  
  optReport <- get("setupReport",envir=nnTK)  #get report options
  checked = optReport$checked
  reportitems = optReport$reportitems
  radiotxt = c( L$on , L$off )
  
  #Provide data to use in report: May be modified if user deselects samples!
  mixLIST <- get("mixDataLIST",envir=nnTK) #get mix list
  mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #get data table from nnTK-environment
  refDataTABLE <- get("refDataTABLE",envir=nnTK) #get data table from nnTK-environment
  mixDataMATCHSTATUS <- get("mixDataMATCHSTATUS",envir=nnTK) #assign to nnTK-environment
  metaDataList <- get("metaDataLIST",envir=nnTK) #list of imported metadata
  consDataTABLE <- get("consDataTABLE",envir=nnTK) #list of imported metadata
  locs <- colnames(mixDataTABLE) #assume at least 1 evidence sample to get loci (decides order etc)
  
  resCompMAC <- get("resCompMAC",envir=nnTK)$MatchMatrix
  resCompLR1 <- get("resCompLR1",envir=nnTK)
  resCompLR2 <- get("resCompLR2",envir=nnTK)
  allMixList = get("allMixList",envir=nnTK) #show all mixtures
  resCompLR = get("resCompLR",envir=nnTK) #used to show match network
  resIBS <- get("resIBS",envir=nnTK)  #store IBS results (to be shown in report)
  resEvidConc <- get("resEvidConc",envir=nnTK)  #store concordant Evidence results (to be shown in report)
  resRMP <- get("resRMP",envir=nnTK)  #Provide random match probabilities
  resMatches <- get("resMatches",envir=nnTK)  #get match list of mixtures
  
  setupThresh = get("setupThresh",envir=nnTK) #get threshold setting 
  setupModel = get("setupModel",envir=nnTK) #get model setting 
  casedir <-  get("setupCase",envir=nnTK)$casepath #get case path
  
  selList = NULL #default is no selection
  if( get("setupAdvanced",envir=nnTK)$selProfiles=="TRUE" ) {
    print("-------USER SELECTION-------") #the user may select a subset of samples to import
    guienv = new.env( parent = emptyenv() ) #create new envornment object. Parent must be empty
    EVIDSall = rownames(mixDataTABLE)
    REFSall = rownames(refDataTABLE)
    assign("selected",list(evids=EVIDSall,refs=REFSall),envir=guienv)
    casesolver::profileSelectorGUI(env=guienv) #calling function to select data
    selList =  get("selected",envir=guienv) #get list of selection (evids,refs) (ONLY THESE SAMPLES WILL BE PRESENT IN RESULTS!)
    
    #GO THROUGH ALL DATA AND REMOVE PROFILES NOT SELECTED:
    EVIDS = selList$evids
    REFS = selList$refs
    EVIDSrm = setdiff(EVIDSall,EVIDS) #get evidence which are removed
    REFSrm = setdiff(REFSall,REFS) #get references which are removed
    
    #FILTERING ON LOADED DATA OBJECTS:
    mixLIST = mixLIST[EVIDS]
    mixDataTABLE = mixDataTABLE[rownames(mixDataTABLE)%in%EVIDS,,drop=FALSE] 
    refDataTABLE = refDataTABLE[rownames(refDataTABLE)%in%REFS,,drop=FALSE] 
    mixDataMATCHSTATUS = mixDataMATCHSTATUS[names(mixDataMATCHSTATUS)%in%EVIDS]
    mixDataMATCHSTATUS[mixDataMATCHSTATUS%in%REFSrm] = L$none #indicate no match if removed
    
    #Helpfunction for LR comparison lists (resCompLR, resCompLR1,resCompLR2)
    getCompList = function(compList) {
      retList = NULL #return list
      if(!is.null(compList)) { #Qualiative LR resulst
        indkeep = compList[,1]%in%EVIDS & compList[,2]%in%REFS #Keep only those with both considered EVIDS and REFS
        compList = compList[indkeep,,drop=FALSE] #keep only those with EVIDS
        if(ncol(compList)>0) retList = compList #copy
      }
      return(retList)
    } #end helpe function

        
    #Helpfunction for match lists (resMatches, allMixList)
    getMatchList = function(matchList) {
      retList = NULL #return list
      if(!is.null(matchList)) { #Qualiative LR resulst
        matchList = matchList[matchList[,1]%in%EVIDS,,drop=FALSE]
        if(nrow(matchList)>0) { #if at least one match in list
          refL = strsplit(matchList[,2],"/") #get list of references
          for(rr in 1:length(refL)) {
            newrefs = refL[[rr]] #get current matching refs
            newrefs = newrefs[newrefs%in%REFS] #keep only refs selected
            matchList[rr,2] = paste0(newrefs,collapse="/")
          }
          retList = matchList
        }
      }
      return(retList)
    } #end helpe function
    
    #FILTER FROM RESULTS ETC
    resMatches = getMatchList(resMatches)
    allMixList = getMatchList(allMixList)
    resCompLR1 = getCompList(resCompLR1)
    resCompLR2 = getCompList(resCompLR2)
    resCompLR = getCompList(resCompLR)
    
    #Other results
    if(!is.null(resCompMAC)) { #MATCH MATRIX
      colkeep = colnames(resCompMAC)%in%REFS
      rowkeep = rownames(resCompMAC)%in%EVIDS
      if( sum(colkeep)+sum(rowkeep) == 0) { #was rotated
        colkeep = colnames(resCompMAC)%in%EVIDS
        rowkeep = rownames(resCompMAC)%in%REFS
      }      
      resCompMAC = resCompMAC[rowkeep,colkeep,drop=FALSE] #UPDATE
    }
    if(!is.null(resEvidConc)) { #EVIDENCE CONCORDANCE
      resEvidConc = resEvidConc[resEvidConc[,2]%in%EVIDS & resEvidConc[,3]%in%EVIDS,,drop=FALSE] #keep only those with both REFS
    }
    
    if(!is.null(resIBS)) { #REF CONCORDANCE
      resIBS = resIBS[resIBS[,2]%in%REFS & resIBS[,3]%in%REFS,,drop=FALSE] #keep only those with both REFS
    }
    
    if(!is.null(resRMP)) { #RANDOM MATCH PROB
      resRMP$evidList = resRMP$evidList[resRMP$evidList[,2]%in%EVIDS,,drop=FALSE] #extract relevant
      resRMP$refList = resRMP$refList[resRMP$refList[,2]%in%REFS,,drop=FALSE] #extract relevant
    }
    #Not uesd: metaDataList, consDataTABLE
    
  }
  
  reportname <- gWidgets2::ginput( L$namingreportfile , text="report", title= L$userinput , icon="question")
  if(reportname=="") return() #report not created
  
  print("Creating report...")
  graphics.off() #close all plots before running..
  brdt1 <- 1 #inner bord type
  brdt2 <- NULL #outer bord type 
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  
  caseID  <- get("caseID",envir=nnTK) #get ID fromenvironment
  if(is.null(caseID)) caseID=0 #return() #make it possible to view report
  path <- paste0(casedir,.sep,caseID)
  if( !file.exists(path) ) path = getwd() #set path to working directory if path not found
  path2 <- paste0(path,.sep,"report")
  dir.create(path2, showWarnings = FALSE) #create folder if not existing
  
  #GENERATE HTML CODE:
  #cssfile <- "http://www.stat.ucl.ac.be/R2HTML/Pastel.css"
  cssfile <- system.file("samples", "R2HTML.css", package="R2HTML")
  htmlf <- R2HTML::HTMLInitFile(path2,filename=reportname,CSSFile=cssfile,Title=paste0("Case ",get("caseID",envir=nnTK)))
  
  #Header:    
  R2HTML::HTML( R2HTML::as.title(paste( L$reportforcase,get("caseID",envir=nnTK))),file = htmlf,HR=1)
  version =  packageVersion("casesolver") #obtain CS version
  
  if(checked[1]=="TRUE") { #show header?
    R2HTML::HTML( paste0("CaseSolver ",L$version," ",version," (euroformix_",packageVersion("euroformix"),").") ,file = htmlf)#,HR=0)
    R2HTML::HTML( R.version.string ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$user,": ",Sys.getenv("USERNAME")) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$created,": ",Sys.time()) ,file = htmlf)#,HR=0)
  }
  
  mixTab <- ssTab <- refTab <- NULL #empty
  if(!is.null(refDataTABLE)) {    #Add ref-table
    refTab <- addRownameTable(refDataTABLE,type=4,L$samplename)
  }
  if(!is.null(mixDataTABLE)) {    #Add evid-tables
    isMixture = rep(TRUE,nrow(mixDataTABLE)) #assume all is mixtures
    isMixture[mixDataMATCHSTATUS!="mixture"] = FALSE #these SS
    isMixture[ match(allMixList[allMixList[,3]=="1",1],rownames(mixDataTABLE)) ] = FALSE #ensure that it becomes SS if assigned as 1 contr.
    
    ssDataTABLE <-  cbind(mixDataMATCHSTATUS,mixDataTABLE)[!isMixture,,drop=FALSE]
    mixDataTABLE <-  mixDataTABLE[isMixture,,drop=FALSE]
    colnames(ssDataTABLE)[1] <- L$matchstatus #"MatchStatus"
    if(!checked[20]) ssDataTABLE = ssDataTABLE[,-1,drop=FALSE] #drop MatchStatus column
    ssTab <-  addRownameTable(ssDataTABLE,type=4,L$samplename)
    mixTab <- addRownameTable(mixDataTABLE,type=4,L$samplename)
  } 
  
  #Add peak heights in datatable:
  allTabLIST <- list() #store each sample in a table-list 
  if(!is.null(mixLIST)) {    #Add evid-tables
    sn <- names(mixLIST) #get sample names
    for(ss in sn) { #for each samples
      maxA <- max( sapply(mixLIST[[ss]],function(x) length(x$adata)) ) #maximum observed alleles (for a given evid)
      allTab <- as.numeric()
      for(loc in locs) { #for each locus
        newrow <- rep("",maxA) #used to insert alllele (height) info
        dat <- mixLIST[[ss]][[loc]]
        nA <- length(dat$adata) #number of alleles
        if(nA>0) { #require at least one allele
          ord <- order(dat$adata)
          suppressWarnings({ #handle warning when converting allele names
            tmp <- as.numeric(dat$adata) #converted alleles
          })
          if(!any(is.na(tmp))) ord <- order(tmp)     
          newrow[1:nA] <- paste0(dat$adata[ord]," (",dat$hdata[ord],")") #store alleles/heights
        }
        allTab <- rbind(allTab, c(loc,newrow) )
      }#end for each locus
      allTab <- t(allTab[,-1])
      colnames(allTab) <- locs
      allTabLIST[[ss]] <- allTab
    } #end for each samples
  } #end if add evid-table
  
  # R2HTML::HTML(as.title(paste0("Data")),file = htmlf,HR=1)
  
  ###############
  ###SHOW DATA###
  ###############
  if(checked[2]=="TRUE") { #show references?
    R2HTML::HTML(R2HTML::as.title( reportitems[2] ),file = htmlf,HR=2)
    if(!is.null(refTab)) {
      #if(nrow(refDataTABLE)>=nLarge) print("THIS MAY TAKE A WHILE!")
      R2HTML::HTML(as.data.frame(refTab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  if(checked[3]=="TRUE") { #show Single source profiles?
    R2HTML::HTML(R2HTML::as.title( reportitems[3] ),file = htmlf,HR=2)
    if(!is.null(ssTab)) {
      R2HTML::HTML(as.data.frame(ssTab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  if(checked[4]=="TRUE") { #show Mix profiles?
    R2HTML::HTML(R2HTML::as.title( reportitems[4] ),file = htmlf,HR=2)
    if(!is.null(mixTab)) {
      R2HTML::HTML(as.data.frame(mixTab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  if(checked[5]=="TRUE") { #Show consensus profiles?
    R2HTML::HTML(R2HTML::as.title( reportitems[5] ),file = htmlf,HR=2)
    if(!is.null(consDataTABLE)) {
      sn <- unique(consDataTABLE[,1]) #get sample names
      consDataOUT <- matrix(ncol=length(locs)+1,nrow=length(sn))
      consDataOUT[,1] <- sn
      for(ss in sn) { #for each samples
        for(loc in locs) { #for each locus
          consDataOUT[which(sn==ss),which(locs==loc)+1] <- consDataTABLE[consDataTABLE[,1]==ss &consDataTABLE[,2]==loc,3]
        }
      }
      colnames(consDataOUT) <- c(L$samplename ,locs)
      R2HTML::HTML(as.data.frame(consDataOUT), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  ###############
  ###SHOW W/PH###
  ###############
  if(checked[6]=="TRUE") { #should SS with peak heights be included?
    R2HTML::HTML(R2HTML::as.title( reportitems[6] ),file = htmlf,HR=2)
    if(!is.null(mixLIST) && !is.null(ssTab)) {
      for(ss in ssTab[,1] ) { #for each single sources
        R2HTML::HTML(R2HTML::as.title( paste0(ss) ),file = htmlf,HR=3)
        if(nrow(allTabLIST[[ss]])==0) { #check if empty
          R2HTML::HTML(paste0( L$none ),file = htmlf)
        } else {
          R2HTML::HTML(as.data.frame(addRownameTable(allTabLIST[[ss]],type=0,L$samplename)), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
        }  
      }
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  if(checked[7]=="TRUE") { #should Mixture with peak heights be included?
    R2HTML::HTML(R2HTML::as.title( reportitems[7] ),file = htmlf,HR=2)
    if(!is.null(mixLIST) && !is.null(mixTab)) {
      for(ss in mixTab[,1] ) { #for each mixtures
        R2HTML::HTML(R2HTML::as.title( paste0(ss) ),file = htmlf,HR=3)
        if(nrow(allTabLIST[[ss]])==0) { #check if empty
          R2HTML::HTML(paste0( L$none ),file = htmlf)
        } else {
          R2HTML::HTML(as.data.frame(addRownameTable(allTabLIST[[ss]],type=0,L$samplename)), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
        }  
      }
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  
  if(checked[8]=="TRUE") { #Show metadata?
    R2HTML::HTML(R2HTML::as.title( reportitems[8] ),file = htmlf,HR=2)
    if(length(metaDataList)>0) { #is empty list?
      for(elem in names(metaDataList)) {
        if(length(metaDataList[[elem]])==0) next #skip if no info
        R2HTML::HTML(R2HTML::as.title( paste0(elem) ),file = htmlf,HR=3)
        R2HTML::HTML(metaDataList[[elem]], file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)       
      }
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  #Helpfunction for insert result table
  easyAdd = function(X,type=2,inclRN=TRUE) {
    R2HTML::HTML(as.data.frame(addRownameTable(X,type=type,L$samplename)), file= htmlf,row.names=inclRN,align="left",innerBorder = brdt1,Border=brdt2)
  }
  
  #################
  ###COMPARISONS###
  #################
  if(any(checked[9:12]=="TRUE")) { #Any comparisons to show?
    R2HTML::HTML(R2HTML::as.title( L$comparisons ),file = htmlf,HR=1)
  }
  
  #Provide match matrix:
  if(checked[9]=="TRUE") {  #Show match matrix?
    R2HTML::HTML(R2HTML::as.title( reportitems[9] ),file = htmlf,HR=2)
    
    #INSERT MATRIX
    if( !is.null(resCompMAC) ) { #if completed
      if(nrow(resCompMAC)>0) {
        ind <- as.numeric(resCompMAC)<setupThresh$MAC #get smaller indices
        resCompMAC[ind] <- "" #show only greater than threshold
        easyAdd(resCompMAC,type=4)
      }
    } else {
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  #Provide match list 1 (qual based):
  if(checked[10]=="TRUE") {
    R2HTML::HTML(R2HTML::as.title( reportitems[10] ),file = htmlf,HR=2)
    
    if(!is.null(resCompLR1) ) { #if completed
      resCompLR1 <- resCompLR1[as.numeric(resCompLR1[,4])>=log10(setupThresh$LRthresh1),,drop=FALSE]
      if(nrow(resCompLR1)>0) {
        easyAdd(resCompLR1,type=0)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else {
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  #Provide match list 2 (quan based):
  if(checked[11]=="TRUE") {
    R2HTML::HTML(R2HTML::as.title( reportitems[11] ),file = htmlf,HR=2)
      
    if(!is.null(resCompLR2) ) { #if completed
      resCompLR2 <- resCompLR2[as.numeric(resCompLR2[,4])>=log10(setupThresh$LRthresh2),,drop=FALSE]
      if(nrow(resCompLR2)>0) {
        easyAdd(resCompLR2,type=0)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else {
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  #Provide match list Final:
  if(checked[12]=="TRUE") { #Show match list?
    R2HTML::HTML(R2HTML::as.title( reportitems[12] ),file = htmlf,HR=2)
    if(!is.null(allMixList) && nrow(allMixList)>0) {
      easyAdd(allMixList,type=0)
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }
  }
  
  if(checked[13]=="TRUE") { #Show match network?
    R2HTML::HTML(R2HTML::as.title( reportitems[13] ),file = htmlf,HR=2)
    okplot <- !is.null(resCompLR)
    if(okplot) {  #only add if OK
      netf <- file.path(path2,"matchnetwork.png") #file of picture 
      png(netf,width = 2000, height = 2000,res=200)
      showMatchNetwork(nnTK,"all",createInteractive=FALSE,selList)
      dev.off()
      R2HTML::HTMLInsertGraph(netf ,file=htmlf, Align = "left", WidthHTML = 1000)
    } else { 
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  ##################
  ###IBS/evidConc###
  ##################
  
  if(checked[16]=="TRUE") { #Show IBS table?
    R2HTML::HTML(R2HTML::as.title( reportitems[16] ),file = htmlf,HR=2)
    if(!is.null(resIBS)) { 
      if(nrow(resIBS)>0) {
        R2HTML::HTML(as.data.frame(resIBS), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else { 
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  if(!is.na(checked[21]) && checked[21]=="TRUE") { #Show convEvid table? (Notice backward compatibility)
    R2HTML::HTML(R2HTML::as.title( reportitems[21] ),file = htmlf,HR=2)
    if(!is.null(resEvidConc)) {
      if(nrow(resEvidConc)>0) {
        R2HTML::HTML(as.data.frame(resEvidConc), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else { 
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  ##############
  ###RMP/RMNE###
  ##############
  
  #store random match prob results (to be shown in report)
  if(any(checked[14:15])) { #Any comparisons to show?
    R2HTML::HTML(R2HTML::as.title( L$randommatchprob ),file = htmlf,HR=1)
  }
  if(checked[14]=="TRUE") { #Show RMNE?
    R2HTML::HTML(R2HTML::as.title( reportitems[14] ),file = htmlf,HR=2)
    if(!is.null(resRMP)) { #if Inclusion probabilities are calculated
      if(nrow(resRMP$evid)>0) {
        R2HTML::HTML(as.data.frame(resRMP$evid), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else { 
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  if(checked[15]=="TRUE") { #Show RMP probs?
    R2HTML::HTML(R2HTML::as.title( reportitems[15] ),file = htmlf,HR=2)
    if(!is.null(resRMP)) { #if random match probabilities are calculated
      if(nrow(resRMP$ref)>0) {
        R2HTML::HTML(as.data.frame(resRMP$ref), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
      } else {
        R2HTML::HTML(paste0( L$nocandidates ),file = htmlf)#,HR=10)
      }
    } else { 
      R2HTML::HTML(paste0( L$notcompleted ),file = htmlf)#,HR=10)
    }
  }
  
  ##############
  ###Settings###
  ##############
  
  canPrintEPG = function(kit) { #Helpfunction to check if can print EPG (kit specified)
    if(is.null(kit) || kit==L$none) return(FALSE)
    kitinfo <- euroformix::getKit(kit)
    return( length(kitinfo)>1 )  #possible to print EPG?
  }
  
  kit0 <- get("setupKit",envir=nnTK)$kitname
  printEPG <- canPrintEPG(kit0)  && !is.null(mixLIST) #possible to print EPG?
  if(checked[17]=="TRUE") { #should settings be shown?
    R2HTML::HTML(R2HTML::as.title(paste0( reportitems[17] )),file = htmlf,HR=1)
    R2HTML::HTML(R2HTML::as.title(paste0( L$threshs," :")),file = htmlf,HR=2)
    R2HTML::HTML( paste0( L$macthreshold ,colonsymbol,setupThresh$MACthresh) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$qualLRthreshold ,colonsymbol,setupThresh$LRthresh1) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$quanLRthreshold ,colonsymbol,setupThresh$LRthresh2) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$minLocSSmatch ,colonsymbol,setupThresh$minLociSS) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$minIBSrelative ,colonsymbol,setupThresh$minIBS) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$probRatioToNext ,colonsymbol,setupThresh$ratio) ,file = htmlf)#,HR=0)
    R2HTML::HTML( paste0( L$probSingleAllele ,colonsymbol,setupThresh$probA) ,file = htmlf)#,HR=0)
    
    R2HTML::HTML(R2HTML::as.title(paste0( L$model, colonsymbol)),file = htmlf,HR=2)
    R2HTML::HTML( paste0( L$popfreq ,colonsymbol,get("setupPop",envir=nnTK)$popfile) ,file = htmlf)#,HR=0)
    
    modeltype=setupModel$modeltype
    if(modeltype%in%c(1,3)) {
      R2HTML::HTML(R2HTML::as.title(paste0( L$qualmodel ,colonsymbol)),file = htmlf,HR=3)
      R2HTML::HTML( paste0( L$dropinprob ,colonsymbol,setupModel$dropinC) ,file = htmlf)#,HR=0)
    } 
    if(modeltype%in%c(2,3)) {
      kit1 = kit0
      if(is.null(kit1) || kit1=="") kit1 = L$none #indicate if not selected kit
      R2HTML::HTML(R2HTML::as.title(paste0( L$quanmodel ,colonsymbol)),file = htmlf,HR=3)
      R2HTML::HTML( paste0( L$kit ,colonsymbol,kit1) ,file = htmlf)#,HR=0)
      R2HTML::HTML( paste0( L$analyticalthreshold, colonsymbol,setupModel$threshT) ,file = htmlf)#,HR=0)
      R2HTML::HTML( paste0( L$degradationmodel, colonsymbol, radiotxt[setupModel$degrad] ) ,file = htmlf)#,HR=0)
      R2HTML::HTML( paste0( L$stuttermodel,colonsymbol, radiotxt[setupModel$stutt]) ,file = htmlf)#,HR=0)
      R2HTML::HTML( paste0( L$dropinprob ,colonsymbol,setupModel$dropinC) ,file = htmlf)#,HR=0)
      R2HTML::HTML( paste0( L$dropinpeakheightlambda, colonsymbol,setupModel$dropinL) ,file = htmlf)#,HR=0)
    } 
  }
  
  ################
  ###ATTACHMENT###
  ################
  if(any(checked[18:19])) { #Any comparisons to show?
    R2HTML::HTML(R2HTML::as.title( L$attachments ),file = htmlf,HR=1)
  }
  whsize <- c(1920*2,1080*2,120*2) #number of pixels and resolution
  sampleType = getSampleType2(mixLIST,kit0) #get sample type
  
  #refL = NULL # list of reference data 
  #if((printEPG || sampleType=="LUS") && any(checked[setdiff(18:21,20)]=="TRUE") ) { #if showing plot
    
  #Plot EPG figures for single sources
  if((printEPG || sampleType=="LUS") && checked[18]=="TRUE") { 
    R2HTML::HTML(R2HTML::as.title(paste0( reportitems[18] )),file = htmlf,HR=2)
    if(nrow(ssDataTABLE)>0) {
      unREF = unique(ssDataTABLE[,1]) #get unique refs
      refL = getRefL(unREF,forceDi=FALSE) # get relevant references
      
      for(i in 1:nrow(ssDataTABLE)) { #for each single source profiles
        evid <- rownames(ssDataTABLE)[i]
        ref <- ssDataTABLE[i,1]
        
        if(evid== L$empty) {
          R2HTML::HTML( L$none ,file = htmlf)
        } else {
          condref <- refL[ref] #extract reference
          if(length(condref)==0) condref = NULL
          
          tryCatch({ suppressWarnings({
            epgf <- file.path(path2,paste0("epg_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
            png(epgf ,width =whsize[1], height = whsize[2],res=whsize[3])
            if(sampleType=="EPG") euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE,threshT=setupModel$threshT)
            if(sampleType=="LUS") euroformix::plotLUS(mixLIST[evid],sn=evid,condref,threshT=setupModel$threshT,LUSsymbol=LUSsymbol)
            dev.off()
            R2HTML::HTMLInsertGraph(epgf,file=htmlf, Align = "left", WidthHTML =  whsize[2]*0.8)
          }) })
        } #end if not empty
        R2HTML::HTML(R2HTML::as.title(paste0("#",i," - ",evid)),file = htmlf,HR=3)
      } #end for each samples
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }#end if
  } #end if plot EPG
  
  #Plot EPG figures for mixtures?
  if((printEPG || sampleType=="LUS") && checked[19]=="TRUE") { 
    R2HTML::HTML(R2HTML::as.title(paste0( reportitems[19] )),file = htmlf,HR=2)
    if(nrow(mixDataTABLE)>0) {
      if(!is.null(resMatches) && nrow(resMatches)>0) { #require match table
        unREF = unique(unlist(strsplit(resMatches[,2],"/"))) #get unique refs
        refL = getRefL(unREF,forceDi=FALSE) # get relevant references
      }
      
      for(i in 1:nrow(mixDataTABLE)) { #for each single source profiles
        evid <- rownames(mixDataTABLE)[i]
        condref = NULL
        if(!is.null(resMatches) && nrow(resMatches)>0) { #require match table
          ind <- resMatches[,1]%in%evid
          if(any(ind)) {
            refs <- unlist(strsplit(resMatches[ind,2],"/")) #get refs
            condref <- refL[refs]
            if(length(condref)==0) condref = NULL
          }
        }
        tryCatch({ suppressWarnings({
          epgf <- file.path(path2,paste0("epg_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
          png(epgf ,width =whsize[1], height = whsize[2],res=whsize[3])
          if(sampleType=="EPG") euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE,threshT=setupModel$threshT)
          if(sampleType=="LUS") euroformix::plotLUS(mixLIST[evid],sn=evid,condref,threshT=setupModel$threshT,LUSsymbol=LUSsymbol)
          dev.off()
          R2HTML::HTMLInsertGraph(epgf,file=htmlf, Align = "left", WidthHTML = whsize[2]*0.8)
          R2HTML::HTML(R2HTML::as.title(paste0("#",i," - ",evid)),file = htmlf,HR=3)
        }) })
      } #end for each samples
    } else {
      R2HTML::HTML(paste0( L$none ),file = htmlf)#,HR=10)
    }#end if
  } #end if plot EPG
  
  print("REPORT SAVED:")
  print(htmlf)
  browseURL(htmlf) #look directly in browser after creating report
} #end create report