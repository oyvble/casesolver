#' @title createReport
#' @description A function for creating report based on CaseSolver object
#' @details The function generated reports based on case analysis. Can be HTML (R2HTML), DOC (rtf), DOCX (officer/flextable)
#' @param nnTK an environment object from stored CaseSolver object
#' @export

#library(casesolver);load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/TutorialDataCaseSolver/CS2_FUSION.Rdata")
#library(casesolver);load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/TutorialDataCaseSolver/cs_fusion.Rdata")
#library(casesolver);load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/TutorialDataCaseSolver/cs_error.Rdata")
#library(casesolver);load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/Issues/19-01209.Rdata")
#source("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolver/R/createReport.R");

createReport = function(nnTK) { #this function loads data and put them in a HTML script

  #################################
  #Obtain options from environment:
  exp = get("setupReportExpTyp",envir=nnTK) #get export options (report and preview)
  formatNames = c("HTML","DOCX","DOC") #names(exp) #
  ext = c("html","docx","doc")  #extension of reportnames
  createHTML = exp$HTML #whether to create report based on HTML
  createDOCX = exp$DOCX #Whether to create report based on DOCX (advanced layout). 
  createDOC = exp$DOC #Whether to create report based on DOC (simple layout)
  
  reportTemplate = "reportTemplate.docx" #file name of report template
  
  #Report type settings  
  opt = get("setupReportExpOpt",envir=nnTK) #get export options (report and preview)
  WHsize = opt$PNG #obtain PNG settings (Number of pixels (width,height) and resolution of png figs)
  if(is.null(WHsize)) WHsize = c( 1920,1080,120 )
  HTMLoptions = list(brdt1=1, brdt2=NULL) #inner/outer bord type 
  
  vals = opt$DOC  #obtain DOC settings (width,height, margin)
  if(is.null(vals)) vals = c(11,8.5,0.1,6,9) #Options for RT format (put in landscape)
  DOCoptions= list(width=vals[1],height=vals[2],margin=vals[3],tableSize=vals[4], fontSize=vals[5])  
   
  vals = opt$DOCX  #obtain DOC settings (width,height, margin)
  if(is.null(vals)) vals = c(6,11) #Options for RT format (put in landscape)
  DOCXoptions= list(tableSize=vals[1],fontSize=vals[2])   

  #Report options (only boolean)  
  showItem = get("setupReportOpt",envir=nnTK) #get export options (report and preview)
  #names(showItem)  "MatchStatus"  "MCMCsettings" "mleLR"        "bayesLR"  "consLR"       "Mx"           "validFailed" 
  
  #obtain  
  defaultReportName = "report"
  LUSsymbol <- "_"
  colonsymbol = ":" #use variable for colon
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  # L = getLanguage( "English")
  if(is.null(L)) {
    print("No language found.")
    return()
  }
  
  formatUse = setNames(c(createHTML[1],createDOCX[1],createDOC[1]),formatNames)
  if(!any(formatUse)) {
    print("No report formats were selected. Returning from function.")
    return()
  } else {
    print("Following report formats will be created:")
    print(paste0(formatNames[formatUse],collapse=","))
  }
  
  if(formatUse[1] && !require(R2HTML,quietly = T)) print("R2HTML must be installed in order to create HTML report.")
  if(formatUse[3] && !require(rtf,quietly = T)) print("rtf must be installed in order to create DOC report.")
  if(formatUse[2]) {
    txt = " must be installed in order to create DOCX report."
    if(!require(officer,quietly = T)) print(paste0("officer",txt))
    if(!require(flextable,quietly = T))print(paste0("flextable",txt))
  }
  
  #helpfunction when extracting reference data
  getRefL = function(refs,forceDi=FALSE) { #return list with same order as for refs
    refL <- casesolver::tabToListRef(tab=refDataTABLE[ match(refs,rownames(refDataTABLE)),,drop=FALSE],setEmpty=FALSE) #FORCING DUP alleles
    return(refL)
  }
  
  #prepare layout (own setting)
  optLay <- get("setupReportLay",envir=nnTK)  #get report layout
  priority = as.integer(optLay$priority)
  reportitems = optLay$reportitems
  names(priority) = reportitems
  radiotxt = c( L$on , L$off )
  
  #Provide data to use in report: May be modified if user deselects samples!
  mixLIST <- get("mixDataLIST",envir=nnTK) #get mix list
  mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #get data table from nnTK-environment
  refDataTABLE <- get("refDataTABLE",envir=nnTK) #get data table from nnTK-environment
  mixDataMATCHSTATUS <- get("mixDataMATCHSTATUS",envir=nnTK) #assign to nnTK-environment
  metaDataList <- get("metaDataLIST",envir=nnTK) #list of imported metadata
  consDataTABLE <- get("consDataTABLE",envir=nnTK) #list of imported metadata
  
  resCompMAC <- get("resCompMAC",envir=nnTK)$MatchMatrix
  resCompLR1 <- get("resCompLR1",envir=nnTK)
  resCompLR2 <- get("resCompLR2",envir=nnTK)
  allMixList = get("allMixList",envir=nnTK) #show all mixtures
  resCompLR = get("resCompLR",envir=nnTK) #used to show match network
  resIBS <- get("resIBS",envir=nnTK)  #store IBS results (to be shown in report)
  resEvidConc <- get("resEvidConc",envir=nnTK)  #store concordant Evidence results (to be shown in report)
  resRMP <- get("resRMP",envir=nnTK)  #Provide random match probabilities
  resMatches <- get("resMatches",envir=nnTK)  #get match list of mixtures
  DCdataTABLE <- get("DClistReport",envir=nnTK) #obtain extracted candidates from deconvolution
  resWOE <- get("resWOEeval",envir=nnTK) #obtain extracted candidates from deconvolution
  
  setupPop = get("setupPop",envir=nnTK) #population settings
  setupRare = get("setupRare",envir=nnTK) #rare allele settings
  setupThresh = get("setupThresh",envir=nnTK) #get threshold setting 
  setupModel = get("setupModel",envir=nnTK) #get model setting 
  setupMarkers = get("setupMarkers",envir=nnTK) #get model setting (per marker) 
  setupMCMC = get("setupMCMC",envir=nnTK) #get MCMC setting (used obtaining bayes/cons LR)
  casedir =  get("setupCase",envir=nnTK)$casepath #get case path
  
  #Obtain locus names from data tables
  locs = NULL
  if(!is.null(mixDataTABLE)) {
    locs <- colnames(mixDataTABLE) #obtain locus names from evidence Data
  } else if(!is.null(refDataTABLE)) {
    locs <- colnames(refDataTABLE) #obtain locus names from reference Data
  } 
  
  #Obtain locus names to use in data tables (from settings)
  locNamesReport <- get("setupReportLocNames",envir=nnTK) #obtain (possibly customized) marker names to be used in report
  if(is.null(locNamesReport)) { #if locus names not stored
    locNamesReport = setNames(locs,locs)
  } 
  locNamesTables = setNames(locNamesReport[locs], locs) #insert modified marker names and ordinary marker names
  insNA = is.na(locNamesTables)
  locNamesTables[insNA] = locs[insNA] #insert ordinary marker names for those missing

  #INSERT MODIFIED MARKER NAMES AS NEW COLNAMES TO ESSENTIAL TABLES
  if(!is.null(mixDataTABLE)) colnames(mixDataTABLE) <-  locNamesTables #update marker names
  if(!is.null(refDataTABLE)) colnames(refDataTABLE) <-  locNamesTables #update marker names

    
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
    DCdataTABLE = DCdataTABLE[rownames(DCdataTABLE)%in%REFS,,drop=FALSE] 

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
  } #end if specific selected profiles
  
  reportname <- gWidgets2::ginput( L$namingreportfile , text=defaultReportName, title= L$userinput , icon="question")
  if(reportname=="") return() #report not created
  
  #prepare data
  print("Preparing data for report...")
  sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types (different tables)
  #Note: using orderTableSort(va1,var2,sort) to obtain order of table
  
  #Prepare evidence data
  mixTab <- ssTab <- NULL #empty
  if(!is.null(mixDataTABLE)) {    #Add evid-tables
    isMixture = rep(TRUE,nrow(mixDataTABLE)) #assume all is mixtures
    isMixture[mixDataMATCHSTATUS!="mixture"] = FALSE #these SS
    isMixture[ match(allMixList[allMixList[,3]=="1",1],rownames(mixDataTABLE)) ] = FALSE #ensure that it becomes SS if assigned as 1 contr.
    
    if(sum(!isMixture)>0) { #if at least one single source
      ssDataTABLE <-  cbind(mixDataMATCHSTATUS,mixDataTABLE)[!isMixture,,drop=FALSE]
      colnames(ssDataTABLE)[1] <- L$matchstatus #"MatchStatus"
      matchStatus = ssDataTABLE[,1] #obtain match status for single sources
      if(!showItem[1]) {
        ssDataTABLE = ssDataTABLE[,-1,drop=FALSE] #drop MatchStatus column if item not to be shown
        matchStatus = NULL
      } 
      ordSS = casesolver::orderTableSort(rownames(ssDataTABLE),matchStatus,sortTypes[1]) #obtain selected sorted order for evids
      ssTab <-  addRownameTable(ssDataTABLE[ordSS,,drop=FALSE],type=4,L$samplename)
    }
    
    if(sum(isMixture)>0) { #if at least one mixture
      mixDataTABLE <-  mixDataTABLE[isMixture,,drop=FALSE]
      ordMIX = casesolver::orderTableSort(rownames(mixDataTABLE),sort=sortTypes[1]) #obtain selected sorted order for evids
      mixTab <- addRownameTable(mixDataTABLE[ordMIX,,drop=FALSE],type=4,L$samplename)
    }
  } 
  
  #Add peak heights in datatable:
  allTabLIST <- list() #store each sample in a table-list 
  if(!is.null(mixLIST)) {    #Add evid-tables
    sn <- names(mixLIST) #get sample names
    for(ss in sn) { #for each samples
      maxA <- max( sapply(mixLIST[[ss]],function(x) length(x$adata)) ) #maximum observed alleles (for a given evid)
      allTab <- as.numeric()
      for(loc in locs) { #for each (ordinary) locus
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
      colnames(allTab) <- locNamesTables #insert locus names for report
      allTabLIST[[ss]] <- allTab
    } #end for each samples
  } #end if add evid-table
  
  #Prepare DC results:
  DCdataList = list()
  ptrn = "-" #pattern used for separating sample name and component
  if(!is.null(DCdataTABLE)) {
    DCsn = rownames(DCdataTABLE) #obtain sample names for DC results
    DCcn= colnames(DCdataTABLE) #obtain column names for DC results
    locUse = locs[locs%in%DCcn] #use only loci present in column names
    
    indExtended = duplicated(DCcn)
    indUse1 = match(locUse,DCcn[!indExtended]) 
    indUse2 = match(locUse,DCcn[indExtended])
    
    for(i in 1:nrow(DCdataTABLE)) {
      
      tmp = strsplit(DCdataTABLE[i,1],ptrn)[[1]] #obtain sample name + component
      sn = paste0(tmp[1:(length(tmp)-1)],collapse=ptrn) #obtain sample name
      #comp = tmp[length(tmp)] #get component (not used)
      
      #Obtain conditional (column 2) and add together with name
      cond0 = DCdataTABLE[i,2]  #conditional reference
      cond = paste0(DCcn[2],"=",cond0) #conditional references
      noc = paste0(DCcn[3],"=",DCdataTABLE[i,3]) #noc assumption
      mx = paste0(DCcn[4],"=",DCdataTABLE[i,4]) #mix prop
      condtxt  = c(sn,noc,mx)
      if(nchar(cond0)>0)  condtxt = c(condtxt,cond)
      condtxt = paste0(condtxt,collapse=" - ")
      sn = paste0(DCsn[i]," (",condtxt,")") #obtain final sample name
      
      newTab = DCdataTABLE[i,indUse1]
      newTab = rbind(newTab,DCdataTABLE[i,indExtended][indUse2]) 
      colnames(newTab) = locNamesTables[locUse] #insert locus names (names for report)
      rownames(newTab) = c("*","**")
      DCdataList[[sn]] = newTab #sample name used as item name 
    }
  }
  DCtext = c( paste0("*",L$DCreporttext1) , paste0("**",L$DCreporttext2))#Create text to add after
  
  #Prepare WoE results (global results, statements, per-marker results, parameter results:
  WOEdataTABLE <- WOEmarkerTABLE <- NULL
  WOEstateList <- WOEparamList <- list()
  if(!is.null(resWOE)) {
    s0 = 2 #number of signif for LRper marker
    WOEdataTABLE = resWOE$resTable #
    nHyps = nrow(WOEdataTABLE) #number of hypotheses
    
    #Subselect columns (chosen in report options)
    colRm = (5:9)[!showItem[3:7]] #columns to remove (based on user selection)
    if(length(colRm)>0) WOEdataTABLE <- WOEdataTABLE[,-colRm,drop=FALSE] #remove columns if any to remove
  
    #Create LR per marker table
    WOEmarkerTABLE = t(sapply(resWOE[1:nHyps],function(x) signif(x$mleLRi,s0)))
    
    #insert rownames and update 
    rownames(WOEmarkerTABLE) <- rownames(WOEdataTABLE) <- 1:nHyps 
    colnames(WOEmarkerTABLE) <- locNamesTables[colnames(WOEmarkerTABLE)] #update with new marker names
    
    #Create statement and param lists
    for(id in 1:nHyps) { #for each hypothesis #id
      hypname = paste0(L$hypothesisset," #",id)
      WOEstateList[[hypname]] = resWOE[[id]]$statement #insert statement
      WOEparamList[[hypname]] = resWOE[[id]]$paramTable #insert parameter table
    }
  } 
  WOEtext = paste0("*",L$WOEreporttext) #Create text to add after
    
  #Prepare reference data (split up in known vs extracted)
  #Separate References and Estimated References (unknown and DC estimated)
  refTabKnown <- refTabExtracted <- NULL #empty by default
  if(!is.null(refDataTABLE)) { 
    refNames = rownames(refDataTABLE) #obtain reference names
    isUnknown = substr(refNames,1,nchar(L$unknown))==L$unknown #indicate which refs that are "unknown"
    isDCed = refNames%in%rownames(DCdataTABLE) #indicate which refs that are "DCed"
    isExtracted = isUnknown | isDCed #indicate which references are extracted

    if(any(!isExtracted)) {
      refTab = refDataTABLE[!isExtracted,,drop=F] #obtain known references
      ord = casesolver::orderTableSort(rownames(refTab),sort=sortTypes[2]) #obtain selected sorted order for refs
      refTabKnown <- addRownameTable(refTab[ord,,drop=FALSE],type=4,L$samplename)
    }
    if(any(isExtracted)) {
      refTab = refDataTABLE[isExtracted,,drop=F] #obtain known references
      ord = casesolver::orderTableSort(rownames(refTab),sort=sortTypes[2]) #obtain selected sorted order for refs
      refTabExtracted <- addRownameTable(refTab[ord,,drop=FALSE],type=4,L$samplename)
    }
  }

  ################################################### 
  #Helpfunctions to ease the insertion result tables#
  ################################################### 
  
  #Helpfunction to add text (including n line shifts)
  insText = function(txt,bold=FALSE,italic=FALSE,n=1) {
    txt2 = paste0(txt,collapse="\n")
    if(formatUse[1]) {
      sapply(txt, function(x) R2HTML::HTML( x, file = htmlf) )
    } 
    if(formatUse[3]) {
      rtf::addText(rtf, txt2 ,bold=bold,italic=italic)
      for(i in 1:n) rtf::addNewLine(rtf)  #print new linest
    }
    if(formatUse[2]) {
      sapply(txt, function(x) officer::body_add_par( docx,x ) ) 
    }
  }
  #helpfunction to add title
  insTitle = function(txt, lvl=1) {
    if(formatUse[1]) {
      R2HTML::HTML(R2HTML::as.title( txt ),file = htmlf,HR=lvl)
    }
    if(formatUse[3]) {
      rtf::addNewLine(rtf) #add a line before insTitle
      rtf::addText(rtf, txt ,bold=TRUE)
      rtf::addNewLine(rtf) #add a line before insTitle
    }
    if(formatUse[2]) { #look at styles_info(docx)
      if(lvl==0) {
        headtype = "Title"
      } else {
        headtype = paste0("heading ",lvl)
      }
      x=officer::body_add_par(docx,txt,style=headtype);
    }
  }

  #HELPFUNCTION TO INSERT WHOLE TABLE
  insTable = function(tab,title=NULL,type=NULL,lvl=2,addRowname=TRUE) {
    if(!is.null(title)) insTitle(title,lvl) #add title name if indicated
    
    if(!is.null(tab) && nrow(tab)>0) {
      #PRE:
      if(!is.null(type)) tab = casesolver::addRownameTable(tab,type=type,L$samplename) #update table by adding rowname to table

      if(formatUse[1]) {
        R2HTML::HTML( as.data.frame(tab), file= htmlf,row.names=addRowname,align="left",innerBorder = HTMLoptions$brdt1,Border=HTMLoptions$brdt2)
      }
      
      tab2 <- tab #copy table
      if(addRowname) { #if add rowname
        tab2 = cbind( rownames(tab),tab) #add rownames table
        colnames(tab2) = c("#",colnames(tab)) 
      }
      tab2 = as.data.frame(tab2)
      if(formatUse[3]) {
        rtf::addTable(rtf, tab2, font.size = DOCoptions$tableSize, col.justify="C",header.col.justify="C")#row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
        rtf::addNewLine(rtf)
      }
      if(formatUse[2]) {
        #style = DOCXoptions$tabStyle
        #x=officer::body_add_table(docx,tab, style = style, alignment="c");
        ft <- flextable::flextable(tab2)
        #ft <- flextable::qflextable(tab)
        ft <- flextable::fontsize(ft, part = "all", size = DOCXoptions$tableSize)
        ft <- flextable::set_table_properties(ft, layout = "autofit")
        flextable::body_add_flextable(docx,value = ft )
      }
    } else {
      insText( L$none ,italic=T,n=2)
    }
  } 

  #HELPFUNCTION TO INSERT TABLES FOR EACH SAMPLE (list)
  insList = function(List, title,selected=NULL, type=0) {
    insTitle(title,2)
    if(is.null(List) || (!is.null(selected) && selected=="NONE") || length(List)==0 ) { #if no elements 
      insText( L$none ,italic=T,n=2)
    } else {
      if(is.null(selected)) selected = names(List) #select all
          
      for(ss in selected) { #for each element
        listElem = List[[ss]] #obtain list element
        insTitle(ss,3) #set sub-title with element name
        
        if(is.null(dim(listElem))) { #if not a table
          insText( listElem ) #insert text
        } else { #if a table
          if( nrow(listElem)==0 ) { #check if empty
            insText( L$none ,italic=T,n=2)
          } else {
            insTable(listElem,type=type) #otherwise print table
          }
        }
      }
    }      
  }
  
  #helpfunction to insert figs  
  insIMG = function(fn, quadratic=FALSE) {

    if(formatUse[1]) {
      HTMLwidth = WHsize[1]
      if(quadratic) HTMLwidth = WHsize[2]
      R2HTML::HTMLInsertGraph(fn ,file=htmlf, Align = "left", WidthHTML = HTMLwidth)
    }
    
    if(formatUse[3]) {
      s1 = DOCoptions$width #WHsize
      s2 = DOCoptions$height #WHsize
      scale = (s1/s2)/(WHsize[1]/WHsize[2])
      s2 = scale*s2 #scale height to make correct ratio
      if(quadratic) s1=s2
      
      rtf::addPng(rtf, fn, s1, s2) #insert image
    }
    if(formatUse[2]) {
      pagedim = docx$sect_dim$page/1500 #obtain page dim (manual scaling)
      s1 = pagedim[1] #WHsize
      s2 = pagedim[2] #WHsize
      scale = (s1/s2)/(WHsize[1]/WHsize[2])
      s2 = scale*s2 #scale height to make correct ratio
      if(quadratic) s1=s2
       x=officer::body_add_img(docx,fn,width = s1, height=s2);
    }
  }
  
  #Extract other info about printing EPG (kitname etc)
  canPrintEPG = function(kit) { #Helpfunction to check if can print EPG (kit specified)
    if(is.null(kit) || kit==L$none) return(FALSE)
    kitinfo <- euroformix::getKit(kit)
    return( length(kitinfo)>1 )  #possible to print EPG?
  }
  
  #Check whether EPGs can be generated (also depends on type of data)
  printEPG <- casesolver::canPrintEPG(nnTK) && !is.null(mixLIST) #possible to print EPG? Only if valid kit specified
  kit0 <- get("setupKit",envir=nnTK)$kitname
  sampleType = casesolver::getSampleType2(mixLIST,kit0) #get sample type
  
###################
###GENERATE REPORT#
###################
  
  graphics.off() #close all plots before running..
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  caseID  <- get("caseID",envir=nnTK) #get ID fromenvironment
  if(is.null(caseID)) caseID=0 #return() #make it possible to view report
  path <- paste0(casedir,.sep,caseID)
  if( !file.exists(path) ) path = getwd() #set path to working directory if path not found
  path2 <- paste0(path,.sep,"report")
  dir.create(path2, showWarnings = FALSE) #create folder if not existing
  reportfn <-  paste0(path2,.sep,reportname,".",ext)#obtain full path of report names

  #obtain CS version:
  version =  packageVersion("casesolver") 
  
  #Create report objects:
  doctxt = paste0("Case ",get("caseID",envir=nnTK)) #Text outside main body
  if(formatUse[1]) {
    #cssfile <- "http://www.stat.ucl.ac.be/R2HTML/Pastel.css"
    cssfile <- system.file("samples", "R2HTML.css", package="R2HTML")
    htmlf <- R2HTML::HTMLInitFile(path2,filename=reportname,CSSFile=cssfile,Title = doctxt)
  }
  if(formatUse[2]) {
    template <- system.file(reportTemplate,  package="casesolver") #obtain inbuilt scheme from casesolver
    docx <- officer::read_docx(path = template) #init new docx object
#  styles_info(docx)  #print stylings
  }
  if(formatUse[3]) { 
    marg0 = DOCoptions$margin
    rtf <- rtf::RTF(reportfn[3], width=DOCoptions$width, height=DOCoptions$height, omi=c(marg0, marg0, marg0, marg0), font.size=DOCoptions$fontSize)  # this can be an .rtf or a .doc
  }
  #insert main header (always first):
  insTitle( paste( L$reportforcase,get("caseID",envir=nnTK)) , 1) #insTitle 

  #################################################################
  #OUTER LOOP TO TRAVERSE (NEED FOR PRIORTY-ARRANGEMENT OF LAYOUT)#
  #################################################################
  nReportItems = length(reportitems)
  ordPriority = order(priority) #items ordered by priority
  for(pos in ordPriority) { #traverse all ordered items and indicate position
    
    if(priority[pos]==0) next #skip if zero

    switch(pos,
           
      #1: HEADER    
      {headTxt = c(paste0("CaseSolver ",L$version,": ",version," (euroformix_",packageVersion("euroformix"),")"),
                    R.version.string,
                    paste0( L$user,": ",Sys.getenv("USERNAME")),
                    paste0( L$created,": ",Sys.time()))
      insText( headTxt ,italic=T)},
        
      #2: References (known)
      insTable(refTabKnown, reportitems[2]),
    
      #3, References (extracted)
      insTable(refTabExtracted, reportitems[3]),

      #4: Single source profiles (alleles)
      insTable(ssTab, reportitems[4]),

      #5: Mix profiles (alleles)
      insTable(mixTab, reportitems[5]),

      #6: Show consensus profiles (alleles)
      {consDataOUT = NULL
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
      }
      insTable(consDataOUT, reportitems[6])},

   
      ###############
      ###SHOW W/PH###
      ###############
  
      #7: Single source profiles (w/PH)
      {selected = "NONE"
      if(!is.null(ssTab)) selected = ssTab[,1] #already sorted
      insList(allTabLIST, reportitems[7], selected)},
  
      #8:  Mixture profiles (w/PH)
      {selected = "NONE"
      if(!is.null(mixTab)) selected = mixTab[,1]  #already sorted
      insList(allTabLIST, reportitems[8], selected)},

      #9; Metadata
      {selected = "NONE"
      if(!is.null(metaDataList) && length(metaDataList)>0 ) selected = names(metaDataList) #must contain elements
      insList(metaDataList,reportitems[9], selected)},

      #################
      ###COMPARISONS###
      #################
      #if(any(checked[9:12])) insTitle(  L$comparisons, 1 )

      #10: Provide match matrix:
      {insTitle( reportitems[10], 2)
      if( !is.null(resCompMAC) ) { #if completed
        if(nrow(resCompMAC)>0) {
          ind <- as.numeric(resCompMAC)<setupThresh$MAC #get smaller indices
          resCompMAC[ind] <- "" #show only greater than threshold
          insTable(resCompMAC,type=4)
        }
      } else {
       insText( L$notcompleted,italic=T)
      }},
    
      #11: Provide match list 1 (qual based):
      {insTitle( reportitems[11], 2)
      if(!is.null(resCompLR1) ) { #if completed
        resCompLR1 <- resCompLR1[as.numeric(resCompLR1[,4])>=log10(setupThresh$LRthresh1),,drop=FALSE]
        ordComp1 = casesolver::orderTableSort(resCompLR1[,1],resCompLR1[,2],sortTypes[3]) #obtain selected sorted order for MatchListQual
        insTable(resCompLR1[ordComp1,,drop=FALSE],type=0)
      } else {
        insText( L$notcompleted,italic=T )
      }},
    
      #12: Provide match list 2 (quan based):
      {insTitle( reportitems[12], 2)
      if(!is.null(resCompLR2) ) { #if completed
        resCompLR2 <- resCompLR2[as.numeric(resCompLR2[,4])>=log10(setupThresh$LRthresh2),,drop=FALSE]
        ordComp2 = casesolver::orderTableSort(resCompLR2[,1],resCompLR2[,2],sortTypes[4]) #obtain selected sorted order for MatchListQuan
        insTable(resCompLR2[ordComp2,,drop=FALSE],type=0)
      } else {
        insText( L$notcompleted,italic=T )
      }},

      #13: Provide match list Final:
      {insTitle( reportitems[13], 2)
      if(!is.null(allMixList) && nrow(allMixList)>0) {
        ordMatches = casesolver::orderTableSort(allMixList[,1],allMixList[,2],sortTypes[5]) #obtain selected sorted order for Matches
        insTable(allMixList[ordMatches,,drop=FALSE],type=0)
      } else {
        insText( L$none )
      }},

      #14: Match network fig is generated as a separet file
      {insTitle( reportitems[14], 2)
      okplot <- !is.null(resCompLR)
      if(okplot) {  #only add if OK
        netf <- file.path(path2,"matchnetwork.png") #file of picture 
        png(netf,width = WHsize[2], height = WHsize[2],res=WHsize[3])
        showMatchNetwork(nnTK,"all",createInteractive=FALSE,selList)
        dev.off()
        insIMG(netf,quadratic = TRUE)
      } else { 
        insText( L$notcompleted )
      }},
    
      ##############
      ###RMP/RMNE###
      ##############
      
      #store random match prob results (to be shown in report)
      #if(any(checked[14:15])) insTitle( L$randommatchprob, 1)
      #15: RMNE
      insTable(resRMP$evid,reportitems[15]),
      
      #16: RMP
      insTable(resRMP$ref,reportitems[16]),
      
      ##################
      ###IBS/evidConc###
      ##################
    
      #17: Evidence concordance 
      insTable(resEvidConc, reportitems[17]),

      #18:  IBS table
      insTable(resIBS, reportitems[18]),
    
      ##################
      #ADVANCED RESULTS#
      ##################
      
      #19: DECONVOLUTION
      {insList(DCdataList, reportitems[19],type=NULL)
      if(length(DCdataList)>0)  insText(DCtext)},  #include text only if results
      
      #20: Weight of evidence results: Table
      insTable(WOEdataTABLE,reportitems[20]), #ADDING #ID FIRST

      #21: Statement results: List
      {insList(WOEstateList,reportitems[21]) #traverses each element in list 
      if(length(WOEstateList)>0) insText(WOEtext)}, #Put woetext at end
      
      #22: Parmaeter results: List
      insList(WOEparamList, reportitems[22],type=NULL), #traverses each element in list 

      #23: Weight of evidence results: LR per marker
      insTable(WOEmarkerTABLE,reportitems[23]),  #ADDING #ID FIRST
      
    ##############
    ###Settings###
    ##############
    
      #24: settings 
      {insTitle( reportitems[24],1 )
      insTitle( L$threshs,2 )
      insText( paste0( L$macthreshold ,colonsymbol,setupThresh$MACthresh) )    
      insText( paste0( L$qualLRthreshold ,colonsymbol,setupThresh$LRthresh1) )
      insText( paste0( L$quanLRthreshold ,colonsymbol,setupThresh$LRthresh2) )
      insText( paste0( L$minLocSSmatch ,colonsymbol,setupThresh$minLociSS) )
      insText( paste0( L$minIBSrelative ,colonsymbol,setupThresh$minIBS) )
      insText( paste0( L$probRatioToNext ,colonsymbol,setupThresh$ratio) )
      insText( paste0( L$probSingleAllele ,colonsymbol,setupThresh$probA) )

      #Obtain model params
      threshT = setupModel$threshT #analytical/detection threshold
      dropinC = setupModel$dropinC #dropin probability
      dropinL = setupModel$dropinL #dropin PH, lambda 
      fst = setupModel$fst #fst correction

      markers = NULL #reset
      if(!is.null(setupMarkers)) {
        vec = function(x) paste0(x,collapse="/")
        markers = vec(setupMarkers[[1]]) #obtain markers
        threshT = vec(setupMarkers[[2]]) #analytical/detection threshold
        dropinC = vec(setupMarkers[[3]]) #dropin probability
        dropinL = vec(setupMarkers[[4]]) #dropin PH, lambda 
        fst = vec(setupMarkers[[5]]) #fst correction
      }
      
      #Include params
      insTitle( L$Modelparameters , 2)
      if(!is.null(markers)) insText( paste0( L$Markers, colonsymbol,markers) )
      insText( paste0( L$analyticalthreshold, colonsymbol,threshT) )
      insText( paste0( L$dropinprob ,colonsymbol,dropinC) )
      insText( paste0( L$dropinpeakheightlambda, colonsymbol,dropinL) )
      insText( paste0(L$fst,colonsymbol,fst) )
      
      #Inlcude population frequency settings
      insTitle(  L$popfreq , 2)
      popfile = basename(setupPop$popfile) #obtain basename of selected population file 
      popfile = strsplit(popfile,"\\.")[[1]][1] #remove extention
      insText( paste0( L$file ,colonsymbol, popfile) )
      insText( paste0( L$AMELincluded ,colonsymbol, ifelse(setupPop$amel=="TRUE",L$yes,L$no) ) )
      insText( paste0( L$Normalized ,colonsymbol, ifelse(as.logical(setupRare$normalize),L$yes,L$no)) )
      if(!is.na(setupRare$minFreq) && setupRare$minFreq!="") insText( paste0( L$minFreq ,colonsymbol, setupRare$minFreq) ) #include minimum freq if set
      
      #Settings or quantitative model:
      insTitle( L$quanmodel , 2)
      kit1 = kit0
      if(is.null(kit1) || kit1=="") kit1 = L$none #indicate if not selected kit
      insText( paste0( L$kit ,colonsymbol,kit1) )
      insText( paste0( L$degradationmodel, colonsymbol, radiotxt[setupModel$degrad] ) )
      insText( paste0( L$BWstuttermodel,colonsymbol, radiotxt[setupModel$stuttBW]) )
      insText( paste0( L$FWstuttermodel,colonsymbol, radiotxt[setupModel$stuttFW]) )
      
      #ADD MCMC SETTINGS:
      if(showItem[2]) {
        insTitle( paste(L$mcmc ,L$settings), 2)
        for(item in names(setupMCMC)) insText( paste0( L[[item]], colonsymbol,setupMCMC[[item]]) ) #include all setings
      }
       
      },  #END SETTINGS
      
    
    ################
    ###ATTACHMENT###
    ################
    #if(any(checked[18:19]))  insTitle( L$attachments , 1)

    #25: Plot EPG figures for single sources
    if( printEPG || sampleType=="LUS" ) { 
      insTitle( reportitems[25], 1)
      if( nrow(ssDataTABLE)>0 ) {
        unREF = unique(ssDataTABLE[,1]) #get unique refs
        refL = getRefL(unREF,forceDi=FALSE) # get relevant references
        
        for(i in 1:nrow(ssDataTABLE)) { #for each single source profiles
          evid <- rownames(ssDataTABLE)[i]
          ref <- ssDataTABLE[i,1]
          
          if(evid == L$empty) {
            insText( L$none )
          } else {
            condref <- refL[ref] #extract reference
            if(length(condref)==0) condref = NULL
            
            tryCatch({ suppressWarnings({
              epgf <- file.path(path2,paste0("epg_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
              png(epgf ,width =WHsize[1], height = WHsize[2],res=WHsize[3])
              if(sampleType=="EPG") euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE,threshT=setupModel$threshT)
              if(sampleType=="LUS") euroformix::plotLUS(mixLIST[evid],sn=evid,condref,threshT=setupModel$threshT,LUSsymbol=LUSsymbol)
              dev.off()
              
              insIMG(epgf) #Insert fig to report
              
            }) })
          } #end if not empty
          insTitle( paste0("#",i," - ",evid), 3)
        } #end for each samples
      } else {
        insText( L$none )
      }#end if
    }, #end if plot EPG
    
    #26: Plot EPG figures for mixtures?
    if((printEPG || sampleType=="LUS")) { 
      insTitle( reportitems[26], 1)
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
            png(epgf ,width =WHsize[1], height = WHsize[2],res=WHsize[3])
            if(sampleType=="EPG") euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE,threshT=setupModel$threshT)
            if(sampleType=="LUS") euroformix::plotLUS(mixLIST[evid],sn=evid,condref,threshT=setupModel$threshT,LUSsymbol=LUSsymbol)
            dev.off()
           
            insIMG(epgf) #Insert fig to report
            insTitle( paste0("#",i," - ",evid), 3)
          }) })
        } #end for each samples
      } else {
        insText( L$none )
      }#end if
    } #end if plot EPG
   ) #end switch-case
  } #end outer for-loop
##################################################################################
  
  print("REPORTS STORED IN:")
  print(path)
  if(formatUse[1] && length(createHTML)==2 ) {
    if(createHTML[2]) browseURL(htmlf) #open file (preview)
  }
  if(formatUse[2]) {
    #officer::body_end_section_landscape(docx)
    print(docx, target = reportfn[2])  #saving to file
    if(length(createDOCX)==2 && createDOCX[2]) browseURL(reportfn[2]) #open file (preview)
  }
  if(formatUse[3]) { 
    rtf::done(rtf)
    if(length(createDOC)==2 && createDOC[2]) browseURL(reportfn[3]) #open file (preview)
  }
} 