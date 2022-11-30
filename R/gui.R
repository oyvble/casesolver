#' @title gui
#' @description gui (GUI) is a GUI wrapper for CaseSolver
#' @details The function starts the CaseSolver software.
#' @param envirfile A file to a saved environment of a project (must contain nnTK)
#' @param envir A saved environment (nnTK)
#' @export

gui = function(envirfile=NULL, envir=NULL) {
  LUSsymbol = "_"
  MPSsymbol = ":" #Added in version 1.5.0
  colonsymbol = ":" #use variable for colon
  #noKit = "NONE" #name of empty kit
  defaultEncoding =  options()$encoding #Set default to local encoder: Makes it possible to use Language-text files in own language
  defaultLanguage = "English"
  
  #size of main window
  mwH <- 500
  mwW <- 1000
  w0 <- 20 #default textwidth (for gWidgets2::gedit)
  nLarge = 10000 #this is large threshold of number of references (dont change table alignment if more than this number)
  
  #type of gwidgets-kit
  options(guiToolkit="tcltk")
  
  #models:
  version =  packageVersion("casesolver") #obtain version
  spc <- 10  #Spacing between widgets
  
  #REPORT ITEMS (may change over versions)
  
  nReportItems = 26 #number of report items (will be inserted when opening program)
  reportitems <- rep("",nReportItems)  #assign emtpy  report items (must correspond to size of reportitems)
  
  #########################################################
  ################Environment variables####################
  #########################################################
  
  #NB: pgkPath and .sep must be changed before compiling R-package!
  pgkPath <- path.package("casesolver", quiet = FALSE) # Get package path.
  #pgkPath <- "config"# Get package path (possible to synchronize)
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  
  setupRead = function(file) scan(file,what=character(),quiet=TRUE,sep="\n") #skip line is new element
  setupWrite = function(vals,file) write(vals,file=file,sep="\n")   #save to file in installation folder
  setupFileThresh <- paste(pgkPath,"configThresh",sep=.sep)   #file for threshold selection
  setupFileModel <- paste(pgkPath,"configModel",sep=.sep)   #file for model selection
  setupFileMarkers <- paste(pgkPath,"configMarkers",sep=.sep)   #file for model selection per markers
  setupFileKit <- paste(pgkPath,"configKit",sep=.sep)  #file for kit selection
  setupFilePop <- paste(pgkPath,"configPop",sep=.sep)  #file for population selection
  setupFileRare <- paste(pgkPath,"configRare",sep=.sep)  #file for options of rare alleles 
  setupFileCase <- paste(pgkPath,"configCase",sep=.sep)   #file for case selection
  setupFileImport <- paste(pgkPath,"configImport",sep=.sep)  #file for import selection
  setupFileReportLay <- paste(pgkPath,"configReportLay",sep=.sep) #used to set report layout
  setupFileReportOpt <- paste(pgkPath,"configReportOpt",sep=.sep) #used to set report options
  setupFileReportExpTyp <- paste(pgkPath,"configReportExpType",sep=.sep) #used to set report export type
  setupFileReportExpOpt <- paste(pgkPath,"configReportExpOpt",sep=.sep) #used to set report export option
  setupFileReportLocNames <- paste(pgkPath,"configReportLocNames",sep=.sep) #used to set marker names in report
  
  setupFileMCMC <- paste(pgkPath,"configMCMC",sep=.sep) #used to set advanced settings
  setupFileSorting <- paste(pgkPath,"configSorting",sep=.sep) #used to set sorting of tables settings
  setupFileAdvanced <- paste(pgkPath,"configAdvanced",sep=.sep) #used to set advanced settings
  setupFileView <- paste(pgkPath,"configView",sep=.sep) #used to set GUI layout
  setupFileLanguage <- paste(pgkPath,"configLanguage",sep=.sep)  #file for Language selection
  setupFileExport <- paste(pgkPath,"mmTK.Rdata",sep=.sep)  #used to run euroformix
  setupFileExport2 <- paste(pgkPath,"mmTK2.Rdata",sep=.sep)  #used to run euroformix
  
  #The files are stored in system settings and loaded when opening the tool:
  
  #Thresholds:
  optF = c(0.8,10,1000,7,14,15,0.99)
  if(file.exists(setupFileThresh)) optF <- setupRead(file=setupFileThresh)
  setupThresh = list(MACthresh=as.numeric(optF[1]),LRthresh1=as.numeric(optF[2]),LRthresh2=as.numeric(optF[3]),minLociSS=as.integer(optF[4]),minIBS=as.integer(optF[5]),ratio=as.numeric(optF[6]),probA=as.numeric(optF[7]))
  
  #Model options (global):
  optF =  c(50,0.05,0.01,0.01,1,2,2,1) #threshT: default value of detection threshold value,dropinC: Dropin probability parameter in model,dropinL: Dropin peak height distribution parameter Lambda. Modeltype={1,2,3}={qual,quan,both}
  if(file.exists(setupFileModel)) optF <- setupRead(file=setupFileModel)
  setupModel= list(threshT=as.numeric(optF[1]),dropinC=as.numeric(optF[2]),dropinL=as.numeric(optF[3]),fst=as.numeric(optF[4]),degrad=as.integer(optF[5]),stuttBW=as.integer(optF[6]),stuttFW=as.integer(optF[7]),modeltype=as.integer(optF[8]))
  
  #Model options per marker (need to traverse a flexible number of items):
  setupMarkers <- NULL #default is no specification (using global)
  nItemsMarkers = 5 #marker names, AT values, Dropin probs, Dropin lambdas, Fsts
  if(file.exists(setupFileMarkers)) {
    optF <- setupRead(file=setupFileMarkers)
    if(length(optF)>0) { #if registered
      nLocs = round(length(optF)/nItemsMarkers) #number of loci in file
      setupMarkers = list()#obtain vector with names and values 
      for(itemInd in 1:nItemsMarkers) { #traverse through each elemen type
        readRange = nLocs*(itemInd-1) + 1:nLocs
        insVal = optF[readRange] #extend existing vector with value
        if(itemInd>1) insVal = as.numeric(insVal)
        setupMarkers[[itemInd]] = insVal
      }
    }
  }
  
  #costumized locus/marker names in report
  setupReportLocNames = NULL
  if(file.exists(setupFileReportLocNames)) {
    optF <- setupRead(file=setupFileReportLocNames)
    if(length(optF)>0) nLocs = round(length(optF)/2) #use number of loci in file
    setupReportLocNames = optF[1:nLocs + nLocs] #obtain edited marker names
    names(setupReportLocNames) = optF[1:nLocs] #obtain Conventional names
  }
  
  #Kit selection
  optF = "" #no kit selected by default
  if(file.exists(setupFileKit)) optF <- setupRead(file=setupFileKit)
  setupKit = list(kitname=optF[1])   
   
  #Population freq (file selection):
  optF = rep("",2) #filename, amelChoice
  if(file.exists(setupFilePop)) optF <- setupRead(file=setupFilePop)
  setupPop = list(popfile=optF[1],amel=optF[2])
   
  #Population freq (rare alleles):
  optF = c("TRUE","") #whether to normalize, min. frequencies
  if(file.exists(setupFileRare)) optF <- setupRead(file=setupFileRare)
  setupRare = list(normalize=optF[1],minFreq=optF[2])
  
  #Path to Case folder (selection)
  optF = "" #empty casepath by default. Points to a folder with CASEID given as folder names. 
  if(file.exists(setupFileCase)) optF <- setupRead(file=setupFileCase)
  setupCase = list(casepath=optF[1]) 
  
  #Path to importData R-file (must contain the R-function importData. Points to a text-file.)
  optF = "" #empty importfile by default (must contain the R-function importData. Points to a text-file.
  if(file.exists(setupFileImport)) optF <- setupRead(file=setupFileImport)
  setupImport = list(importfile=optF[1])
  
  #Advanced options
  optF = c(4,3,3,"TRUE","FALSE","FALSE","FALSE")
  if(file.exists(setupFileAdvanced)) optF <- setupRead(file=setupFileAdvanced)
  setupAdvanced = list(maxC1=as.integer(optF[1]),maxC2=as.integer(optF[2]),nDone=as.integer(optF[3]),useMinK1=optF[4],compSS=optF[5],isSNP=optF[6],selProfiles=optF[7])   
  
  #Data view (vertical or horizontal)
  optF = c("FALSE")
  if(file.exists(setupFileView)) optF <- setupRead(file=setupFileView)
  setupView = list(importHorizontal=optF[1])
  
  #selected language
  optF = defaultLanguage #,encoding=defaultEncoding) #English langugage used by default.  
  if(file.exists(setupFileLanguage)) optF <- setupRead(file=setupFileLanguage)
  setupLanguage = list(language=optF[1])#,encoding=optF[2]) 
  
  #Option for sorting of tables:
  #Tables in import: 1) Evid (SingleSources and Mixtures) [EPGs will follow same order] and  2) Refs (known and extended)
  #Tables from comparison: 3) MatchList(Qual), 4) MatchList(Quan), 5) Matches
  optF = rep(1, 5) #defualt sorting
  if(file.exists(setupFileSorting)) optF <- setupRead(file=setupFileSorting)
  setupSorting = as.integer(optF)
  
  #Options for MCMC
  optF = c(2000,2,0.05,1)
  if(file.exists(setupFileMCMC)) optF <- setupRead(file=setupFileMCMC)
  setupMCMC = list(niter=as.integer(optF[1]),delta=as.numeric(optF[2]),quantile=as.numeric(optF[3]),seed=as.integer(optF[4]))   
  
  #Report:
  #Layout option
  optF = 1:nReportItems #order considered (0 means not used)
  if(file.exists(setupFileReportLay)) optF <- setupRead(file=setupFileReportLay)
  setupReportLay = list(priority=as.integer(optF)) #Gives priority of how layout should be given
  
  #Other options (whether to include in report or not)
  optNames = c("MatchStatus","MCMCsettings","mleLR","bayesLR","consLR","Mx","validFailed","verbalLR","headerTime")
  setupReportOpt <- optF <- setNames( rep(TRUE,length(optNames)),optNames) #WHether to show different options
  if(file.exists(setupFileReportOpt)) {
    optF2 <- setupRead(file=setupFileReportOpt)
    setupReportOpt[seq_along(optF2)] <- as.logical(optF2) #handling special case
  }
    
  #Export options: Export/Preview
  optF = c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE) #default values to insert
  if(file.exists(setupFileReportExpTyp)) optF <- as.logical( setupRead(file=setupFileReportExpTyp) )
  setupReportExpTyp = list(HTML=optF[1:2],DOCX=optF[3:4],DOC=optF[5:6]) 
  
  #Other export options (Table, sizes), DOC
  optF = c(1920,1080,120, 11,8.5,0.1,6,9,  6,11)
  if(file.exists(setupFileReportExpOpt)) optF <- setupRead(file=setupFileReportExpOpt)
  names(optF) = c("Width (px)","Height (px)","Resolution (ppi)",   #image options
                 "Width (inch)","Height (inch)","Margin (inch)","Table size","Text font size", #DOC options
                 "Table size","Text font size") #DOCX options
  setupReportExpOpt = list(PNG=optF[1:3], DOC=optF[4:8],DOCX=optF[9:10]) 
  
  if(is.null(envir)) { #Need to create invironment
    #Always done: create new envornment object. Parent must be empty
    nnTK = new.env( parent = emptyenv() ) 
    
    #Toolbar options: can be changed any time by using toolbar
    assign("setupThresh",setupThresh,envir=nnTK) 
    assign("setupModel",setupModel,envir=nnTK) 
    assign("setupMarkers",setupMarkers,envir=nnTK)  
    assign("setupKit",setupKit,envir=nnTK) 
    assign("setupPop",setupPop,envir=nnTK)
    assign("setupRare",setupRare,envir=nnTK)
    assign("setupCase",setupCase,envir=nnTK) 
    assign("setupImport",setupImport,envir=nnTK) 
    assign("setupAdvanced",setupAdvanced,envir=nnTK) 
    assign("setupView",setupView,envir=nnTK) 
    assign("setupLanguage",setupLanguage,envir=nnTK) 
    assign("setupMCMC",setupMCMC,envir=nnTK) 
    assign("setupSorting",setupSorting,envir=nnTK) 
    
    #Report settings
    assign("setupReportLay",setupReportLay,envir=nnTK)  #layout options
    assign("setupReportExpTyp",setupReportExpTyp,envir=nnTK)  #Export types
    assign("setupReportExpOpt",setupReportExpOpt,envir=nnTK)  #export format
    assign("setupReportOpt",setupReportOpt,envir=nnTK)  #other options (also whether to rotate tables)
    assign("setupReportLocNames",setupReportLocNames,envir=nnTK)  #User may modify marker names
    
    #initializing environment variables
    assign("workdir",NULL,envir=nnTK) #assign working directory to nnTK-environment
    
    #Considered CaseID:
    assign("caseID",NULL,envir=nnTK) #Used for report
    assign("version",version,envir=nnTK) #Include version which is used (useful for later)
    
    #imported data:
    assign("popFreq",NULL,envir=nnTK) #assign to nnTK-environment
    assign("mixDataTABLE",NULL,envir=nnTK) #Table of evidence profiles (only alleles)
    assign("refDataTABLE",NULL,envir=nnTK) #Table of ref profiles (only alleles)
    assign("mixDataMATCHSTATUS",NULL,envir=nnTK) #vector with match status for evidence profiles
    assign("mixDataLIST",NULL,envir=nnTK)  #list of evidence profiles (alleles and heights)
    assign("metaDataLIST",NULL,envir=nnTK) #list of imported metadata
    assign("consDataTABLE",NULL,envir=nnTK) #list of imported consensus data
    
    #table click: Program remembers last settings
    assign("clicktableLast",NULL,envir=nnTK) #used for table click when visualizing epgs
    
    #results: stored results after calculations
    assign("resCompMAC",NULL,envir=nnTK)  #store MAC results from comparison
    assign("resCompLR1",NULL,envir=nnTK)  #store LR results from qual based comparison (sorted)
    assign("resCompLR2",NULL,envir=nnTK)  #store LR results from quan comparison (sorted)
    assign("resCompLR",NULL,envir=nnTK)  #store LR results from qual/quan comparison (unsorted) 
    assign("resMatches",NULL,envir=nnTK)  #store match-results from comparison (those with LR>threshold)
    assign("allMixList",NULL,envir=nnTK)  #store match-results from comparison (those with LR>threshold) together with all mixtures
    assign("DClist",NULL,envir=nnTK)  #store list with deconvoluted reference candidates
    assign("DClistReport",NULL,envir=nnTK)  #store list with confirmed deconvoluted reference candidates (includes probs)
    assign("resRMP",NULL,envir=nnTK)  #store random match prob results 
    assign("resIBS",NULL,envir=nnTK)  #store IBS compare results (used in report)
    assign("storedFitHp",NULL,envir=nnTK) #store mle-fitted objects under Hp (for quantitative LR matchlist)
    assign("resEvidConc",NULL,envir=nnTK)  #store concordance results
    
    #Object for WoE calculations
    assign("modelSettings",NULL,envir=nnTK)  #currently not in use: would include following objects: popFreq,kit,fst,lambda,threshT,xiBW,xiFW,pC (updated when saved)
    assign("resWOEeval",NULL,envir=nnTK)  #store match-results from weight of evidence evaluations (LR-per marker, PH-validation)
    
    nnTK0 = nnTK #store environment to make older projects backward compatible
  } else if(is.null(envirfile)) {
     nnTK = envir #copy environment
  }
    
  #} else { #IF FILE IS GIVEN
  if(!is.null(envirfile)) {
    load(envirfile) #loading environment
    
    #Inform the user that the old project file is not supported by current version of casesolver
    outdated=FALSE
    projVersion = nnTK$version #obtain project version get("version",envir=nnTK)
    if(is.null(projVersion)) {
      outdated = TRUE
    }  else {
      projVersion2 = unlist(projVersion)[1:2] #strsplit(projVersion,"\\."))[1:2]
      version2 = unlist(version)[1:2]
      if( is.null(projVersion) || version2[1]!=projVersion2[1] || version2[2]!=projVersion2[2]) outdated = TRUE
    }
    if(outdated) { #notify if stored project version is outdated
      versionTxt = "The loaded project file was outdated:\n"
      versionTxt = paste0(versionTxt,"Project version: ",projVersion)
      versionTxt = paste0(versionTxt,"\n\nPlease create a new project file with current CaseSolver version.\n\nWarning: The program may not work properly with an outdated project file.")
      gWidgets2::gmessage(versionTxt,"Outdated project file", icon="info")
     
      #Try make loaded project as backward compatible as possible: Use 'default settings' if non-compatible objects
      reqObjs = names(nnTK0) #obtain (all) required objects
      oldObjs = names(nnTK)  #objects from old version
      for(obj in reqObjs) {  #traverse all objects
        if(!obj%in%oldObjs) {
          assign(obj , nnTK0[[obj]] ,envir=nnTK)  #store default object if not found
        } else {
          reqItems = names(nnTK0[[obj]]) #obtain required items
          if(is.null(reqItems)) next #jump to next object
          
          oldItems = names(nnTK[[obj]]) #obtain required items
          if(length(reqItems)!= length(oldItems) || !all(reqItems==oldItems)) {
            assign(obj , nnTK0[[obj]] ,envir=nnTK)  #store default object if missmatch found
          }
        }
      }
    } #end if outdated
    #assign("setupRare",setupRare,envir=nnTK)
  } #end if project file was restored
  
  #################################################################################
  ###########################GUI HELPFUNCTIONS#####################################
  #################################################################################
  
  #This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
  mygfile <- function(text,type,filter=list(),initf=NULL) { #Possible bug: text ignored when type="selectdir"
   file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
   Encoding(file ) <- defaultEncoding #Set to local encoder L$encoding #set special encoding (latin1 is default)
   return(file)
  }
  
  #Helpfunction to check if file is OK
  fileNotOK = function(x) { return( is.na(x) || length(x)==0)  }
   
  #Helpfunction to save a table to text file
  saveTable = function(tab,sep="txt") {
    tabfile  <- mygfile(text= paste( L$save , L$table) ,type="save") #csv is correct format!
    if(fileNotOK(tabfile)) return()
    
    if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
    if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
    if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=";",row.names=FALSE) 
    print(paste("Table saved in ",tabfile,sep=""))
  } #end file
  
  #Helpfunction to set pop freq to environment 
  setPopFreq = function(change=FALSE,giveMessage=TRUE) { #helpfunction to read popFreq from file and set to environment
    if(!change && !is.null( get("popFreq",envir=nnTK))) return(TRUE) #return if already set
    
    opt <- get("setupPop",envir=nnTK) 
    tryCatch( {
      popFreq <- getFreqs(opt$popfile)
      AMEL <- c(0.75,0.25) #Assuming 50-50 Male/Femal population, abuse on Y/Y possibility
      names(AMEL) <- c("X","Y")
      if( opt$amel=="TRUE" && !any(grepl("AM",names(popFreq))) ) popFreq$AMEL =AMEL
      assign("popFreq",popFreq,envir=nnTK) #assign popFreq to nnTK-environment
      #if(verbose) print(popFreq) #print first time
    }, error = function(e) return(FALSE) )
    if(is.null(get("popFreq",envir=nnTK))) { #if popFreq is missing
      if(giveMessage) gWidgets2::gmessage( paste( L$msg.setPopFreq, L$setup,">",L$popfreq  ) ) 
      return(FALSE) 
    }
    return(TRUE)
  } 
  
  #Helpfunction to show table in GDF (editable cells)
  showGDFtable = function(title,table) {
    setwin2 <- gWidgets2::gwindow( title ,visible=FALSE) 
    guitab <- gWidgets2::gdf(items=table,container = setwin2) 
    gWidgets2::visible(setwin2) <- TRUE
    gWidgets2::focus(setwin2) <- TRUE
  }
  
  #Helpfunction to ensure focus of main window
  setFocus = function() { gWidgets2::focus(mainwin) <- TRUE }

  ###########################FILE#####################################
  f_setwd = function(h,...) {
    dirfile = mygfile(text= paste( L$select , L$folder) ,type="selectdir")
    if(fileNotOK(dirfile)) return()
    setwd(dirfile)
    assign("workdir",dirfile,envir=nnTK) #assign working directory
  }
  
  f_openproj = function(h,...) {
    filterList = list()
    filterList[[ L$proj  ]] = list(patterns=list("*.Rdata")) #set file extension pattern 
    filterList[[ L$all  ]] = list(patterns=list("*")) #set file extension pattern 
    projfile = mygfile(text= paste( L$open , L$proj) ,type="open", filter=filterList)
    if(fileNotOK(projfile)) return()
    gWidgets2::dispose(mainwin)
    casesolver::gui(projfile) #open environment file in program
  }
  
  f_saveproj = function(h,...) {
    projfile = mygfile(text= paste( L$save , L$proj) ,type="save")
    if(fileNotOK(projfile)) return()
    if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata")
    print("Size of stored objects (in MB):") #prints size of each stored object
    print(sapply(nnTK,object.size)/1e6) #prints size of each stored object
    save(nnTK,file=projfile,compress="xz",eval.promises=FALSE,precheck=FALSE,compression_level=2)
    print(paste("Project saved in ",projfile,sep=""))
  }
  
  f_quitproj = function(h,...) {
    ubool <- gWidgets2::gconfirm( L$msg.saveproj ,title= L$quitprog ,icon="info")
    if(ubool) {
      f_saveproj(h)
    } else { 
      print("Program terminated without saving")
    }
    gWidgets2::dispose(mainwin) #remove window!
  }
  
  ###########################REPORT SETTINGS#####################################
  #Selection of report layout
  f_reportlay = function(h,...) { #GUI function to set report layout
    opt <- get("setupReportLay",envir=nnTK) #get layout settings
    setwin <- gWidgets2::gwindow( L$reportlayout ,visible=FALSE)
    group = gWidgets2::ggroup(spacing=5,container= setwin,horizontal=FALSE) 
    checkbox = gWidgets2::glayout(spacing=0,container=(group),horizontal=FALSE) 
    reportitems2 = paste( L$show , reportitems) #put "show" in front of each reportitems
    nItems = length(reportitems) #number of items
    itemRange = 0:nItems
    checkbox[1,1] = gWidgets2::glabel("Report object",container=checkbox)
    checkbox[1,2] = gWidgets2::glabel("Priority",container=checkbox)
    for(i in 1:nItems) {
      val = opt$priority[i] #obtain value
      checkbox[i+1,1] <- gWidgets2::gcheckbox( reportitems2[i] ,checked=val>0,container=checkbox)
      checkbox[i+1,2] <- gWidgets2::gcombobox(items=itemRange,selected=which(val==itemRange),editable=TRUE,container=checkbox)
      gWidgets2::size(checkbox[i+1,2]) = 4
    } 
    #gWidgets2::add(tabval,checkbox)#,expand=TRUE,fill=TRUE)
    savebut <- gWidgets2::gbutton( L$save ,use.table=FALSE,container=group,handler = function(h, ...) { 
      for(i in 1:nItems) {
       val = 0 #default is not checked
       checked = gWidgets2::svalue(checkbox[i+1,1]) #indicate if checked
       if(checked) val= as.integer(gWidgets2::svalue(checkbox[i+1,2])) #obtain value
       opt$priority[i] = val #update priority variable
      } 
      assign("setupReportLay",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt$priority),file=setupFileReportLay)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    })
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #Selection of report options
  f_reportopt = function(h,...) { #GUI function to set report option
    opt <- get("setupReportOpt",envir=nnTK) #get layout settings
    #obtain names of options (INSERT HERE ONLY, NOT USED IN REPORT)
    items = paste(L$show,names(opt)) #c(L$matchstatus) #"MatchStatus"  
    setwin <- gWidgets2::gwindow( paste(L$report, L$options)  ,visible=FALSE)
    tabtmp <- gWidgets2::ggroup(horizontal = FALSE,container=setwin)
    
    #Option layout
    checkBox = gWidgets2::glayout(spacing=0,container= tabtmp)  
    for(j in 1:length(items)) {
      checkBox[j,1] = gWidgets2::glabel(names(items)[j],container=checkBox)
      checkBox[j,2] = gWidgets2::gcheckbox(items[j],checked = as.logical(opt[j]),container=checkBox)
    }
    
    #Save-button:
    savebut = gWidgets2::gbutton(L$save,container=tabtmp, handler = function(h) {
      for(j in 1:length(items)) {
         val = as.logical( gWidgets2::svalue(checkBox[j,2]) )
         if(!is.na(val)) opt[j]= val #insert value if not NA
      }
      assign("setupReportOpt",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileReportOpt)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    })
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #Selection of (exported) report types 
  f_reportexptyp = function(h,...) { #GUI function to set report export options
    opt <- get("setupReportExpTyp",envir=nnTK) #get layout setting
    setwin <- gWidgets2::gwindow( paste(L$report, L$export)  ,visible=FALSE, width=300,height=200)
    checkbox = gWidgets2::glayout(spacing=5,container=(setwin),horizontal=FALSE) 
    
    #Showing selection tables:
    types = names(opt) #c("HTML","DOCX","DOC") #report types
    items = c(L$report, L$preview) 
    for(i in 1:length(types)) checkbox[i+1,1] = gWidgets2::glabel(types[i],container=checkbox)
    for(j in 1:length(items)) checkbox[1,j+1] = gWidgets2::glabel(items[j],container=checkbox)
    for(i in 1:length(types)) {
      type = types[i]
      for(j in 1:length(items)) {
       val=opt[[type]][j] #obtain value
       checkbox[i+1,j+1] <- gWidgets2::gcheckbox("",checked=val,container=checkbox)
      }
    }
    #Save-button:
    checkbox[1,1] = gWidgets2::gbutton(L$save,container=checkbox, handler = function(h) {
      for(i in 1:length(types)) {
       type = types[i]
       for(j in 1:length(items)) {
         opt[[type]][j]=  gWidgets2::svalue(checkbox[i+1,j+1]) #obtain value
       }
      }
      assign("setupReportExpTyp",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileReportExpTyp)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    })
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #export options
  f_reportexpopt = function(h,...) { #GUI function to set report option
    opt <- get("setupReportExpOpt",envir=nnTK) #get layout settings
    types = names(opt)
    setwin <- gWidgets2::gwindow( paste(L$report,  L$export, L$options)  ,visible=FALSE)
    tabtmp <- gWidgets2::ggroup(horizontal = FALSE,container=setwin)
    
    #Option layout
    checkList = list() #layout list
    for(i in 1:length(types)) {
      type = types[i]
      items = opt[[type]]
      checkList[[i]] = gWidgets2::glayout(spacing=0,container= gWidgets2::gframe(type, container=tabtmp))  
      
      for(j in 1:length(items)) {
        checkList[[i]][j,1] = gWidgets2::glabel(names(items)[j],container=checkList[[i]])
        checkList[[i]][j,2] = gWidgets2::gedit(items[j],container=checkList[[i]])
      }
    }
    #Save-button:
    savebut = gWidgets2::gbutton(L$save,container=tabtmp, handler = function(h) {
      for(i in 1:length(types)) {
        type = types[i]
        items = opt[[type]]
        for(j in 1:length(items)) {
          val = as.numeric( gWidgets2::svalue(checkList[[i]][j,2]) )
          if(!is.na(val)) opt[[type]][j]= val #insert value if not NA
        }
      }
      assign("setupReportExpOpt",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileReportExpOpt)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    } )
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #helpfunction to modify marker names for report tables
  f_reportlocnames = function(h,...) {
    locNameVec <- get("setupReportLocNames",envir=nnTK) #obtain stored object (vector of marker names)
  
    #Obtain marker names from current imported case
    evidTable = get("mixDataTABLE",envir=nnTK) #Table of evidence profiles (only alleles)
    refTable = get("refDataTABLE",envir=nnTK) #Table of ref profiles (only alleles)
    
    caseLocNames <- NULL #toupper( names(locNameVec)) #default is markers stored in setupFile
    if(!is.null(evidTable)) { #if evid table found
      caseLocNames = toupper(colnames(evidTable))
    } else if(!is.null(refTable)) { #if else ref table found
      caseLocNames = toupper(colnames(refTable))
    }
    if(is.null(caseLocNames) && is.null(locNameVec)) return() #if no marker info found
    
    #POSSIBLY MODIFY THE locNAme vector    
    if(is.null(locNameVec)) {
      locNameVec = setNames(caseLocNames,caseLocNames) #use markers from imported case
    } else if(!is.null(caseLocNames)) { #If something was previously stored and case info is found
      locNameVec2 = setNames(caseLocNames,caseLocNames) #use markers from imported case
      commonLocs = intersect(caseLocNames,names(locNameVec)) #find common
      locNameVec2[commonLocs] = locNameVec[commonLocs] #update marker names from settings (not from case)
      locNameVec = locNameVec2 #update
    }
    locNames = names(locNameVec) #obtain names
    
    #Create window for selection
    setwin <- gWidgets2::gwindow( paste(L$Markernames)  ,visible=FALSE, width=100, height=500)
    grouplay <- gWidgets2::ggroup(spacing=2,container=(setwin),horizontal = FALSE, use.scrollwindow = TRUE)  #set group layout
    tabsel = gWidgets2::glayout(spacing=10,container=(grouplay),horizontal = TRUE)  #set grid  (will contain buttons)
    tabsel[1,1] = gWidgets2::gbutton(text=L$restore,container=tabsel, handler = function(h,...) {
      for(rowind in 1:length(locNameVec))  gWidgets2::svalue(tabval[rowind,2]) <- locNameVec[rowind] #recover name in edit box
    }) 
    tabsel[1,2] = gWidgets2::gbutton(text=L$save,container=tabsel, handler = function(h,...) {
      for(rowind in 1:length(locNameVec))  locNameVec[rowind] <- gWidgets2::svalue(tabval[rowind,2]) #obtain names from edit box
      assign("setupReportLocNames",locNameVec,envir=nnTK)  #assign user-values
      setupWrite(c(names(locNameVec),locNameVec),file=setupFileReportLocNames)    #save to file in installation folder
      gWidgets2::dispose(setwin) #remove window
    }) 
    
    #Box with names
    tabval = gWidgets2::glayout(spacing=0,container=(grouplay))  #create grid layout
    w0 <- 15 #width of textbox
    for(rowind in 1:length(locNameVec)) {
      tabval[rowind,1] <- gWidgets2::glabel(text=locNames[rowind],container=tabval)  #insert conventional marker name
      tabval[rowind,2] <- gWidgets2::gedit(text=locNameVec[rowind],width=w0,container=tabval)  #insert conventional marker name
    }
    gWidgets2::focus(setwin) = TRUE #set top
    gWidgets2::visible(setwin) = TRUE #show window
  }
  
  
  ###########################SETTINGS#####################################
  
  #GUI function to selecting threshold values
  f_threshsel=  function(h,...) { 
    opt <- get("setupThresh",envir=nnTK) 
    setwin <- gWidgets2::gwindow(paste0(  L$threshsettings ),visible=FALSE)
    tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
    tabval[1,1] <- gWidgets2::glabel(text= paste0( L$macthreshold ),container=tabval) #MATCH THRESHOLD
    tabval[1,2] <- gWidgets2::gedit(text=opt$MACthresh,width=w0,container=tabval)
    tabval[2,1] <- gWidgets2::glabel(text= paste0( L$qualLRthreshold ),container=tabval) #QUAL LR THRESHOLD
    tabval[2,2] <- gWidgets2::gedit(text=opt$LRthresh1,width=w0,container=tabval)
    tabval[3,1] <- gWidgets2::glabel(text= paste0( L$quanLRthreshold ) ,container=tabval) #QUAN LR THRESHOLD
    tabval[3,2] <- gWidgets2::gedit(text=opt$LRthresh2,width=w0,container=tabval)
    tabval[4,1] <- gWidgets2::glabel(text= paste0( L$minLocSSmatch ) ,container=tabval) #MINIMUM NUMBER OF LOCI TO MATCH WITH SINGLE SOURCE
    tabval[4,2] <- gWidgets2::gedit(text=opt$minLociSS ,width=w0,container=tabval)
    tabval[5,1] <- gWidgets2::glabel(text= paste0( L$minIBSrelative ) ,container=tabval) #REQUIRED MATCHING ALLELES BETWEEN REFS
    tabval[5,2] <- gWidgets2::gedit(text=opt$minIBS ,width=w0,container=tabval)
    tabval[6,1] <- gWidgets2::glabel(text= paste0( L$probRatioToNext ) ,container=tabval) #RATIO TO NEXT GENO-PROBABILITY THRESHOLD
    tabval[6,2] <- gWidgets2::gedit(text=opt$ratio ,width=w0,container=tabval)
    tabval[7,1] <- gWidgets2::glabel(text= paste0( L$probSingleAllele ) ,container=tabval) #ALLELE-PROBABILITY THRESHOLD
    tabval[7,2] <- gWidgets2::gedit(text=opt$probA ,width=w0,container=tabval)
    tabval[8,1] <- gWidgets2::gbutton( L$save, container=tabval,handler = function(h, ...) { 
      opt$MACthresh <- as.numeric(gWidgets2::svalue(tabval[1,2])) 
      opt$LRthresh1 <- as.numeric(gWidgets2::svalue(tabval[2,2])) 
      opt$LRthresh2 <- as.numeric(gWidgets2::svalue(tabval[3,2])) 
      opt$minLociSS <- as.numeric(gWidgets2::svalue(tabval[4,2])) #number of loci for being a match
      opt$minIBS <- as.numeric(gWidgets2::svalue(tabval[5,2])) #number of common alleles for being a relative
      opt$ratio <- as.numeric(gWidgets2::svalue(tabval[6,2])) #ratio of probability (deconvolution candidates)
      opt$probA <- as.numeric(gWidgets2::svalue(tabval[7,2])) #ratio of probability (deconvolution candidates)
      assign("setupThresh",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileThresh)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    })
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #GUI for for selecting model settings
  f_modelsel=  function(h,...) {
   modtypetxt <- c( L$qualmodel , L$quanmodel , L$both )
   opt <- get("setupModel",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0(  L$modelsettings ),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
  
   tabval[1,1] <- gWidgets2::glabel(text= paste0(L$modeltypes," (",L$compare,")") ,container=tabval) #Model type(s)
   tabval[1,2] <- gWidgets2::gcombobox(items=modtypetxt,selected=opt$modeltype,horizontal=TRUE,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text= L$analyticalthreshold ,container=tabval) #Detection threshold
   tabval[2,2] <- gWidgets2::gedit(text=opt$threshT,width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text= L$dropinprob ,container=tabval) #Dropin probability
   tabval[3,2] <- gWidgets2::gedit(text=opt$dropinC,width=w0,container=tabval)
   tabval[4,1] <- gWidgets2::glabel(text= L$dropinpeakheightlambda ,container=tabval) #Dropin peak height Lambda (EFM)
   tabval[4,2] <- gWidgets2::gedit(text=opt$dropinL,width=w0,container=tabval)
   tabval[5,1] <- gWidgets2::glabel(text= L$fstsetting,container=tabval) #Dropin peak height Lambda (EFM)
   tabval[5,2] <- gWidgets2::gedit(text=opt$fst,width=w0,container=tabval)
   tabval[6,1] <- gWidgets2::glabel(text= L$degradationmodel ,container=tabval) #Degradation model (EFM)
   tabval[6,2] <- gWidgets2::gradio(items=radiotxt,selected=opt$degrad,horizontal=TRUE,container=tabval)
   tabval[7,1] <- gWidgets2::glabel(text= L$BWstuttermodel,container=tabval) #Stutter model (EFM)
   tabval[7,2] <- gWidgets2::gradio(items=radiotxt,selected=opt$stuttBW,horizontal=TRUE,container=tabval)
   tabval[8,1] <- gWidgets2::glabel(text= L$FWstuttermodel,container=tabval) #Stutter model (EFM)
   tabval[8,2] <- gWidgets2::gradio(items=radiotxt,selected=opt$stuttFW,horizontal=TRUE,container=tabval)
   
   #helpfunction to set marker specific settings
   f_setpermarker = function(h) {
     isok = setPopFreq(giveMessage=FALSE) #set  population frequency from last selected file
     if(!isok) {
       print("Specify a frequency file to continue.")
       return()
     }
     gWidgets2::dispose(setwin) #close setup window
     print("Please be patient (this could take a while)...")
     casesolver::setMarkerSettings(nnTK) # Obtain updated object
     
     #Store object to project and to file:     
     write(unlist(get("setupMarkers",envir=nnTK)),file=setupFileMarkers)    #save to file in installation folder
     gWidgets2::gmessage("Marker settings was successfully stored.")
   }
   
   tabval[2,3] <- gWidgets2::gbutton(text= L$setpermarker ,container=tabval, handler=f_setpermarker) #Detection threshold
   tabval[3,3] <- gWidgets2::gbutton(text= L$setpermarker ,container=tabval, handler=f_setpermarker) #Drop-in
   tabval[4,3] <- gWidgets2::gbutton(text= L$setpermarker ,container=tabval, handler=f_setpermarker) #dropin lambda
   tabval[5,3] <- gWidgets2::gbutton(text= L$setpermarker ,container=tabval, handler=f_setpermarker) #fst
   
   tabval[9,1] <- gWidgets2::gbutton( L$save , container=tabval,handler = function(h, ...) { 
    opt$modeltype <- which(gWidgets2::svalue(tabval[1,2])==modtypetxt) #set index 1,2,3
    opt$threshT <- as.numeric(gWidgets2::svalue(tabval[2,2]))
    opt$dropinC <- as.numeric(gWidgets2::svalue(tabval[3,2])) 
    opt$dropinL <- as.numeric(gWidgets2::svalue(tabval[4,2]))
    opt$fst <- as.numeric(gWidgets2::svalue(tabval[5,2]))
    opt$degrad <- which(gWidgets2::svalue(tabval[6,2])==radiotxt)  #set index 1,2
    opt$stuttBW <-  which(gWidgets2::svalue(tabval[7,2])==radiotxt)  #set index 1,2 
    opt$stuttFW <-  which(gWidgets2::svalue(tabval[8,2])==radiotxt)  #set index 1,2 
  
    assign("setupModel",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileModel)    #save to file in installation folder
    
    #STORE setting INFO IN separate modelSettings Object (not used)
    assign("modelSettings",NULL,envir=nnTK)  #includes following objects: popFreq,kit,fst,lambda,threshT,xiBW,xiFW,pC (updated when saved)
    
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
  }
  
  #Function for selecting kit
  f_kitsel=  function(h,...) { #GUI function to select kit
   items0 = c(L$none,euroformix::getKit()) #no kits also possible
   opt <- get("setupKit",envir=nnTK) 
   setwin <- gWidgets2::gwindow( paste( L$select ,L$kit) ,visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text= paste0( L$selected ," ", L$kit,colonsymbol) ,container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$kitname,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text= paste( L$select , L$kit),container=tabval)
   tabval[2,2] <- gWidgets2::gcombobox(items=items0, width=100, selected=0, editable = FALSE, container = tabval)
   tabval[3,1] <- gWidgets2::gbutton( L$save , container=tabval,handler = function(h, ...) { 
    opt$kitname <-gWidgets2::svalue(tabval[2,2]) #get selected
    assign("setupKit",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileKit)    #save to file in installation folder
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
  }
  
  #GUI for importing and modifying population frequency settings
  f_popsel=  function(h,...) {
   opt <- get("setupPop",envir=nnTK) 
   popfn <-  opt$popfile #population file name
   rareOpt = get("setupRare",envir=nnTK) 
   
   #Helpfunction for storing rare allele settings from GUI
   storeRareSettings = function(dispose=TRUE) {
     rareOpt$normalize = gWidgets2::svalue(grid2[1,1])
     rareOpt$minFreq = gWidgets2::svalue(grid2[3,1])
     assign("setupRare",rareOpt,envir=nnTK)   #store to environment
     setupWrite(unlist(rareOpt),file=setupFileRare)    #save to file in installation folder if successful
     if(dispose) gWidgets2::dispose(setwin) #remove subwindow
   }
   
   setwin <- gWidgets2::gwindow( paste( L$select ,L$popfreq) ,visible=FALSE)
   grp = gWidgets2::ggroup(spacing=4,container=(setwin),horizontal = FALSE)
   grid1 = gWidgets2::glayout(spacing=0,container=grp) 
   grid1[1,1] <- gWidgets2::glabel(text= paste0( L$selected ," ", L$popfreq,colonsymbol),container=grid1)
   grid1[1,2] <- gWidgets2::glabel(text= basename( popfn ) ,container=grid1)
   gWidgets2::tooltip(grid1[1,2]) <- popfn #highligh full file name when hovered
   grid1[2,1] <- gWidgets2::gcheckbox(text=paste( L$include , L$AMEL ),checked = ifelse(opt$amel=="TRUE",TRUE,FALSE),container=grid1)
   gWidgets2::tooltip(grid1[2,1]) <- "Including AMEL into the probabilistic model (ad-hoc approach). Needs to be ticked before selecting file." 
   grid1[3,1] <- gWidgets2::gbutton( paste( L$select , L$popfreq) , container = grid1,handler = function(h, ...) { 
    ff <- mygfile(paste( L$select , L$file ),type="open")
    if(length(ff)==0) return()
    opt$popfile <- ff
    opt$amel <-gWidgets2::svalue(grid1[2,1])
    assign("setupPop",opt,envir=nnTK)  #assign user-value to opt-list
    ok <- setPopFreq(change=TRUE)	#Assume that another freq file has been selected
    if(ok) {
      setupWrite(unlist(opt),file=setupFilePop)    #save to file in installation folder if successful
      storeRareSettings(TRUE) #store rare settings first
      #gWidgets2::dispose(setwin) #remove subwindow
      #f_popsel(); #update gui window again after selecting new folder
    }
   })
   
    #BUTTON TO SHOW FREQUENCIES (OWN WINDOW)
    grid1[4,1] <- gWidgets2::gbutton( paste( L$show , L$popfreq ), container = grid1,handler = function(h, ...) { #SHOW FREQS
      ok = setPopFreq() #import population frequency from last selected file
      if(!ok) return()
      popL = get("popFreq",envir=nnTK)
      maxA <- max(sapply(popL,length)) #maximum observed
      locs <- names(popL)
      nL <- length(locs)
      nB <- 4 #block size
      tab <- matrix("",nrow=nL*nB,ncol=maxA)
      for(ll in 1:nL) { #for each loci 
        loc <- locs[ll]
        inds <- 1:length(popL[[loc]])
        tab[(ll-1)*nB + 1,1] <- loc 
        tab[(ll-1)*nB + 2,inds] <- names(popL[[loc]])
        tab[(ll-1)*nB + 3,inds] <- popL[[loc]] #frequencies
      }
      colnames(tab) <- paste0(1:ncol(tab)) #set index (column numbers)
      showGDFtable( L$popfreq, tab)
    })
    
    grid1[5,1] <- gWidgets2::gbutton( L$removeselected, container = grid1,handler = function(h, ...) { 
      opt = list(popfile="",amel="")
      assign("setupPop",opt,envir=nnTK)  
      setupWrite(unlist(opt),file=setupFilePop) #make empty
      assign("popFreq",NULL,envir=nnTK) #assign popFreq to nnTK-environment
      gWidgets2::dispose(setwin) #remove subwindow
      f_popsel() #open again
    })
    
    #Other options (about freq for rare alles)
    grid2 = gWidgets2::glayout(spacing=0,container=gWidgets2::gframe( L$rarealleles ,container=grp)) 
    grid2[1,1] <- gWidgets2::gcheckbox( L$normalizeimpute ,checked = as.logical(rareOpt$normalize),container=grid2)
    gWidgets2::tooltip(grid2[1,1]) <- "Normalize allele frequencies after inserting rare alleles (not in allele freq population file)"
    grid2[2,1] <- gWidgets2::glabel(text= paste0( L$minFreq ,colonsymbol) ,container=grid2)
    gWidgets2::tooltip(grid2[2,1]) <- "Select minimum allele frequency value to be used for rare alleles (leaving it empty will make it use minimium observed)"
    minFreq = rareOpt$minFreq
    if(is.na(minFreq)) minFreq = ""
    grid2[3,1] <- gWidgets2::gedit( minFreq,container=grid2)
    grid2[4,1] <- gWidgets2::gbutton(L$save,container=grid2, handler=function(h) storeRareSettings())
    gWidgets2::tooltip(grid2[4,1]) <- "Selecting allele frequency file will also save the settings"
    
    gWidgets2::visible(setwin) <- TRUE
    gWidgets2::focus(setwin) <- TRUE
  }#end popfreq function
  
  #Function for selecting importData function file (R-script)
  f_importsel=  function(h,...) {
   opt <- get("setupImport",envir=nnTK) 
   setwin <- gWidgets2::gwindow( paste( L$select , L$importfun) ,visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text= paste0( L$selected ," ", L$importfun,colonsymbol) ,container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$importfile,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text= paste( L$select , L$importfun),container=tabval)
   tabval[2,2] <- gWidgets2::gbutton(  L$select , container = tabval,handler = function(h, ...) { 
     ff <- mygfile(paste( L$select , L$file ),type="open")
     if(length(ff)==0) return()
     opt$importfile <- ff
     assign("setupImport",opt,envir=nnTK)  #assign user-value to opt-list
     setupWrite(unlist(opt),file=setupFileImport)    #save to file in installation folder
     gWidgets2::dispose(setwin) #remove subwindow
     #f_importsel(); #update gui window again after selecting new folder
   })
   gWidgets2::visible(setwin) <- TRUE
  }
  
  #Function for selecting case directory
  f_casedirsel =  function(h,...) {
    opt <- get("setupCase",envir=nnTK) 
    setwin <- gWidgets2::gwindow( paste( L$select , L$pathcasefolders ) ,visible=FALSE)
    tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
    tabval[1,1] <- gWidgets2::glabel(text= paste0( L$selected ," ", L$pathcasefolders,colonsymbol),container=tabval)
    tabval[1,2] <- gWidgets2::glabel(text=opt$casepath,container=tabval)
    tabval[2,1] <- gWidgets2::glabel(text= paste( L$select , L$pathcasefolders) ,container=tabval)
    tabval[2,2] <- gWidgets2::gbutton( L$select , container = tabval,handler = function(h, ...) { 
      opt$casepath <- mygfile( paste( L$select , L$folder ),type="selectdir")
      
      assign("setupCase",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileCase)    #save to file in installation folder
      gWidgets2::dispose(setwin) #remove subwindow
      gWidgets2::dispose(mainwin) #remove main window
      gui(); #update restarted gui window again after selecting new folder 
    })
    gWidgets2::visible(setwin) <- TRUE
  }
  
    #The user can change advanced model settings (nDone,maxContributors)
    f_advancedoptions = function(h,...) { 
      opt <- get("setupAdvanced",envir=nnTK) 
      setwin <- gWidgets2::gwindow( paste( L$advanced , L$options ) ,visible=FALSE)
      tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
      tabval[1,1] <- gWidgets2::glabel(text= L$maxcontrqual ,container=tabval) #"Maximum contributors in QualLR (LRmix)"
      tabval[1,2] <- gWidgets2::gedit(text=opt$maxC1,width=w0,container=tabval) 
      tabval[2,1] <- gWidgets2::glabel(text= L$maxcontrquan ,container=tabval) #Maximum contributors in QuanLR (EFM)
      tabval[2,2] <- gWidgets2::gedit(text=opt$maxC2,width=w0,container=tabval)
      tabval[3,1] <- gWidgets2::glabel(text= L$numoptim ,container=tabval) #"Number of required optimizes (EFM)"
      tabval[3,2] <- gWidgets2::gedit(text=opt$nDone,width=w0,container=tabval)
      tabval[4,1] <- gWidgets2::glabel(text= L$useoneaslowestnoc ,container=tabval) #Start with one contributor when estimating contrs
      tabval[4,2] <- gWidgets2::gcheckbox(text="",checked=opt$useMinK1=="TRUE",container=tabval)
      tabval[5,1] <- gWidgets2::glabel(text= L$showssinmatchlist ,container=tabval) #"Compare single sources"
      tabval[5,2] <- gWidgets2::gcheckbox(text="",checked=opt$compSS=="TRUE",container=tabval)
      tabval[6,1] <- gWidgets2::glabel(text= L$useSNPmodule ,container=tabval) #Use SNP module
      tabval[6,2] <- gWidgets2::gcheckbox(text="",checked=opt$isSNP=="TRUE",container=tabval)
      tabval[7,1] <- gWidgets2::glabel(text= L$profileselector ,container=tabval) #User can select profile when import/report
      tabval[7,2] <- gWidgets2::gcheckbox(text="",checked=opt$selProfiles=="TRUE",container=tabval)
      
      tabval[8,1] <- gWidgets2::gbutton( L$save , container=tabval,handler = function(h, ...) { 
        opt2 = list() #avoid wrong order 
        opt2$maxC1 <- as.numeric(gWidgets2::svalue(tabval[1,2]))  #max number of contributors in LRmix model
        opt2$maxC2 <- as.numeric(gWidgets2::svalue(tabval[2,2]))  #max number of contributors in EFM model
        opt2$nDone <- as.numeric(gWidgets2::svalue(tabval[3,2]))  #iterations in optimizer
        opt2$useMinK1 <- as.character(gWidgets2::svalue(tabval[4,2])==TRUE)  #should K=1 be used as start contr?
        opt2$compSS <- as.character(gWidgets2::svalue(tabval[5,2])==TRUE)  #Should single source profiles be compared?
        opt2$isSNP <- as.character(gWidgets2::svalue(tabval[6,2])==TRUE)  #Should SNP module be used (all evid samples are mixtures)
        opt2$selProfiles  <- as.character(gWidgets2::svalue(tabval[7,2])==TRUE)   #User can select profile when import/report
        assign("setupAdvanced",opt2,envir=nnTK)  #assign user-value to opt-list
        setupWrite(unlist(opt2),file=setupFileAdvanced)    #save to file in installation folder
        gWidgets2::dispose(setwin)
      } )
      gWidgets2::visible(setwin) <- TRUE
    }
  
  #The user can change MCMC settings (niter,delta,seed)
  f_mcmcoptions = function(h,...) { 
    opt <- get("setupMCMC",envir=nnTK) 
    setwin <- gWidgets2::gwindow( paste( L$mcmc , L$options ) ,visible=FALSE)
    tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
    
    #Traversing each option:
    for(elem in names(opt)) {
      rowind = which(names(opt)==elem) #get rowindex
      tabval[rowind,1] <- gWidgets2::glabel(text= L[[elem]] ,container=tabval) #"Maximum contributors in QualLR (LRmix)"
      tabval[rowind,2] <- gWidgets2::gedit(text=opt[[elem]],width=w0,container=tabval) 
    }
  
    #Last button is save (storing and writing settings)
    tabval[length(opt)+1,1] <- gWidgets2::gbutton( L$save , container=tabval,handler = function(h, ...) { 
      opt2 = list() #avoid wrong order 
      for(elem in names(opt)) {
        rowind = which(names(opt)==elem) #get rowindex
        opt2[[elem]] <- as.numeric(gWidgets2::svalue(tabval[rowind,2]))  #max number of contributors in LRmix model
      }
      assign("setupMCMC",opt2,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt2),file=setupFileMCMC)    #save to file in installation folder
      gWidgets2::dispose(setwin)
    } )
    gWidgets2::visible(setwin) <- TRUE
  }
  
  #The user can select language:
  f_selLanguage =  function(h,...) { #GUI function to select language
    items0 = casesolver::getLanguage() #get list of available languages
    opt <- get("setupLanguage",envir=nnTK) 
    #if(is.null(opt)) opt = list(language=defaultLanguage,encoding=defaultEncoding) #set default if not found
    setwin <- gWidgets2::gwindow( paste( L$select , L$language ),visible=FALSE)
    tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
    tabval[1,1] <- gWidgets2::glabel(text= paste0( L$selected ," ", L$language,colonsymbol) ,container=tabval)
    tabval[1,2] <- gWidgets2::glabel(text=opt$language,container=tabval)
    tabval[2,1] <- gWidgets2::glabel(text= paste0( L$select ," ", L$language,colonsymbol) ,container=tabval)
    tabval[2,2] <- gWidgets2::gcombobox(items=items0, width=100, selected=0, editable = FALSE, container = tabval)
    tabval[3,1] <- gWidgets2::gbutton( L$save , container=tabval,handler = function(h, ...) { 
      opt$language <- gWidgets2::svalue(tabval[2,2]) #get selected
      #opt$encoding <- defaultEncoding #set default encoding (could extract encoding for current language)
      assign("setupLanguage",opt,envir=nnTK)  #assign user-value to opt-list
      setupWrite(unlist(opt),file=setupFileLanguage)    #save to file in installation folder
      gWidgets2::dispose(setwin) #close windown
      
      #RESTARTING PROGRAM
      gWidgets2::dispose(mainwin) #shut down main window
      gui(envir=nnTK); #start an with same session (recognized)
    } )
    gWidgets2::visible(setwin) <- TRUE
  }
  
  
  ################################################# 
  ##########FUNCTIONALITIES HELPFUNCTIONS########## 
  #################################################
  
  #Helpfunction for resaving to setupSorting object
  resave_Sorting = function(index,value) {
    sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
    sortTypes[index] = value
    assign("setupSorting",sortTypes,envir=nnTK)
  }
  
  #MAIN FUNCTION TO IMPORT DATA
  f_importData = function(h,...) { #wrapper function which calls other functions: importData and getStructuredData
    #REMOVE PREV. RESULTS WHEN NEW IMPORT:
    if(!is.null(get("DClist",envir=nnTK)) || !is.null(get("resCompMAC",envir=nnTK))) {
      gWidgets2::gmessage( L$msg.restartfirst ) #message that user should restart before proceed with new case
    }
    caseID = gWidgets2::svalue( tabimportA[1,2] ) #get ID from table
    assign("caseID",caseID,envir=nnTK) #store ID in environment
    fn <- list.files(path=paste0(casedir,.sep,caseID), pattern="",full.names=TRUE) #get full names 
    
    print("-----------------------------")
    print("-------IMPORTING DATA--------")
    #NB: "importData" is a tailored function which must be created such that it provides a list with Evidence and Reference table:
    #Evidences: colnames(data$mix)=("SampleName","Marker","Alleles","Heights"). Separated with "/" in text format.
    #References: colnames(data$ref)=("SampleName","Marker","Alleles"). Separated with "/" in text format.
    
    #obtain importData from R-script and use to load data
    Rfile = get("setupImport",envir=nnTK)$importfile
    if(is.na(Rfile) ||  Rfile=="") {
      print("Couldn't find importData file. Returning..")
      return()
    }
    #This approach enables the R-function to have other name than importData
    importData <- source(Rfile, new.env(parent = globalenv() ),echo=FALSE)$value #obtain function with source 
    if(is.null(importData)) {
      gWidgets2::gmessage( paste( L$msg.importfunerror, L$settings,">",L$set,L$importfun ) )
      return()
    }
    data <- list(mix= matrix(nrow=0,ncol=4),ref= matrix(nrow=0,ncol=3))  #Formar: SampleName,Marker,Allele,(Height)
    markers <- numeric()
    metalist <- list() #contains list of table-elements
    consdata <- NULL #default value
  
   for(ff in fn) { #for each files: 
  #ff=fn[2]
     #if( file.info(ff)$isdir ) next #skip if it was a folder
     
     #IMPORTING DATA FROM USER-SPECIFIED FUNCTION (MUST BE NAMED "importData")
    tryCatch({ 
     data2 <- importData(ff) #import data for selected case. Structure of markers must be given inside this function and returned by "markers".
     data$mix <-  rbind(data$mix,data2$mix) #add data to table
     data$ref <-  rbind(data$ref,data2$ref) #add data to table
     markers <- data2$markers #get marker order from costumized importData file
     consdata <- rbind(consdata,data2$cons) #add data to table (consensus data)
  
     #Add metadata (assumed to be matrix/dataframes:
     tmplist <- data2$meta
     if(length(tmplist)==0) next #skip if no elements
     if(length(metalist)==0) { #if no list elements
       metalist <- tmplist 
     } else { #list elements found
       if(all(names(metalist)==names(tmplist))) { #check if containing same elements
         for(elem in names(metalist)) {
          if(length(tmplist[[elem]])>0)  { #check if not empty
            if( is.matrix(tmplist[[elem]]) ) { #BLOCK MODIFIED (v1.8.1)
              metalist[[elem]] <- rbind(metalist[[elem]],tmplist[[elem]]) #add to matrix
            } else {
              metalist[[elem]] <- c(metalist[[elem]],tmplist[[elem]]) #append to vector
            }
          }
         }
       } else {
        print("Metadata did not contain same list elements")        
       }
     }
    }, error = function(e) e) 
   }
   data$ref <- unique(data$ref) #consider only uniques
   data$mix <- unique(data$mix) #consider only uniques
   print("------------------------------")
   
   #Store non-profile data (CAN BE SHOWN IN REPORT)
   assign("metaDataLIST",metalist,envir=nnTK) #assign to nnTK-environment
   assign("consDataTABLE",consdata,envir=nnTK) #list of imported consensus data
   
   if( get("setupAdvanced",envir=nnTK)$selProfiles=="TRUE" ) {
     print("-------USER SELECTION-------") #the user may select a subset of samples to import
     guienv = new.env( parent = emptyenv() ) #create new envornment object. Parent must be empty
     assign("selected",list(evids=unique(data$mix[,1]),refs=unique(data$ref[,1])),envir=guienv)
     profileSelectorGUI(env=guienv) #calling function to select data
     selList =  get("selected",envir=guienv) #get list of 
     data = list( mix=data$mix[data$mix[,1]%in%selList$evids,,drop=FALSE] , ref=data$ref[data$ref[,1]%in%selList$refs,,drop=FALSE] ) #update data object
   }
   
   #In case of no data (RETURN WITH MESSAGE):
   if(length(markers)==0 || (nrow(data$mix)==0 && nrow(data$ref)==0) ) {
     gWidgets2::gmessage( L$msg.nodata ) #user need to select at least one profile to proceed
     return()
   }
     
   print("-------STRUCTURING DATA-------")
   datalist <- casesolver::getStructuredData(data,ln=toupper(markers),minLoc=get("setupThresh",envir=nnTK)$minLociSS) #get Data in both List-format and Table-format (mixDataTABLE,refDataTABLE,mixDataLIST)
   datalist$mixDataMATCHSTATUS = changeUnknownName(datalist$mixDataMATCHSTATUS) #update unknown names
   rownames(datalist$refDataTABLE) = changeUnknownName(rownames(datalist$refDataTABLE)) #update unknown names
   
   if(!is.null(datalist$mixDataMATCHSTATUS) && length(datalist$mixDataMATCHSTATUS)>0 && !is.null(get("setupAdvanced",envir=nnTK)$isSNP) && get("setupAdvanced",envir=nnTK)$isSNP=="TRUE") {
     datalist$mixDataMATCHSTATUS[1:length(datalist$mixDataMATCHSTATUS)] = L$mixture #assign as mixture
   } 
  
   ##Notice: The alleles in a loci should be ordered: Hence samples with same alleles can be detected   
   #STORE DATA (BOTH TYPES) -> EASY AVAILABLE THROUGH ENVIRONMENT
   assign("mixDataTABLE",datalist$mixDataTABLE,envir=nnTK) #assign to nnTK-environment
   assign("refDataTABLE",datalist$refDataTABLE,envir=nnTK) #assign to nnTK-environment
   assign("mixDataMATCHSTATUS",datalist$mixDataMATCHSTATUS,envir=nnTK) #assign to nnTK-environment
   assign("mixDataLIST",datalist$mixDataLIST,envir=nnTK) #assign to nnTK-environment
  # assign("refDataLIST",datalist$refDataLIST,envir=nnTK) #NOT USED ANYMORE!
  
   updateTables() #update datables (default sorting first time)
   refreshTabLIST() #update mixture-list (default sorting first time)
   setFocus() #set focus
  } #end import function
  
  
  #A do-all function for selected profiles(substitutes Export and Deconvolute):
  f_markprofs = function(h,...) { #helpfunction to operate on selected profiles (View/Export/Deconvolve/Delete)
    mixSelID <- refSelID <- NULL
    tryCatch( { mixSelID =  as.integer(gsub("#","",gWidgets2::svalue(mixTabGUI))) }, error=function(e) print("No EVIDS in table"))
    tryCatch( { refSelID =  as.integer(gsub("#","",gWidgets2::svalue(refTabGUI))) }, error=function(e) print("No REFS in table"))
    if(length(mixSelID)==0 && length(refSelID)==0) {
  	  gWidgets2::gmessage( L$msg.selectprofile ) #user need to select at least one profile to proceed
  	 	return()
    }
    mixL <- get("mixDataLIST",envir=nnTK)[mixSelID] #get selected mixtures
    refTab <- get("refDataTABLE",envir=nnTK)[refSelID,,drop=FALSE] #get reference table
    refL <- casesolver::tabToListRef(tab=refTab,setEmpty=FALSE) #FORCING 1-alleles to be homozygous (necessary for EFM)	
    rm(refTab);gc()
    evids <- names(mixL)
    refs <- names(refL)
    
    #Begin GUI
    selwin <- gWidgets2::gwindow( L$functionalities ,width=100,height=100,visible=FALSE)
    tabtmp <- gWidgets2::glayout(container=selwin)
    
    #Data selection:
    tabSel = gWidgets2::glayout(spacing=0,container=(tabtmp[1,1] <- gWidgets2::gframe( L$data.select ,container=tabtmp)))  
    tabSel[1,1] = gWidgets2::glabel(text= L$evidences ,container=tabSel)
    tabSel[1,2] = gWidgets2::glabel(text= L$references ,container=tabSel)
    tabSel[2,1] <- gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
    tabSel[2,2] <- gWidgets2::gcheckboxgroup(items=refs,container=tabSel,checked=TRUE)
    
    tabFun = gWidgets2::glayout(spacing=spc/2,container=(tabtmp[1,2] <- gWidgets2::gframe( L$functionalities ,container=tabtmp)))  
    
    getSelected = function() { #helpfunction to get selected
     selEvid <- selRef <-  NULL
     if(!is.null(evids) && length(evids)>0 ) selEvid <- gWidgets2::svalue(tabSel[2,1]) #get selected evids
     if(!is.null(refs) && length(refs)>0) selRef <- gWidgets2::svalue(tabSel[2,2]) #get selected refs
     return( list(Evid=selEvid,Ref=selRef) )
    }
    
    #FUNCTION 1: View data
    tabFun[1,1] = gWidgets2::gbutton(text= paste( L$show,L$data ),container=tabFun, handler=function(h,...) { 
     selL <- getSelected()
     if( all(sapply(selL,length)==0) ) return() #return if no selected profiles
     viewdata(mixL[selL[[1]]],refL[selL[[2]]]) 
    }) 
    
    #FUNCTION 2: Deconvolute
    tabFun[2,1] = gWidgets2::gbutton(text= L$deconvolve ,container=tabFun, handler=function(h,...) { 
     selL <- getSelected() #get selected evid/refs
     if(length(selL[[1]])==0) return()  #return if no selected evid profile
     createHypDCWindow(selL[[1]],selL[[2]],nC=length(selL[[2]])+1) #create window for specify model for DC  
    }) 
    
    #FUNCTION 3: Open data in EFM:
    tabFun[3,1] = gWidgets2::gbutton(text=paste( L$openin ,"EuroForMix"),container=tabFun, handler=function(h,...) { 
      load(setupFileExport) #load emtpy object with ESX17 kit popfreq in euroformix 
      selL <- getSelected() #get selected evid/refs
      if(length(selL[[1]])>0) mmTK$mixData <- mixL[selL[[1]]]
      if(length(selL[[2]])>0) mmTK$refData <- refL[selL[[2]]]
      
      #Obtain popfreq and selected kit: 
      mmTK$popFreq = nnTK$popFreq #obtain frequencies if any
      mmTK$selPopKitName = rep(NA,2)
      if(!is.null(nnTK$setupKit$kitname)) mmTK$selPopKitName[1] = nnTK$setupKit$kitname #insert selected kit
      if(!is.null(nnTK$setupPop$popfile)) mmTK$selPopKitName[2] = basename(nnTK$setupPop$popfile) #insert selected kit
      
      #Obtain settings from GUI:
      #Obtain global settings:
      mmTK$optSetup$thresh0 <- nnTK$setupModel$threshT
      mmTK$optSetup$pC0 <- nnTK$setupModel$dropinC
      mmTK$optSetup$lam0 <- nnTK$setupModel$dropinL
      mmTK$optSetup$fst0 <-  nnTK$setupModel$fst
      
      #Specify marker-based settings if given (load it into EFM object)
      if(!is.null(nnTK$setupMarkers)) {
        locNames = nnTK$setupMarkers[[1]] #obtain loci names
        mmTK$optMarkerSetup = list( 
          threshv = setNames(nnTK$setupMarkers[[2]],locNames),
          pCv = setNames(nnTK$setupMarkers[[3]],locNames),
          lamv = setNames(nnTK$setupMarkers[[4]],locNames),
          fstv = setNames(nnTK$setupMarkers[[5]],locNames))
      }
      save(mmTK,file=setupFileExport2,compress="xz") #save envir object to file (Rdata)
      euroformix::efm(envirfile=setupFileExport2) #run efm with saved file
    }) #end export to EFM
    
    
    #FUNCTION 4: Delete selected profiles from GUI
    tabFun[4,1] = gWidgets2::gbutton(text= paste(L$deletefrom, L$gui) ,container=tabFun, handler=function(h,...) { 
     	selL <- getSelected()
     	if( length(selL[[1]])==0 && length(selL[[2]])==0) return() #return if none selected
    	txt = L$msg.deleteprofiles #delete following profiles?
    	sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types (all tables except of matchMatrix)
    	
    	#Obtain match lists (matchmatrix,qualLR,quanLR,matchlist)
    	resCompMAC = get("resCompMAC",envir=nnTK)
    	resCompLR1 = get("resCompLR1",envir=nnTK)
    	resCompLR2 = get("resCompLR2",envir=nnTK)
    	resCompLR  = get("resCompLR",envir=nnTK) #overview of combined results
    	resMatches = get("resMatches",envir=nnTK) #this is final match results
    	allMixList = get("allMixList",envir=nnTK) #list of all evidence profiles
    	#storedFitHp = get("storedFitHp",envir=nnTK)

    	#if removing any evid profiles: 
    	anyIsDeleted = FALSE
    	if( length(selL[[1]])>0 ) {  
        bool <- gWidgets2::gconfirm(paste0(txt,"\n",paste0(selL[[1]],collapse="\n")))
        if(bool) {
          if( !is.null(get("resCompMAC",envir=nnTK)) || !is.null(get("DClist",envir=nnTK)) ) {
            gWidgets2::gmessage("You can't delete evidence profiles after comparison or deconvolution!")
          } else { #If deletion was possible
            anyIsDeleted = TRUE
            mixTab = get("mixDataTABLE",envir=nnTK) #get evid table
            mixStatus = get("mixDataMATCHSTATUS",envir=nnTK) #get status
            mixList = get("mixDataLIST",envir=nnTK) #get evid list table
            allevids = names(mixList) #rownames(tabT) #get name of all evidence
            keepevids = setdiff(allevids,selL[[1]])
            keepevidsInd = which(allevids%in%keepevids) #obtain index for keeping
            mixStatus = mixStatus[keepevidsInd]
            mixTab = mixTab[keepevidsInd,,drop=FALSE] #update table
            
            #Remove evid from comparisons
            if(FALSE) {
              if(!is.null(resCompLR1)) resCompLR1 = resCompLR1[resCompLR1[,1]%in%keepevids,,drop=FALSE]
              if(!is.null(resCompLR2)) resCompLR2 = resCompLR2[resCompLR2[,1]%in%keepevids,,drop=FALSE]
              if(!is.null(resCompLR)) resCompLR  = resCompLR[resCompLR[,1]%in%keepevids,,drop=FALSE]
              if(!is.null(resMatches)) resMatches  = resMatches[resMatches[,1]%in%keepevids,,drop=FALSE]
              if(!is.null(allMixList)) allMixList = allMixList[allMixList[,1]%in%keepevids,,drop=FALSE]
              if(!is.null(resCompMAC)) {
                resCompMAC$MatchList = resCompMAC$MatchList[resCompMAC$MatchList[,1]%in%keepevids,,drop=FALSE]
                if(any( colnames(resCompMAC$MatchMatrix)%in%selL[[1]]) ) { #if columns was evids
                  resCompMAC$MatchMatrix= resCompMAC$MatchMatrix[,colnames(resCompMAC$MatchMatrix)%in%keepevids,drop=FALSE]
                } else {  #if rows was evids
                  resCompMAC$MatchMatrix= resCompMAC$MatchMatrix[rownames(resCompMAC$MatchMatrix)%in%keepevids,,drop=FALSE]
                }
              }
              #storedFitHp #REFRESH DECONVOLUTION TAB?
            } #end
            
            #REASSIGNING VALUES
            assign("mixDataTABLE",mixTab,envir=nnTK) #store evid table 
            assign("mixDataMATCHSTATUS",mixStatus,envir=nnTK) #store match status again 
            assign("mixDataLIST", get("mixDataLIST",envir=nnTK)[rownames(mixTab)],envir=nnTK) #store evid list
            updateTables(type="mix",sort=sortTypes[1]) #updates evid tables again 
          }
        } 
    	} #end if remove evid profiles
    	
    	#if removing any ref profiles
    	if(length(selL[[2]])>0) {  
        bool <- gWidgets2::gconfirm(paste0(txt,"\n",paste0(selL[[2]],collapse="\n")))
        if(bool) {
          anyIsDeleted = TRUE
          refTab = get("refDataTABLE",envir=nnTK) #get table with reference data
          allrefs = rownames(refTab)
          keeprefs = setdiff(allrefs,selL[[2]])
          mixStatus = get("mixDataMATCHSTATUS",envir=nnTK) #get status
          refTab = refTab[allrefs%in%keeprefs,,drop=FALSE] #update reference table
          mixStatusIndremove = mixStatus%in%selL[[2]] #get bool of those to remove
          mixStatus[mixStatusIndremove] = "" #NA #set NA as empty if removed
          assign("mixDataMATCHSTATUS",mixStatus,envir=nnTK) #store match status again 
          assign("refDataTABLE",refTab,envir=nnTK) #store ref table 
          
          #Helpfunction to remove references in table
          removeRef = function(tab,refsrm) {
            rmIdx = NULL #obtain indices where reference is removed
            for(ref in refsrm) rmIdx = c(rmIdx,grep(ref, tab[,2]))
            for(i in unique(rmIdx)) {
              refs = strsplit(tab[i,2],"/")[[1]]  #obtain refs
              refs = refs[!refs%in%refsrm]
              tab[i,2] = paste0(refs,collapse="/")
            }
            return(tab)
          }
          #Update comparison tables:
          if(!is.null(resCompLR1)) resCompLR1 = resCompLR1[resCompLR1[,2]%in%keeprefs,,drop=FALSE]
          if(!is.null(resCompLR2)) resCompLR2 = resCompLR2[resCompLR2[,2]%in%keeprefs,,drop=FALSE]
          if(!is.null(resCompLR)) resCompLR  = resCompLR[resCompLR[,2]%in%keeprefs,,drop=FALSE]
          if(!is.null(allMixList)) allMixList = removeRef(allMixList,selL[[2]])
          if(!is.null(resMatches)) resMatches = removeRef(resMatches,selL[[2]])
          if(!is.null(resCompMAC)) {
            resCompMAC$MatchList = resCompMAC$MatchList[resCompMAC$MatchList[,2]%in%keeprefs,,drop=FALSE]
            if(any( rownames(resCompMAC$MatchMatrix)%in%selL[[2]])) { #if columns was evids
              resCompMAC$MatchMatrix = resCompMAC$MatchMatrix[rownames(resCompMAC$MatchMatrix)%in%keeprefs,,drop=FALSE]
            } else {  #if rows was evids
              resCompMAC$MatchMatrix= resCompMAC$MatchMatrix[,colnames(resCompMAC$MatchMatrix)%in%keeprefs,drop=FALSE]
            }
          }
       	  #Update tables     	  
       	  updateTables(type="ref",sort=sortTypes[2]) #updates tables again 
       	  if(any(mixStatusIndremove)) updateTables(type="mix",sort=sortTypes[1]) #updates mix tables again 
      	} #end bool
    	}
    	
    	if(anyIsDeleted) { #restore tables if any deleted
      	#Save match lists (matchmatrix,qualLR,quanLR,matchlist) and update
      	assign("resCompMAC",resCompMAC,envir=nnTK)
      	assign("resCompLR1",resCompLR1,envir=nnTK)
      	assign("resCompLR2",resCompLR2,envir=nnTK)
      	assign("resCompLR",resCompLR,envir=nnTK) #overview of combined results
      	assign("resMatches",resMatches,envir=nnTK) #this is final match results
      	assign("allMixList",allMixList,envir=nnTK) #list of all evidence profiles
      	#assign("storedFitHp",storedFitHp,envir=nnTK)
      	refreshTabMATRIX(1)
      	refreshTabLIST1(sort=sortTypes[3]) #update match list (quanLR)
      	refreshTabLIST2(sort=sortTypes[4]) #update match list (quanLR)
      	refreshTabLIST(sort=sortTypes[5]) #update match list (final)
    	}
    }) #end delete funciton
    
    #FUNCTION 5: Modify match s selected profiles from GUI
    tabFun[5,1] = gWidgets2::gbutton(text= paste(L$select,L$matchstatus) ,container=tabFun, handler=function(h,...) { 
      selL <- getSelected()
      sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
      #if removing any evid profiles: 
      evidSel = selL[[1]] #obtain selected evidence
      if( length(evidSel)==0 ) return()
      #if( !is.null(get("resCompMAC",envir=nnTK)) || !is.null(get("DClist",envir=nnTK)) ) {
      #gWidgets2::gmessage("You can't delete evidence profiles after comparison or deconvolution!")
      refTab = get("refDataTABLE",envir=nnTK) #get table with reference data
      refNames = rownames(refTab) #obtain reference names
      mixStatus = get("mixDataMATCHSTATUS",envir=nnTK) #get status for evids
      #Note: samples with "mixture" are separated from the others in report

      #Create dropdown meny with choices:
      items = c(L$none, L$mixture, refNames)
      itemselwin <- gWidgets2::gwindow( paste(L$select,L$matchstatus),visible=FALSE)
      itemselgrid =  gWidgets2::glayout(container=itemselwin)
      itemselgrid[1,1] = gWidgets2::glabel("Evidence(s):",container=itemselgrid) #export peak heights?
      itemselgrid[1,2] = gWidgets2::glabel( paste0(evidSel,collapse="/") ,container=itemselgrid) #export peak heights?
      itemselgrid[2,1] = gWidgets2::gcombobox(items,container=itemselgrid) #export peak heights?
      itemselgrid[2,2] = gWidgets2::gbutton(text=L$save,container=itemselgrid, handler=function(h,...){
        matchStatusSel = gWidgets2::svalue(itemselgrid[2,1]) #Obtain value from dropdown list
        if(matchStatusSel==items[1]) matchStatusSel = "" #make empty
        if(matchStatusSel==items[2]) matchStatusSel = "mixture" #insert as mixture (translated after)
        mixStatus[names(mixStatus)%in%evidSel] = matchStatusSel #modify matchstatus
        assign("mixDataMATCHSTATUS",mixStatus,envir=nnTK) #store match status again 
        updateTables(type="mix",sort=sortTypes[1]) #updates evid tables again 
        gWidgets2::dispose(itemselwin) #close selection window
      }) #Modify matchstatus?
      gWidgets2::visible(itemselwin)=TRUE
      gWidgets2::size(itemselwin) = c(300,30)
    }) 
    
    #Export reference file panels
    tabExp = gWidgets2::glayout(spacing=spc/2,container=(tabtmp[1,3] <- gWidgets2::gframe( L$fileexport ,container=tabtmp)))  
    	tabExp[2,1] = gWidgets2::gcheckbox(text=paste0( L$export ,"\n", L$withPH ),container=tabExp,checked=TRUE) #export peak heights?
    	tabExp[1,1] = gWidgets2::gbutton(text= paste0( L$export ,"\n", L$evidences ) ,container=tabExp,handler=function(h,...) {  #"Store\nevidence(s)
      selEvid <- gWidgets2::svalue(tabSel[2,1]) #get selected evids
      mixL2 = mixL[selEvid]
     	 if(length(mixL2)>0)  saveTable( sampleListToTable(mixL2,PH=gWidgets2::svalue(tabExp[2,1])) )
    })
    tabExp[1,2] = gWidgets2::gbutton(text= paste0( L$export ,"\n", L$references ) ,container=tabExp,handler=function(h,...) {  #"Store\nreference(s)"
     selRef <- gWidgets2::svalue(tabSel[2,2]) #get selected refs
     refL2 = refL[selRef]
     	 if(length(refL2)>0) saveTable(sampleListToTable(refL2))
    })
    gWidgets2::visible(selwin)=TRUE
  } #end function
  
  
  f_calcRMP = function(h,...) {  #Function to calculate RMP for each references and RMNE for all evidence
    isok = setPopFreq(giveMessage=TRUE) #import population frequency from last selected file (not required)
    if(!isok) return()
    resRMP = casesolver::calcRMP(nnTK)
    colnames(resRMP$evid) <- c(".", L$samplename , L$RMNE )
    colnames(resRMP$ref) <- c(".", L$samplename , L$RMP ) #insert column names
    assign("resRMP",resRMP,envir=nnTK)  #store random match prob results (to be shown in report)
  
    setwin <- gWidgets2::gwindow(paste0( L$randommatchprob ),visible=FALSE) 
    tabval <- gWidgets2::ggroup(container=gWidgets2::gframe( paste0( L$evidences,paste0(rep("\t",4),collapse=""), L$references ),container=setwin,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
    guitab1 <- gWidgets2::gtable(items=resRMP$evid,container = tabval) 
    guitab2 <- gWidgets2::gtable(items=resRMP$ref,container = tabval) 
    gWidgets2::add(tabval,guitab1,expand=TRUE,fill=TRUE) 
    gWidgets2::add(tabval,guitab2,expand=TRUE,fill=TRUE) 
    gWidgets2::visible(setwin) <- TRUE 
  } #end function
  
  #Function to calculate concordance between evidences: Much of the code is similar as IBS
  f_calcEvidConc = function(h,...) {  
    setPopFreq(giveMessage=FALSE) #import population frequency from last selected file (not required)
     
    out = NULL
    tryCatch({
      out = casesolver::calcEvidConc(nnTK)
    },error = function(e)  
       gWidgets2::gmessage( L$msg.largesampleerror ) #throw large sample eror  
    )
    if(!is.null(out)) {
     colnames(out) <- c( L$concordance , L$comparison, "." )# "Concordance","Comparison")
     assign("resEvidConc",out,envir=nnTK)  #store random match prob results (to be shown in report)
    
     #show candidates in GUI:
     out = cbind(1:nrow(out),out) #include rank
     colnames(out)[1] = L$rank 
     setwin <- gWidgets2::gwindow(paste0( L$concordantevidences ) ,visible=FALSE) 
     tabval <- gWidgets2::ggroup(container=gWidgets2::gframe( L$comparisons ,container=setwin,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE,horizontal = FALSE) #evidence,ref dataframe
     guitab1 <- gWidgets2::gtable(items=out,container = tabval) 
     gWidgets2::add(tabval,guitab1,expand=TRUE,fill=TRUE) 
     gWidgets2::visible(setwin) <- TRUE 
    } else {
     gWidgets2::gmessage( L$msg.nocandidates )
     assign("resEvidConc",matrix(nrow=0,ncol=3),envir=nnTK)  #store empty table 
    }
  }
  
  
  f_calcIBS = function(h,...) {  #Function to calculate IBS between references
    setPopFreq(giveMessage=FALSE) #import population frequency from last selected file (not required)
    tabIBS = casesolver::calcIBS(nnTK,nLarge,L$mixture) #calculate IBS
    
    if(!is.null(tabIBS)) {
     colnames(tabIBS) <- c( L$nummismatch, L$comparison, ".", L$IBS, L$nummarkers  ) #note the added column names
     assign("resIBS",tabIBS,envir=nnTK)  #store random match prob results (to be shown in report)
    
     #show candidates in GUI:
     tabIBS = cbind(1:nrow(tabIBS),tabIBS) #include rank
     colnames(tabIBS)[1] = L$rank 
     setwin <- gWidgets2::gwindow( L$concordantreferences ,visible=FALSE)  #"Similarity of references"
     tabval <- gWidgets2::ggroup(container=gWidgets2::gframe( L$comparisons ,container=setwin,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE,horizontal = FALSE) #evidence,ref dataframe
     guitab1 <- gWidgets2::gtable(items=tabIBS,container = tabval) 
     gWidgets2::add(tabval,guitab1,expand=TRUE,fill=TRUE) 
     gWidgets2::visible(setwin) <- TRUE 
    } else {
     gWidgets2::gmessage( L$msg.nocandidates )
     assign("resIBS",matrix(nrow=0,ncol=3),envir=nnTK)  #store empty IBS 
    }
  }
  
  f_calcIBSdist = function(h,...) {  #Function to calculate IBD between references
    ok = setPopFreq() #import population frequency from last selected file
    if(!ok) return() 
    
    #Ask for user permission to run (because of extensive calculation)
    userin = gWidgets2::gconfirm(  L$msg.calculationconfirmwarning ,icon="warning")
    if(userin) {
     dist <- casesolver::getIBSdistr(get("popFreq",envir=nnTK)) #simulate random match probas
     tab <- 1-cumsum(dist[,2]) #consider cumulative probs
     names(tab) <- dist[,1]
     barplot(tab,main="Random allele sharing under unrelatedness",xlab="x=number of shared alleles",ylab="Prob. sharing>=x") 
     #abline(h=0.05,lty=2)
     print(tab)
    }
  }
  
  #Helpfunction to create TopPlots (after model fit)
  makePlotTop = function(type,contFit,dc,kitname=NULL) {
   if(require(plotly)) {
     switch(type,
      "EPG" = euroformix::plotTopEPG2(contFit,dc,kit=kitname), #,AT=threshT)
      "LUS" = euroformix::plotTopMPS2(contFit,dc,grpsymbol=LUSsymbol),
      "MPS" = plotTopMPS2(contFit,dc, grpsymbol=MPSsymbol) )
     
   } else {
     switch(type,
      "EPG" = euroformix::plotTopEPG(contFit,dc,kitname=kitname),#threshT=threshT) 
      "LUS" = euroformix::plotTopLUS(contFit,dc,LUSsymbol=LUSsymbol),
      "MPS" = gWidgets2::gmessage("Please install plotly to show plot!") )
   }
  } #end function
  
  
  viewdata = function(mixL,refL=NULL,printout=TRUE) { #helpfunction to visualize data
    kitname = casesolver::getEnvirKit(nnTK)
    threshT=get("setupModel",envir=nnTK)$threshT
    type = casesolver::getSampleType2(mixL, kitname) #get sample types (EPG/MPS/LUS)
    
    #print out to console:
    if(printout) {
     for(mixid in names(mixL)) { #print out evidence profiles 
         print("--------------------------------") 
         print(mixid) 
         out = rbind(sapply(mixL[[mixid]],function(x) paste0(x$adata,collapse="/")),sapply(mixL[[mixid]],function(x) paste0(x$hdata,collapse="/")))
         rownames(out) = c("Alleles","Heights")
         print(out) 
     }  
     out = numeric()
     for(refid in names(refL)) { #print out evidence profiles 
         out = rbind( out, sapply(refL[[refid ]],function(x) paste0(x$adata,collapse="/")) )
     }      
     if(length(out)>0) {
       rownames(out) = names(refL)
       print("--------------------------------") 
       print(out) 
     }
    } #end if printout
    if( length(mixL)==0) return() #return if no evidence profiles given
    
    #Show plots in browser if plotly is installed (PRETTY)
    if( require(plotly) ) { 
      switch(type,
             "EPG" =  euroformix::plotEPG2(mixData=mixL,refData=refL,kit=kitname,AT=threshT),
             "LUS" =  euroformix::plotMPS2(mixData=mixL,refData=refL,AT=threshT,grpsymbol=LUSsymbol),
             "MPS" =  euroformix::plotMPS2(mixData=mixL,refData=refL,AT=threshT,grpsymbol=MPSsymbol))
      
    } else { #Otherwise we show in Rgui
      if(type=="MPS") {
        gWidgets2::gmessage("Please install plotly to show plot!")
      } else {
        dev.new(width=25, height=10)
        if(type=="LUS") {
          for(mixid in names(mixL)) { 	 
            euroformix::plotLUS(mixData=mixL[[mixid]],sn=mixid,refData=refL,LUSsymbol=LUSsymbol,threshT=threshT) 
            if(length(mixL)>1) dev.new() #plot in new window
          }
        } else {
          euroformix::plotEPG(mixL,refcond=refL,kitname=kitname,threshT=threshT)
        }
        dev.new()
        op <- par(no.readonly = TRUE)
        dev.off()
        par(op)
        setFocus()
      } #end if not MPS
    }
  } #end viewdata function
  
  
  #Function executed when double clicked one of the data tables
  clicktable = function(h,...) { #UPDATED V1.5 to HANDLE MPS FORMAT
  last <- get("clicktableLast",envir=nnTK) #assign to nnTK-environment
  mixL <- get("mixDataLIST",envir=nnTK)
  
  if(h$action=="mix") { #IF evid profile selected
   id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
   last <- list(mix=id) #store id, refs deleted
   assign("clicktableLast",last,envir=nnTK) #
   viewdata(mixL[id])
  }
  if(h$action=="ref" && !is.null(last$mix)) { #IF ref profile selected
   id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
   allid <- c(last$ref,id) #remember the other refs also
   last <- list(mix=last$mix,ref=allid) #store ids
   assign("clicktableLast",last,envir=nnTK) #
   refL <- casesolver::tabToListRef(tab=get("refDataTABLE",envir=nnTK)[allid,,drop=FALSE],setEmpty=FALSE) #FORCING DUP alleles
   viewdata(mixL[last$mix],refL)
  } 
  } #end function
  
  f_calcQuanLRall = function(h,...) { #function to run LR for all comparisons having LR>LRthresh(Qual)
   sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
   
   getMatchesLR(type="quan") #run EFM 
   refreshTabLIST2(sort=sortTypes[4]) #update QUAN LR table with results
  
   #Create matchlist (Final step)
   createMatchlist(modtype=2) #update matchlist with results from QUAN LR
   refreshTabLIST(sort=sortTypes[5]) #update tables  
   gWidgets2::svalue(nb) <- 5 #go to overview when done
  }
  
  #Function executed when double clicked on the matchlist (LRmix): Specify #contr and calc LR.
  clickmatchlistQUAL = function(h,...) {
    LRtab <- get("resCompLR1",envir=nnTK) #get match list (sorted)
    if(is.null(LRtab)) return()
    matchlist <- get("resMatches",envir=nnTK) #get match list in mix-table
    id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
    evid = LRtab[id,1]
    ref = LRtab[id,2]
    lrval <- format(as.numeric(LRtab[id,4]),digits=4) #get LR for ref
    matchind = which(matchlist[,1]==evid) #get index of matchlist
    condRefs = NULL
    if(length(matchind)>0) {
     condRefs= unlist(strsplit(matchlist[matchind,2],"/"))
     condRefs = setdiff(condRefs,ref)
     if(length(condRefs)==0) {
       condRefs = NULL
       lrval2 = NULL
     } else {
       LRtab0 <- get("resCompLR",envir=nnTK) #get default match list
       subtab = LRtab0[ LRtab0[,1]==evid ,,drop=FALSE] #get relevant evidence 
       lrval2 <- format( as.numeric(na.omit(subtab[match(condRefs,subtab[,2]),4])) ,digits=4) #Order fixed in v2.1.0
     }
    } else {
      lrval2 = NULL #NONE is also possible
    }
    #User may specify number of contributors and conditionals
    #OPEN WINDOW TO SPECIFY HYP:
    suppressWarnings({
    #evids=evid;nC=as.integer(LRtab[id,5]);lrvals=lrval;condRefs=condRefs
     createHypLRWindow(evids=evid,ref=ref,condRefs=condRefs,nC=as.integer(LRtab[id,5]),lrvals=lrval,lrvals2=lrval2)  #creates a window to specify hypotheses
    })
  }
  
  #Function executed when double clicked on the matchlist (EFM): Show expected peak heights
  clickmatchlistQUAN = function(h,...) {
    id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
    tab <- get("resCompLR",envir=nnTK) #get match list (unsorted but truncated)
    tab2 <- get("resCompLR2",envir=nnTK) #get match list (sorted)
    ind <- which(tab2[id,1]==tab[,1] & tab2[id,2]==tab[,2]) #get correct index of unsorted list
    if(length(ind)==0) {
      tab1 <- get("resCompLR1",envir=nnTK) #get match list (sorted and not truncated)
      ind <- which(tab2[id,1]==tab1[,1] & tab2[id,2]==tab1[,2]) #get correct index of unsorted list
    }
    contFit <- get("storedFitHp",envir=nnTK)[[ind]]  #get stored model fit under Hp
    kitname = getEnvirKit(nnTK) #get selected kit
    
    suppressWarnings({
     dc <- euroformix::deconvolve(contFit,maxlist=1) #get top candidate profiles
    })
    
    #Again: determine whether it is EPG/MPS(LUS)/(strings): Check if "_" is used. Otherwise check with all alleles are strings
    type = casesolver::getSampleType2(contFit$model$samples,kitname) #get sample type
    makePlotTop(type,contFit,dc,kitname) #create plot
  } #end if clickmatchlistQUAN 
  
  #Helpfunction to get refList for specific given refs (single alleles are always filled twice, assumed homozygouz)
  getRefL = function(refs) { #return list with same order as for refs
    reftab=get("refDataTABLE",envir=nnTK)
    refL <- casesolver::tabToListRef(tab=reftab[match(refs,rownames(reftab)),,drop=FALSE],setEmpty=TRUE) #Genotypes with 1-alleles are ignored
    rm(reftab);gc()
    return(refL)
  }
  
  #FUNCTION WHICH PERFORMS DC (uses settings in GUI)
  doDC = function(nC,evids=NULL,refs=NULL,showPlot=TRUE,useplotly=TRUE,addedProfiles=NULL) {
   useplotly <- useplotly && require(plotly) #must be installed
   
   refData <- condOrder <- NULL
   if(!is.null(refs) && length(refs)>0) {
     refData = getRefL(refs) # get list of reference data
     condOrder = 1:length(refs) #hypothesis is to condition on all references
   }
   evidData <- get("mixDataLIST",envir=nnTK)[evids] #evidence to consider
   evids = paste0(evids,collapse="/") #collapse multiple evidence names
   condrefs = paste0(refs,collapse="/")
   
   suppressWarnings({ 
     contFit <- casesolver::calcQuanMLE(evidData,refData,condOrder,nC,nnTK,verbose=TRUE) #get fitted object
     dc <- euroformix::deconvolve(contFit,maxlist=1) #get top candidate profiles
   })
   if(showPlot) {
     kitname = casesolver::getEnvirKit(nnTK) #get kitname
     type = casesolver::getSampleType2(evidData,kitname) #get sample type
  
     #INSERTING ADDED PROFILE BY MANIPULATE FITTED MLE-model
     if(!is.null(addedProfiles)) {
       locs = names(dc$toprankGi) #obtain locus names
       for(loc in locs) {
         if(is.null( contFit$model$refData)) contFit$model$refData = list() #must create list if not exist
         for(comp in names(addedProfiles)) {
           if(is.null( contFit$model$refData[[loc]])) contFit$model$refData[[loc]] = list() #must create list if not exist
           alleles = addedProfiles[[comp]][[loc]]$adata
           contFit$model$refData[[loc]][[comp]] = alleles #insert aleles
         }
       }
       #conditional index to insert
       if(!is.null(contFit$model$condOrder)) {
         contFit$model$condOrder = c(contFit$model$condOrder, max(contFit$model$condOrder)+1)
       } else {
         contFit$model$condOrder = 1
       }
     }
     
     tryCatch({
       makePlotTop(type,contFit,dc,kitname)
      }, error = function(e) print(e))
   } #end if showPlot
   if(!is.null(addedProfiles)) return() #stopfunction if there was added profiles to show in plot (special case)
   
   #INSERTING CANDIDATE DECONVOLED PROFILES:
   ratio <- get("setupThresh",envir=nnTK)$ratio #get ratio-threshold
   probA <- get("setupThresh",envir=nnTK)$probA #get probability of allele - threshold
   
   locs0 <- names(dc$toprankGi) #get loci (from DC)
   contrs = colnames(dc$toprankGi[[1]]) #get contributors
   candtab <- matrix(nrow=0,ncol=2*length(locs0)+4) #list of candidates with genotypes (and probabilities)
   colnames(candtab) <- c(L$Component, L$Conditionals,L$NOC,L$MixProp,locs0,locs0)
   nR = length(refData) #number of conditional refs
   for(cind in 1:length(contrs)) { #for each contributors
     compn <-  contrs[cind] #paste0(evid,"_C",uind) #component name
     addRef = FALSE #Should the profile be added? (Must have at least 1 deduced allele for an unknown component)
     mxhat <- contFit$fit$thetahat2[cind] #get mixture proportion
     newrow <- rep(NA,2*length(locs0))
     for(loc in locs0) { #for each locus 
       insind <- which(locs0==loc) #insert index for genotypes
       insind2 <- length(locs0) + insind #insert index for probabilities (added last)
       if(is.null(dc$toprankGi[[loc]])) next #skip if marker not found
       cand <- dc$toprankGi[[loc]][,cind] #get candidate  
       candRatio = as.numeric(cand[3]) #obtain 'ratio to next' for candidate
       if(!is.na(cand[3]) &&  candRatio<ratio) { #if not a likely genotype
         ind <- which(dc$table4[,1]==compn & dc$table4[,2]==loc)[1] #find top ranked single allele
         candProbA = as.numeric(dc$table4[ind,4]) #get allele prob for candidate
         if( candProbA >= probA ) {
           newrow[insind] <- dc$table4[ind,3] #insert allele if prob>probA
           #newrow[insind2] <- candProbA #insert allele prob
         } 
       } else { #else insert genotype candidate
         newrow[insind] <- cand[1] #insert genotype
       }
       newrow[insind2] <- signif(candRatio,2) #insert ratio, rounded, to marker (always)
       if( cind <= nR && length(refData[[cind]][[loc]]$adata)==0 && !is.na(newrow[insind]) ) addRef  = TRUE #indicate that ref prof. should be added
     } #end for each loci
     if(all(is.na(newrow[1:length(locs0)]))) next #skip if not deduced genotypes (Notice change from v1.8 when adding probabilities)
     if( cind <= nR && !addRef) next #skip if ref should not be added
     if(cind <= nR ) compn = names(refData)[cind] #use ref name instead if conditioned on
     newrow <- c(paste0(evids,"-",compn),condrefs,nC,signif(mxhat,2),newrow)
     candtab <- rbind(candtab, newrow)
   } 
   return(candtab)
  } #end doDC
  
  #helpfunction to change "Unknown" to language specific name
  changeUnknownName = function(x) {
   unknowntxt = "Unknown" #this is default name
   isUnknown=substr(x,1,nchar(unknowntxt))==unknowntxt #recognize positions
   x = gsub(unknowntxt,L$unknown,x) #update vector
   return(x)
  }
  
  #Helpfunction to calc LR (fit quan model for both hyps)
  #REPLICATES NOT SUPPORTED!
  fitEFMHYPs = function(nC,evids,ref,condref=NULL) {
    #Prepare data:
    #popFreq=get("popFreq",envir=nnTK);
    samples <- get("mixDataLIST",envir=nnTK)[evids] #consider lists
    #uselocs <-  names(popFreq)[toupper(names(popFreq))%in%toupper(names(samples[[1]]))] #loci to consider
    #if(is.null(xi)) uselocs <- setdiff(uselocs,"AMEL") #EXCLUDE AMEL IF STUTTER CONSIDERD (QUICK SOLUTION) 
    #samples = lapply(samples,function(x) x[uselocs])
    allrefs = c(ref,condref) #get all refs. POI is first ref
    refData = getRefL(allrefs) # get list of reference data (1-alleles are put as empty)
  
    condhp <- 1:length(allrefs) #conditional order must be increasing order
    condhd <- condhp - 1 #don't condition on POI under Hd
  
    #RUN CALCULATIONS:
    fitMLE = function(condOrder) {
      mlefit <- casesolver::calcQuanMLE(samples,refData,condOrder,nC,nnTK,verbose=TRUE)
      return(mlefit)
    }
    #Calculate under Hp:
    fithp = fitMLE(condhp)
    print("Estimate under Hp:")
    print(prettyNum(fithp$fit$thetahat2))
  
    #Calculate under Hd:
    fithd = fitMLE(condhd)
    print("Estimate under Hd:")
    print(prettyNum(fithd$fit$thetahat2))
    mleLR = (fithp$fit$loglik - fithd$fit$loglik)/log(10) #get estimated log10LR (Quan based)
  
    return(list(fithp=fithp,fithd=fithd,mleLR=mleLR))
  } #end function
  
  #FUNCTION TO SPECIFY HYPOTHESIS AND CALCUALTE SINGLE LR
  createHypLRWindow = function(evids,ref,nC,condRefs=NULL,lrvals=NULL,lrvals2=NULL) {
    Krange = 1:4 #range of number of contributors
    selind = which(nC==Krange)
    if(length(selind)==0) selind = max(Krange)
    selwin <- gWidgets2::gwindow(paste0( L$quanmodel ), visible=TRUE)
    tabSel <- gWidgets2::glayout(container=selwin)
    tabSel[1,1] = gWidgets2::glabel(text= paste0( L$evidence ,colonsymbol) ,container=tabSel)
    tabSel[1,2] = gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
    tabSel[2,1] = gWidgets2::glabel(text= paste0( L$POI ,colonsymbol),container=tabSel)
    tabSel[2,2] = gWidgets2::glabel(text=ref,container=tabSel)
    tabSel[3,1] = gWidgets2::glabel(text= paste0( L$numcontr , colonsymbol) ,container=tabSel) 
    tabSel[3,2] = gWidgets2::gcombobox(items=Krange,selected=selind,editable=TRUE,container=tabSel) 
    # gWidgets2::gedit(text=nC,container=tabSel)
    tabSel[4,1] = gWidgets2::glabel(text= paste0( L$condtionON, colonsymbol) ,container=tabSel)
    tabSel[4,2] = gWidgets2::gcheckboxgroup(items=condRefs,container=tabSel,checked=FALSE) #don't select by default
    if(!is.null(lrvals)) { #if lr values given
     tabSel[1,3] = gWidgets2::glabel(text= L$logLR ,container=tabSel)
     tabSel[2,3] = gWidgets2::glabel(lrvals,container=tabSel)
    }
    if(!is.null(lrvals2)) { #if lr values given
     tabSel[4,3] = gWidgets2::glabel(lrvals2,container=tabSel)
    }
    tabSel[5,1] = gWidgets2::gbutton(text= L$calculate ,container=tabSel,handler=function(h,...) {
      nC2 <- as.integer(gWidgets2::svalue(tabSel[3,2])) #get selected number of contributors
      if(nC2>get("setupAdvanced",envir=nnTK)$maxC2) gWidgets2::gmessage( paste( L$msg.morecontrwarning, L$advanced , L$options  ) ,title= L$warning )
    
      condsel <- NULL #conditional references
      if(!is.null(condRefs) && length(condRefs)>0) {
       condsel <- gWidgets2::svalue(tabSel[4,2]) #get selected ref to condition on
      }
      mixsel <-  gWidgets2::svalue(tabSel[1,2]) #get selected samples
      if(length(mixsel)==0) {
        gWidgets2::gmessage( L$msg.selectprofile )
        return()
      }
      gWidgets2::dispose(selwin) #remove window 
      ret = fitEFMHYPs(nC=nC2,evids=mixsel,ref=ref,condref=condsel) 
    
      #user can choose whether to replace results
      txt = paste0( L$msg.calculatedLR ,"\n", L$logLR ,"=",signif(ret$mleLR,digits=4),"\n\n", L$msg.useresult )
      ubool <- gWidgets2::gconfirm(txt,title= L$msg.useresult ,icon="info")
      if(!ubool) return() #return from function 
    
      #STORE MATCH HERE (not in function directly)  
      matchlist0 = get("resCompLR",envir=nnTK)  #get stored results from Qual/Quan comparison (unsorted, but truncated)
      matchlist1 = get("resCompLR1",envir=nnTK)  #get stored qual based results from comparison (sorted wrt QualLR)
      matchlist2 = get("resCompLR2",envir=nnTK)  #get stored quan based results from comparison (sorted wrt QuanLR)
      hpfitlist = get("storedFitHp",envir=nnTK)  #get stored model results under Hp
    
      condREF = FALSE #bool if cond. ref is included
      hasCondREF = FALSE #whether condtional ref columns is in table
      if(!is.null(condsel) && length(condsel)>0) condREF = TRUE #add cond.refs to table?
      if(!is.null(matchlist2) && L$condref%in%colnames(matchlist2)) hasCondREF = TRUE 
      if(!is.null(matchlist2) && condREF && !hasCondREF) { #need to extend the table
        matchlist2 = cbind(matchlist2,"") #add an extra column
        colnames(matchlist2)[ncol(matchlist2)] = L$condref #"condRef" #add name
      } 
       
      #NEED TO REMOVE ALREADY EXISTING COMPARISON AND EXCHANGE WITH NEW RESULTS:
      ind = which(matchlist0[,1]==mixsel & matchlist0[,2]==ref)  #get index in matchlist0 (hp stored in this index)
      if(length(ind)==0) { #NB: this variant was not computet before
        ind = which(matchlist1[,1]==mixsel & matchlist1[,2]==ref)  #get index in matchlist0 (hp NOT stored in this index)
        newrow = matchlist1[ind,]
      } else { #Found
        newrow = matchlist0[ind,-6] #don't include model type
      }
      if(condREF || hasCondREF) newrow = c(newrow,"") #add column if condRef considered OR before
      newrow[4] = signif(ret$mleLR,digits=4) #round
      newrow[5] = nC2 #insert number of contributors
      if(condREF) newrow[6] = paste0(condsel,collapse="/") #show conditionalRefs in table
      if(is.null(matchlist2) && condREF) names(newrow)[6] = L$condref #"condRef" #need to add name
      add = FALSE #should be added to list
      insInd= 1 #give index of inserted row in matchlist2  (also used to fast recognize inserted)
      if(!is.null(matchlist2)) {
        ind2 = which(matchlist2[,1]==mixsel & matchlist2[,2]==ref) #find index in matchlist2 (LR placed in this index)
        if(length(ind2)>0) {
         matchlist2[ind2,] = newrow #insert row
         insInd = ind2 #index of inserted row 
        } else {
         add = TRUE
         insInd = nrow(matchlist2)+1 #added as last row
        }
      } else {
         add = TRUE      
      }
      if(add)   matchlist2 = rbind(matchlist2,newrow) #add to table if not already there
      ord = order(as.numeric(matchlist2[,4]),decreasing=TRUE) #get order of table wrt LR vals
      matchlist2 = matchlist2[ord,,drop=FALSE] #sort table wrt LR values
      assign("resCompLR2",matchlist2,envir=nnTK)  #store sorted matchlist2 
    
      #New in v1.2.1: Update resCompLR object 
      matchlist0[ind,4] <- ret$mleLR
      matchlist0[ind,5] <- nC2 
      matchlist0[ind,6] <- "quan" #update model type
      assign("resCompLR",matchlist0,envir=nnTK)  
    
      #STORE HPFIT RESULTS
      hpfitlist = get("storedFitHp",envir=nnTK)  #get already stored objects
      if(is.null(hpfitlist)) hpfitlist = replicate(nrow(matchlist1),list()) #init. list if first time
      hpfitlist[[ind]] <- ret$fithp  #insert on right index
      assign("storedFitHp",hpfitlist,envir=nnTK)  #store object
    
      #UPDATE TABLES:
      sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
      refreshTabLIST2(which(ord==insInd),sort=sortTypes[4]) #update tables with results, with marked on selected one
      createMatchlist(modtype=3) #final step is to update MIXTURES table
      refreshTabLIST(evidsel=mixsel,sort=sortTypes[5]) #update match-tables: mark on considered evidence  
    
      gWidgets2::svalue(nb) <- 4 #go to quan LR result tab when done
    }) #end button
  }
  
  #FUNCTION TO SPECIFY HYPOTHESIS AND PERFORMS A SINGLE DC
  createHypDCWindow = function(evids,refs,nC,lrvals=NULL) { 
    Krange = 1:4 #range of number of contributors
    selind = which(nC==Krange)
    if(length(selind)==0) selind = max(Krange)
    selwin <- gWidgets2::gwindow(paste0( L$deconvolution,"/", L$expPHplot ), visible=TRUE)
    tabSel <- gWidgets2::glayout(container=selwin)
    tabSel[1,1] = gWidgets2::glabel(text=paste0( L$evidences ,colonsymbol),container=tabSel)
    tabSel[1,2] = gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
    tabSel[2,1] = gWidgets2::glabel(text=paste0( L$numcontr ,colonsymbol),container=tabSel)
    tabSel[2,2] = gWidgets2::gcombobox(items=Krange,selected=selind,editable=TRUE,container=tabSel) 
    # gWidgets2::gedit(text=nC,container=tabSel)
    tabSel[3,1] = gWidgets2::glabel(text=paste0( L$condtionON,colonsymbol),container=tabSel)
    tabSel[3,2] = gWidgets2::gcheckboxgroup(items=refs,container=tabSel,checked=TRUE)
    if(!is.null(lrvals)) { #if lr values given
     tabSel[2,3] = gWidgets2::glabel(text= L$logLR ,container=tabSel) #log10 LR
     tabSel[3,3] = gWidgets2::glabel(lrvals,container=tabSel)
    }
    tabSel[4,1] = gWidgets2::gbutton(text= L$calculate ,container=tabSel,handler=function(h,...) {
      nC2 <- as.integer(gWidgets2::svalue(tabSel[2,2])) #get selected references
      refsel <- NULL #conditional references
      if(length(refs)>0) {
       refsel <- gWidgets2::svalue(tabSel[3,2]) #get selected refs
      }
      mixsel <- gWidgets2::svalue(tabSel[1,2]) #get selected samples
      if(length(mixsel)==0) {
        gWidgets2::gmessage( L$msg.selectprofile )
        return()
      }
      ok = setPopFreq()
      if(!ok) return()
      gWidgets2::dispose(selwin) #remove window
      candtab <- doDC(nC=nC2,evids=mixsel,refs=refsel) #perform DC   
      if(nrow(candtab)>0) { #if more than one new 
       dclist <- get("DClist",envir=nnTK) #get stored DC-list
       dclist <- rbind(dclist,candtab) #add candidates
       assign("DClist",dclist,envir=nnTK) #get stored DC-list  
       refreshDCLIST() #refresh DC-list
       gWidgets2::svalue(nb) <- 5 #go to DC-tab
      }
    }) #end button
  } #end function createHypDCWindow
  
  #Functions executed when double clicked on the Mixture list: 
  #1) Show expected peak heights (given as a plot)
  #2) Show Deconvoluted candidates (added to deconvoluted reference list)
  clickmixlist = function(h,...) {
    id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
    tab <- get("allMixList",envir=nnTK)  #load from envir
    if(is.null(tab)) return() #return if no elements
    evid <- tab[id,1] #get evidence
    refs <- unlist(strsplit(tab[id,2],"/"))  #get refs
    #Get corresponding LR values
    LRtab <- get("resCompLR",envir=nnTK)  #load from envir
    LRtab <- LRtab[LRtab[,1]==evid,,drop=FALSE]
    lrvals <- format(as.numeric(LRtab[match(refs,LRtab[,2]),4]),digits=4)
    nC <- tab[id,3] #get number of contributors (should be editable)
    createHypDCWindow(evid,refs,nC,lrvals) #create window for specify model for DC
  }
  
  #Functions executed when double clicked on the DC list: 
  #2) Show Deconvoluted candidates (added to deconvoluted reference list)
  addDCprofile = function(h,...) {
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   suppressWarnings({
     id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
   })
   if(is.na(id)) id <- as.integer(gsub("#","",gWidgets2::svalue(DClistGUI)))
   
   if(length(id)>1) {
   	gWidgets2::gmessage( L$msg.selectoneprofile )
     return()
   }
   answ <- gWidgets2::gconfirm(paste0( L$msg.addDCprofile ,"\n",DClist[id,1]))
   if(answ) { #if extracting DCed candidate
     DCrow = DClist[id,] #obtain DC row to extract
     f_addref(h=list(action=DCrow)) #open edit window of references
  
     #Store extracted DC-result to report object:
     refTab = get("refDataTABLE",envir=nnTK)
     DCtab = rbind(get("DClistReport",envir=nnTK),DCrow) #create DC-tab
     rownames(DCtab)[nrow(DCtab)] = rownames(refTab)[nrow(refTab)] #obtain last extracted ref-name and insert to row
     
     #UPDATING DClistReport object (used for reporting):
     assign("DClistReport",DCtab,envir=nnTK)
     #gWidgets2::svalue(nb) <- 1 #go to data-tab     
   }
  }
  
  #Helpfunction to show LRper marker when clicked
  showLRperMarker = function(h,...) {
    resWOEeval = get("resWOEeval",envir=nnTK) #obtain object
    if(is.null(resWOEeval)) return()
    suppressWarnings({
      ids <- as.integer(gsub("#","",gWidgets2::svalue(h$obj))) #svalue gives name of button if pressed, otherwise its the row in table
    })
    if(is.na(ids)) ids <- as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI)))
    
    s0 = 2 #signif level
    for(id in ids) {
      hypsettxt = paste0(L$lrpermarker,": Hyp. set #",id)
      tab <- cbind(round(resWOEeval[[id]]$mleLRi,s0))  #obtain LR per marker
      colnames(tab) = L$LR
      showGDFtable(hypsettxt, casesolver::addRownameTable(tab,samplename = L$Marker))
    }
  }
  
  #function which takes all matches (with LR>threshold) and create a list to double click on (showing confirming under all conded)
  createMatchlist = function(modtype) { #directly after calculations are done
  #modtype: 0=MAConly(noLR), 1=All Qual LR, 2=All Quan LR, 3=Original Qual LR, but some updated with Quan LR
  threshLR <- get("setupThresh",envir=nnTK)$LRthresh1 #QualLR used by default
  tab <- get("resCompLR1",envir=nnTK) #QualLR used by default
  if(modtype==2)  { 
    threshLR <- get("setupThresh",envir=nnTK)$LRthresh2 #QuanLR used if type 2
    tab <- get("resCompLR2",envir=nnTK)
  }
  if(modtype==0)  { #if only MAC were used
    threshMAC <- get("setupThresh",envir=nnTK)$MACthresh #MAC used if type 0
    tab <- get("resCompMAC",envir=nnTK)$MatchList
    tab <- tab[as.numeric(tab[,3])>=threshMAC,,drop=FALSE] #combinations to consider (all above MAC treshold AND above threshold used in Compare)
  } else { #IF LR was calculated (MAC also stored here)
    score <- as.numeric(tab[,4])
    tab <- tab[score>=log10(threshLR),,drop=FALSE] #combinations to consider (all above LR treshold)
  }
  if(nrow(tab)==0) return()
  if(modtype==3) {
    threshLR2 <- get("setupThresh",envir=nnTK)$LRthresh2 #QualLR used if type 1
    tab2 <- get("resCompLR2",envir=nnTK) #must create a concensus table of Qual/Quan based (using both thresholds)
    if(nrow(tab2)>0) {
     for(rr in 1:nrow(tab2)) { #for each row
       checkind = which( tab[,1]==tab2[rr,1] & tab[,2]==tab2[rr,2]) #find corresponding comparison
       if(length(checkind)==0) {
          if(as.numeric(tab2[rr,4])>=log10(threshLR2)) tab = rbind(tab,tab2[rr,])  #if not found AND it is above threshold, we add it to the list      
       } else {
        	if( as.numeric(tab2[rr,4])<log10(threshLR2)) {
         	 tab = tab[-checkind,,drop=FALSE] #remove from list
       	} else {
           tab[checkind,4:5] = tab2[rr,4:5] #update tab if still keeped
        	}
       }
     } #end for each rr
    } #if any QUAN based LRs
  } #end model type =3
  unEvid <- unique(tab[,1]) #get unique evidence
  unRef <- unique(tab[,2]) #get unique references
  if(length(unEvid)>0) { #if any matches
   resMatches <- numeric() #create match table
   for(evid in unEvid) { #for each evidence we condition on all evidence 
    evidind <- tab[,1]==evid
    refs <- tab[evidind,2] #get references
    if(modtype==0) {
     nC <- sapply(get("mixDataLIST",envir=nnTK)[evid],function(elem) { ceiling(max(sapply(elem,function(x) length(x$adata)))/2) })
    } else {
     nC <- as.numeric(tab[evidind,5][1]) #get number of contributors (equal for all)
    }  
    resMatches <- rbind(resMatches, c(evid,paste0(refs,collapse="/"),nC) )
   }
  } else {
   resMatches =  matrix(nrow=0,ncol=3)
  }
  colnames(resMatches) <- c( L$evidence , L$references , L$numcontr )
  assign("resMatches",resMatches,envir=nnTK) 
  }
  
  #Show matches in a graph:
  #UPDATED IN v1.5.0: Uses plotly to get interactive plot
  f_showMatchNetwork = function(h,...) {
   require(igraph)
  
   createInteractive = FALSE
   if(!is.null(h)) createInteractive = TRUE
   action = "all" #indicate action of plot ("all, onlymix, onlyss)
   if(!is.null(h)) { #IF BUTTON CLICKED: 
     if(h$action=="onlymix") action = "onlymix" #tab = tab[tab[,5]!="1",,drop=FALSE]  # Check if comparing only against Mixture
     if(h$action=="onlyss") action = "onlyss" #tab = tab[tab[,5]=="1",,drop=FALSE]  # Check if comparing only against Mixture
   }
   casesolver::showMatchNetwork(nnTK,action,createInteractive)
   return(TRUE)
  }
  
  ############
  #COMPARISON#
  ####################################
  #STEP 1) CALCULATE MATCHING ALLELES#
  #Function executed when clicking "Comparison"
  getMatchesMAC = function(locs=NULL) { #Compare only alleles (given selection of loci
  DBmix <- get("mixDataTABLE",envir=nnTK) #consider lists
  DBref <- get("refDataTABLE",envir=nnTK) #consider lists
  DBmixmatch <- get("mixDataMATCHSTATUS",envir=nnTK) #consider lists
  
  #Add evidence profiles as "mixture" OR considered as unknown 
  #ONLY EVALUATE PROFILES WHICH HAS NOT MATCHED A REFERENCE (Mixtures and EMpty)
  indUse = DBmixmatch=="mixture" | DBmixmatch==""  #not L$mixture ?
  #if(get("setupAdvanced",envir=nnTK)$compSS=="TRUE")  indUse = indUse | TRUE # grepl("Unknown ",DBmixmatch) #search only those with matchstatus "Unknown"
  DBmix <- DBmix[indUse,,drop=FALSE] #All evidences must be mixtures OR assigned to unknown
  
  if(!is.null(locs)) {
   keep = toupper(colnames(DBmix))%in%toupper(locs)
   DBmix <- DBmix[,keep,drop=FALSE] #use relevant loci
   DBref <- DBref[,keep,drop=FALSE] #use relevant loci
  }
  
  #Perform the Matching allele counting algorithm-> Output is score for all comparisons
  matchMAC <- calcMACcomparison(DBmix=DBmix,DBref=DBref,threshMAC=get("setupThresh",envir=nnTK)$MACthresh) #calculating the score=normalized number of allele match counting for all combinations
  #  colnames(samplename
  
  #Post-process: Remove candidates that are same profile (early stage): BASEDO ON MATCH STATUS???
  if(!is.null(matchMAC) && nrow(matchMAC$MatchList)>0) { #if any comparison results
   DBmixmatchConsider = DBmixmatch[grepl( L$unknown ,DBmixmatch)] #profiles to consider
   DBmixmatchConsiderEvid = names( DBmixmatchConsider)
   rmind = rep(FALSE,nrow(matchMAC$MatchList)) #indices to remove
   for(i in 1:length(DBmixmatchConsider) ) {
    ind = DBmixmatchConsiderEvid[i] == matchMAC$MatchList[,1] &  DBmixmatchConsider[i] == matchMAC$MatchList[,2] 
    rmind[ind] = TRUE #index to remove 
   }
   matchMAC$MatchList = matchMAC$MatchList[!rmind,,drop=FALSE] #remove matches
  } #end If comparison results
  assign("resCompMAC",matchMAC,envir=nnTK)  #store match-matrix in environment 
  } #end MAC comparison
  
  ############################################
  #STEP 2) CALCULATE LR FOR REMAINING MATCHES#
  getMatchesLR = function(type="quan") { #Calculate LR for individuals in MatchList 
  #type={"qual","quan"} 
  Clist <- get("resCompMAC",envir=nnTK)$MatchList  #list to consider for calculating LR
  
  if(type=="quan" && !is.null(get("resCompLR1",envir=nnTK)) ) {
   Clist <- get("resCompLR1",envir=nnTK)  #list to consider for calculating LR (based on qualLR)
   Clist <- Clist[as.numeric(Clist[,4])>=log10(get("setupThresh",envir=nnTK)$LRthresh1),,drop=FALSE] #keep only variants above thrshold AND COLUMNS in MAC
  }
  DBmix <- get("mixDataLIST",envir=nnTK)[unique(Clist[,1])] #get relevant evidence
  DBref <- get("refDataTABLE",envir=nnTK)
  DBref = DBref[rownames(DBref)%in%unique(Clist[,2]),,drop=FALSE] #get only relevant references
  
  mod <- casesolver::getModelSettings(nnTK)  #object for model settings
  
  #matchlist=Clist;popFreq=mod$popFreq;pC=mod$pC; maxC=get("setupAdvanced",envir=nnTK)$maxC1;useMinK1=as.logical(get("setupAdvanced",envir=nnTK)$useMinK1);normalize=mod$normalize;minFreq=mod$minFreq
  suppressWarnings({
   if(type=="qual")  matchLRres <- casesolver::calcQualLRcomparison(DBmix,DBref,matchlist=Clist,popFreq=mod$popFreq,pC=mod$pC, maxC=get("setupAdvanced",envir=nnTK)$maxC1,useMinK1=as.logical(get("setupAdvanced",envir=nnTK)$useMinK1),normalize=mod$normalize,minFreq=mod$minFreq)
   if(type=="quan") {
     
      nContr = NULL #Number of contributors to use (default is the max(nA)/2 rule)
      if(ncol(Clist)==5) nContr = Clist[,5] #number of contributors given in last column if qualLR calculated
  
      #use "Rule of three" for EFM model when applied to SNPs
      if(!is.null(get("setupAdvanced",envir=nnTK)$isSNP) && get("setupAdvanced",envir=nnTK)$isSNP=="TRUE")  nContr = rep("3",nrow(Clist))
  #matchlist=Clist[,1:3,drop=FALSE];popFreq=mod$popFreq;kit=mod$kit;xiBW=mod$xiBW;xiFW=mod$xiFW;pC=mod$pC;lambda=mod$lambda;threshT=mod$threshT;nDone=mod$nDone;maxC=get("setupAdvanced",envir=nnTK)$maxC2;normalize=mod$normalize;minFreq=mod$minFreq
      matchLRres <- casesolver::calcQuanLRcomparison(DBmix,DBref,matchlist=Clist[,1:3,drop=FALSE],popFreq=mod$popFreq,kit=mod$kit,xiBW=mod$xiBW,xiFW=mod$xiFW,pC=mod$pC,lambda=mod$lambda,threshT=mod$threshT,nDone=mod$nDone,maxC=get("setupAdvanced",envir=nnTK)$maxC2,nContr=nContr, normalize=mod$normalize,minFreq=mod$minFreq) 
   }
  })
  matchlist <- matchLRres$MatchList
  matchlist2 = cbind(matchlist,type) #Added in v1.2.1
  assign("resCompLR",matchlist2,envir=nnTK)  #store results from comparison
  
  if(type=="qual") {
   #sort the matchlist with respect to the LRs:  
   LRcol <- which(colnames(matchlist)=="log10LR") #get column where LR is
   LRval <- as.numeric(matchlist[,LRcol])
   ord <- order( LRval,decreasing=TRUE)
   matchlist[,LRcol] <- round(LRval,2) #round to 2 dec
   #matchlist[ord,] #sort list by LR
   assign("resCompLR1", matchlist[ord,,drop=FALSE],envir=nnTK)  #store sorted matchlist 
  }
  if(type=="quan") {
   assign("storedFitHp",matchLRres$storedFitHp,envir=nnTK)  #store model results under Hp
  
   #sort the matchlist with respect to the LRs:  
   LRcol <- which(colnames(matchlist)=="log10LR") #get column where LR is
   LRval <- as.numeric(matchlist[,LRcol])
   ord <- order( LRval,decreasing=TRUE)
   matchlist[,LRcol] <- round(LRval,2) #round to 2 dec
   #matchlist[ord,] #sort list by LR
   assign("resCompLR2", matchlist[ord,,drop=FALSE],envir=nnTK)  #store sorted matchlist 
  }
  } #end LR comparison
  
  
  #Function giving window for editing alleles for new references (returns ref-name)
  f_addref = function(h,...) { 
    refT <- get("refDataTABLE",envir=nnTK)
    refTNames = rownames(refT)
    if(nrow(refT)>0) {
     locs <- colnames(refT) #get loci names from ref-table
    } else {
     locs <- colnames(get("mixDataTABLE",envir=nnTK)) #get loci names from mix-table (in case of no refs)
    }
    refind = 1:nrow(refT) #index of refs to consider (may be empty also)
    if( length(refTNames)>=nLarge) { #the user must give a segment of refs
      ret = gWidgets2::ginput( L$msg.indexinput , text = "1-10", title = "User input", icon = "info")
      ret = unlist(strsplit(ret,","))
      refind = numeric()
      for(rng in ret) {
       tmp = as.integer(unlist(strsplit(rng,"-")))
       if(length(tmp)==1) refind = c(refind , tmp)
       if(length(tmp)==2) refind = c(refind , tmp[1]:tmp[2])
      }
      refind = refind[refind>=1 & refind<=nrow(refT) ]  #require valid range
      refT <- refT[refind,,drop=FALSE]
    }
    def <- "A/B" #default genotype
    defN <- L$name #default name
    refN <- rownames(refT) #refnames
    newTab <- cbind(refN,refT)
    colnames(newTab)[1] <- L$samplename
    
    newline <- c(defN,rep(def,length(locs)))
    if(!is.null(h$action)) { #opened with deconvolution
      tmp <- h$action #get sample name
      newline[1] <- tmp[1] #get sample name (component)
      dat <- tmp[5:length(tmp)] #get suggested deconvolved elem
      newline[-1] <- dat[match(locs,names(dat))] #obtain correct order (notice that only 1st match element is used)
      newline[is.na(newline)] <- ""
    }
    newTab <- rbind(newTab,newline) #add empty/new row
    
    #CREATE GUI WINDOW (OTHER PROCESS WILL STOP)
    flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
    setwin <- gWidgets2::gwindow(paste0( L$data.editrefs ),visible=FALSE) 
    gWidgets2::addHandlerUnrealize( setwin, handler = function(h,...) { #
      tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
      return(NULL) #return selected profiles when window is exit 
    }  ) #call quit function
    
    tabval = gWidgets2::ggroup(spacing=0,container=(setwin),horizontal=FALSE)  
    guitab <- gWidgets2::gdf(items=newTab,container = tabval) 
    gWidgets2::add(tabval,guitab,expand=TRUE,fill=TRUE) 
    
    gWidgets2::gbutton(text= L$save ,spacing=0,container=tabval,handler = function(h, ...) { 
      if(nrow(refT)==0) {
        newref <- t(as.character(unlist(guitab[])))
      } else {
        newref <- as.matrix(guitab[])
      }
      nL <- nrow(newref) #rowindex of new ref
      addRef <- newref[nL,] #reference to add
      newref = newref[-nL,,drop=FALSE] #keep only tab without added ref
      delindsOld <- newref[,1]=="" #index of all refs to delete (When names are set to "")
    
      if(any(delindsOld)) { #if deleting previous stored refs
       answ <- gWidgets2::gconfirm(paste0( L$msg.suredelete ,":\n",paste0(refN[delindsOld],collapse="/"),"?"))
       if(answ) { #if agree then refs are deleted
        refT <- refT[!delindsOld,,drop=FALSE] #update ref-table
        newref = newref[!delindsOld,,drop=FALSE] #BUG fixed in v1.4.0
        refN = refN[!delindsOld]
       }
      } 
    
      #CHECK FOR CHANGES (not added ref)
      if(nrow(refT)>0) { #must have at least one ref
       isChanged = which(refT != newref[,-1],arr.ind=TRUE)
       if(length(isChanged)>0) {
        text = paste0(refN[isChanged[,1]]," (",colnames(refT)[isChanged[,2]],")")
        answ <- gWidgets2::gconfirm(paste0( L$msg.sureapplychanges ,":\n",paste0(text,collapse="/"),"?"))
        if(answ) { #apply changes:
         refT <- newref[,-1,drop=FALSE]  #update ref-table
         rownames(refT) = refN  #add rownames
        } else {
         newref[,-1] = refT #replace with prev
        }
       }
      }
      sn = addRef[1]
      Anew2 = NULL
      if(!sn%in%c("",defN)) { #only include ref if a new sample name were given
       Anew <- addRef[-1] #new alleles
       Anew[Anew==def] <- "" #Defaults are set as empty
       if(sn%in%refTNames) {
        gWidgets2::gmessage( L$msg.refnametaken )
       } else {
         if( any(Anew!="")) Anew2 = Anew
       }
      } #end if insert ref
    
      #In case the user gave a segment of refs
      if( length(refTNames)>=nLarge) {
         refT2 = get("refDataTABLE",envir=nnTK) #get full table 
         if(length(isChanged)>0) { #if anything was changed
          refT2[refind[!delindsOld],] = refT
         }
         if(any(delindsOld)) { #if any deleted
           refT2 = refT2[-refind[delindsOld],,drop=FALSE] #indices to be removed
         } 
         refT = refT2 #update
         rm(refT2);gc()
      } #end selected segment
    
      #Add to table if something to add:
      if(!is.null(Anew2)) {
       refT <- rbind(refT,Anew2) #add to existing table
       rownames(refT)[nrow(refT)] <- sn
       
       
      }    
      assign("refDataTABLE",refT,envir=nnTK) #store table 
      gWidgets2::dispose(setwin) #close window
      tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
      
      sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
      updateTables(type="ref",sort=sortTypes[2]) #updates tables again (with same sorting)
    }) #end button
    gWidgets2::visible(setwin) <- TRUE
    tcltk::tkwait.variable(flag) #important to not quit window before broken
  }  #end add ref
  
  #FUNCTION TO STORE TABLE:
  f_exporttable = function(h,...) { 
    tab <- NULL
    if(h$action=="mac")    tab <- casesolver::addRownameTable(get("resCompMAC",envir=nnTK)$MatchMatrix,type=1,L$samplename)
    if(h$action=="qual")   tab <- get("resCompLR1",envir=nnTK)
    if(h$action=="quan")   tab <- get("resCompLR2",envir=nnTK)
    if(h$action=="final")  tab <- get("allMixList",envir=nnTK)
    if(h$action=="dc")     tab <- get("DClist",envir=nnTK)
    if(h$action=="woe")    tab <-  get("resWOEeval",envir=nnTK)$resTable
    saveTable(tab,"csv")
  }
  
  #Function to deconvolve for all mixtures (conditions on matching references)
  f_doDCall = function(h,...) {
    mixList <- get("allMixList",envir=nnTK) #get match list (from function)
    if(is.null(mixList)) return() #return if no list found
    nS = nrow(mixList) #get number of mixtures 
    if(nS==0) return()
    candtabs = numeric() 
    for(ss in 1:nS) {
      evid = mixList[ss,1]
      refs = unlist(strsplit(mixList[ss,2],"/"))
      if(length(refs)==0) refs = NULL 
      nC = as.integer(mixList[ss,3])
    
      #UPDATED v1.4.1 to fix possible crash.
      if(nC<=length(refs)) {
       print(paste0("No unknowns to deconvolve for sample ",evid)) 
      } else {
       print(paste0("Deconvolving sample #",ss,": ",evid)) 
       candtab <- doDC(nC,evid,refs,showPlot=FALSE,useplotly=FALSE)
       candtabs = rbind(candtabs,candtab) #add candidate table
       print("")
       print(paste0(round(ss/nS* 100), "% DC calculation complete...")) 
      } 
    }
    if(length(candtabs)>0) { #if more than one new 
     dclist <- get("DClist",envir=nnTK) #get stored DC-list
     dclist <- rbind(dclist,candtabs) #add candidates
     assign("DClist",dclist,envir=nnTK) #get stored DC-list  
     refreshDCLIST() #refresh DC-list
     gWidgets2::svalue(nb) <- 6 #go to DC-tab after calculation
    }
  } #end f_doDCall
  
  #Function to do WOE as a final step from matchlist
  f_doWOEcall = function(h,...) {
   
   #CHECK AND SET FREQUENCY FILE BEFORE
   ok = setPopFreq(giveMessage=TRUE) #import population frequency from last selected file
   if(!ok) return() #retun from function if not frequencies are set.
  
   #perform WOE calculation (returns from function when closed/finished)
   casesolver::calcWOEhyps(nnTK) 
   #Update with WoE table when done evaluated
   #length( get("resWOEeval",envir=nnTK) )
   resList = get("resWOEeval",envir=nnTK) #obtain results
   #object.size(resWOEeval)/1e6 #size of object
   
   if(length(resList)==0) return() #return if no elements
   extractrow = function(x) {
     s0 = 2
     evidtxt = paste0(x$evid,collapse="/")
     condtxt = paste0(x$cond ,collapse="/")
     validtxt = paste0(x$nFailedHp,"/",x$nFailedHd)
     c(evidtxt,x$poi,condtxt,x$NOC,round(x$mleLR,s0),round(x$bayesLR,s0),round(x$consLR,s0),round(x$MxPOI,s0),validtxt)
   }
   resTable = t(sapply(resList,  extractrow))
   colnames(resTable) = c(L$evidences, L$POI, L$Conditionals, L$NOC, L$LRmle, L$LRbayes, L$LRcons, L$MxPOI, L$nFailed)
   resList$resTable =resTable #insert table
   
   print("Done with Weight-of-evidence calculations!")
   assign("resWOEeval",resList,envir=nnTK) #obtain results
   refreshWOELIST()
   gWidgets2::svalue(nb) <- 7  #change panel
   setFocus() #refocus
  }
  
  
  ##################################################################################################
  ########### Program starts #######################################################################
  ##################################################################################################
  
  
  #SELECT LANGUAGE (before anything else sat)
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  
  #software name:
  softname <- paste("CaseSolver",L$version,version)
  
  #DEFINE REPORT ITEMS BASED ON LANGUAGE:
  reportitems = rep(NA,nReportItems)
  reportitems[1] = L$header #"Header"
  reportitems[2] = L$references  #"References"
  reportitems[3] = L$extractedprofiles #Extracted profiles"
  reportitems[4] = paste0( L$singlesources ," (", L$alleles,")")  #"Single sources (alleles)"
  reportitems[5] = paste0( L$mixtures ," (", L$alleles,")") #"Mixtures (alleles)"
  reportitems[6] = paste0( L$consensusprofiles ," (", L$alleles,")") #"Consensus (alleles)"
  reportitems[7] = paste( L$singlesources , L$withPH ) #"Single sources w/PH"
  reportitems[8] = paste( L$mixtures , L$withPH ) #"Mixtures w/PH"
  reportitems[9] = L$metadata  # "Metadata"
  
  #Comparisons:
  reportitems[10] = L$matchmatrix #"Match matrix"
  reportitems[11] = paste0( L$matchlist ," (", L$qualLR,")")  #"MatchList (Qual)"
  reportitems[12] = paste0( L$matchlist ," (", L$quanLR,")")  #"MatchList (Quan)"
  reportitems[13] = paste( L$finalmatchlist )  #"FinalMatchList"
  reportitems[14] = paste( L$matchnetwork )  #""MatchNetwork"
  reportitems[15] = paste0( L$RMNE ," (", L$evidence,")" ) #"RMNE for evidence" 
  reportitems[16] = paste0( L$RMP ," (", L$reference,")" ) #"RMP for references"
  reportitems[17] = paste( L$concordantevidences ) #"Concordant Evidences"
  reportitems[18] = paste0( L$IBS ," (", L$reference,")" ) #"IBS for references"
  
  #Quantifications:
  reportitems[19] = L$deconvoluted #Deconvoluted
  reportitems[20] = L$woeResult #WoE table
  reportitems[21] = L$woeStatements #Statements
  reportitems[22] = L$woeParams #Estimated model parameters (listed per hypothesis)
  reportitems[23] = L$lrpermarker #LR per marker (woe)
  
  #attachments:
  reportitems[24] = L$settings #"Settings"
  reportitems[25] = paste0(  L$EPG, " (", L$singlesources ,")") #"EPG figures for single sources"
  reportitems[26] = paste0(  L$EPG, " (", L$mixtures ,")") # "EPG figures for mixtures"
  #reportitems[25] = L$expPHplot #"Expected PH plots"
  if(length(reportitems)!=nReportItems) stop("Number of report items were wrong.")
  
  #storing report item names to be used in report:
  optReport <- get("setupReportLay",envir=nnTK)  #get report options
  optReport$reportitems = reportitems #store name of reportnames
  assign("setupReportLay",optReport,envir=nnTK) #store back to environment
  
  #Get ON-OFF choices (used in f_createreport and f_modelsel)
  radiotxt = c( L$on , L$off )
  
  ###############
  #start layout:#
  ###############
  mblst = list() 
  mblst[[ L$file ]] = list(  
   gWidgets2::gaction( paste( L$set , L$workdir ) ,handler=f_setwd),
   gWidgets2::gaction( paste( L$open , L$proj ) ,handler=f_openproj),
   gWidgets2::gaction( paste( L$save , L$proj ) ,handler=f_saveproj),
   gWidgets2::gaction( L$quit ,handler=f_quitproj,icon="close")
  )
  mblst[[ L$setup ]] =list( 
   gWidgets2::gaction( paste( L$set , L$threshsettings ) ,handler=f_threshsel),
   gWidgets2::gaction( paste( L$set , L$modelsettings )  ,handler=f_modelsel),
   gWidgets2::gaction( paste( L$select ,L$kit) ,handler=f_kitsel),
   gWidgets2::gaction( paste( L$select ,L$popfreq) ,handler=f_popsel),
   gWidgets2::gaction( paste( L$select , L$importfun) ,handler=f_importsel),
   gWidgets2::gaction( paste( L$select , L$pathcasefolders ) ,handler=f_casedirsel)
  )
   mblst[[ L$report ]] =list( 
   gWidgets2::gaction( paste( L$set , L$reportlayout ) ,handler=f_reportlay),
   gWidgets2::gaction( paste( L$set , L$report, L$options) ,handler=f_reportopt),
   gWidgets2::gaction( paste( L$set , L$report, L$export, "Types" ) ,handler=f_reportexptyp),
   gWidgets2::gaction( paste( L$set , L$report, L$export, L$options) ,handler=f_reportexpopt),
   gWidgets2::gaction( paste( L$set , L$Markernames) ,handler=f_reportlocnames)
  )
  mblst[[ L$advanced ]] = list( 
   gWidgets2::gaction( paste( L$advanced , L$options ) ,handler= f_advancedoptions),
   gWidgets2::gaction( paste( L$mcmc , L$options ) ,handler= f_mcmcoptions),
   gWidgets2::gaction( L$randomibs ,handler=f_calcIBSdist), #calculates random IBS
   gWidgets2::gaction( paste( L$select , L$language ) ,handler=f_selLanguage) #select language
  )
  
  
  #change working directory to the one stored in nnTK-environment
  wd=get("workdir",envir=nnTK) #assign working directory to nnTK-environment
  if(!is.null(wd)) {
  tryCatch( { setwd(wd) }, error = function(e) print("Warning: Workdirectory could not be changed!") )
  }
  
  #############
  #Main window#
  #############
  mainwin <- gWidgets2::gwindow(softname, visible=FALSE)#, width=mwW,height=mwH)
  gWidgets2::addHandlerUnrealize( mainwin, handler = function(h,...) { return( !gWidgets2::gconfirm( L$msg.quit ) )}  ) #call quit function
  
  gWidgets2::gmenu(mblst,container=mainwin)
  nb = gWidgets2::gnotebook(container=mainwin)
  tabimport = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= L$data ) #tab2: (imports all files)
  tabmatchmatrix = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= paste( L$matchmatrix ) ) #comparison results
  tabmatchlist1 = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= paste0( L$matchlist," (", L$qualLR,")" )) #results from Qualitative model (LRmix)
  tabmatchlist2 = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= paste0( L$matchlist," (", L$quanLR,")" )) #results from Quantitative model (EFM)
  tabmixtures = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= L$matches ) #comparison results (collapsed)
  tabdeconv = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= L$deconvoluted ) #Deconvolution results
  tabWOE = gWidgets2::ggroup(horizontal=FALSE,spacing=spc/2,container=nb,label= L$woe ) #WoE results\
  
  ####################################################
  ###############Tab 1: Import Data:##################
  ####################################################
  
  #TAB layout
  layhor0 = as.logical(get("setupView",envir=nnTK)$importHorizontal) #get configured layout of tables
  tabimportA = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe( paste( L$select , L$caseid ) ,container=tabimport)) #kit and population selecter
  tabimportC = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe( L$functionalities ,container=tabimport)) #Tasks button
  tabimportB = gWidgets2::gpanedgroup(horizontal=layhor0,container=gWidgets2::gframe( paste0( L$data.evidref ),container=tabimport,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
  
  txtEmptycasedir = paste( L$msg.emptycasedir, L$setup, ">", L$set , L$pathcasefolders ) #create a error message for not finding case folders
  #Choose box and import button
  casedir <-  get("setupCase",envir=nnTK)$casepath
  if(is.na(casedir) || casedir=="") gWidgets2::gmessage( txtEmptycasedir )
  casefolds <- list.dirs(casedir,recursive=FALSE, full.names = TRUE)  #keep full path
  casefolds2 <- list.dirs(casedir,recursive=FALSE, full.names = FALSE)  #extract only folder CHanged in v1.7
  
  caseid <- get("caseID",envir=nnTK) 
  if(is.null(caseid)) {
    caseid <- 0 #set a default index
  } else {
    caseid <- which(casefolds2==caseid) #must be index
  }
  tabimportA[1,1] = gWidgets2::gbutton(text= L$changeview ,container=tabimportA,handler=function(h,...) { #CHANGE VIEW
    opt = get("setupView",envir=nnTK) #get options
    if(layhor0) {
      opt$importHorizontal = "FALSE"
    } else {
      opt$importHorizontal = "TRUE"
    }
    assign("setupView",opt,envir=nnTK)  #store to environment
    setupWrite(unlist(opt),file=setupFileView)    #save to file in installation folder
    gWidgets2::dispose(mainwin) #shut down window
    #gc() #empty garbage memory
    gui(envir=nnTK) #restart GUI session with project environment
  })
  
  #BUTTON FOR IMPORTING DATA (calls f_importData)
  if( length(casefolds2)>0 && length(caseid)>0  ) {#v1.4.1: ADDED if to make possible to exchange proj files
    tabimportA[1,2] <- gWidgets2::gcombobox(items=casefolds2, selected = caseid  , editable = TRUE, container = tabimportA) 
    tabimportA[1,3] <- gWidgets2::gbutton(text= L$import ,container=tabimportA,handler=f_importData) #IMPORT FUNCTION!!!
    gWidgets2::tooltip(tabimportA[1,3]) <- L$tip.data.import #add tooltip only if button is visible
  }
  
  mixTabGUI <- gWidgets2::gtable(items="",multiple = TRUE,container = tabimportB, handler=clicktable,action="mix")
  refTabGUI <- gWidgets2::gtable(items="",multiple = TRUE,container = tabimportB, handler=clicktable,action="ref")
  #gWidgets2::add(tabimportB,mixTabGUI,expand=TRUE,fill=TRUE)
  #gWidgets2::add(tabimportB,refTabGUI,expand=TRUE,fill=TRUE)
  gWidgets2::svalue(tabimportB) = 0.5 #set as 50-50 by default
  
  #TABLE SORTING
  sortstring1 <- paste( L$data.sortevid, c("#", L$name , L$matchstatus )) #EVID
  tabimportA[1,4] = gWidgets2::gcombobox(items=sortstring1,selected=1,container=tabimportA,horizontal = TRUE, handler=
  function(h,...) { 
   sortval = which(sortstring1== gWidgets2::svalue(tabimportA[1,4])) #sort value
   resave_Sorting(1,sortval) #resaving to setupSorting object
   updateTables(sort=sortval,type=h$action) 
  },action="mix")
  sortstring2 <-  paste( L$data.sortref, c("#", L$name )) #REFERENCE
  tabimportA[1,5] = gWidgets2::gcombobox(items=sortstring2,selected=1,container=tabimportA,horizontal = TRUE, handler=
  function(h,...) { 
   sortval = which(sortstring2== gWidgets2::svalue(tabimportA[1,5]))
   resave_Sorting(2,sortval) #resaving to setupSorting object
   updateTables(sort=sortval,type=h$action) 
  },action="ref")
  
  #UPDATE IN v1.5: Possible to add references afterwards
  tabimportA[1,6] <- gWidgets2::gbutton(text= L$importref ,container=tabimportA,handler = function(h, ...) { #IMPORT REFS
    ff <- mygfile( paste( L$select , L$file) ,type="open")
    if(length(ff)==0) return()
    tab = euroformix::tableReader(ff) #read table and insert to GUI format: [SampleName,Marker,Allele 1, Allele 2]
    sn = unique(tab[,1]) 
    ln = toupper(unique(tab[,2]))
    if(length(sn)==0) return()
  
    refT <- get("refDataTABLE",envir=nnTK) #store table 
    cn = toupper(colnames(refT))
    if(length(cn)==0) return()
    newsn = c(rownames(refT),sn)
  
    for(ss in sn) { #for each sample
      newT = rep("",length(cn)) 
      for(loc in ln) { #for each loci
        indins = which( cn==loc )
        av = unlist(tab[tab[,1]==ss & tab[,2]==loc,3:4])
      	 if(length(av)==1) av = rep(av,2)
      	 if(length(av)==2) {
         newT[indins] <- paste0(av,collapse="/")  #factor(sn, levels= c(sn,levels(guitab[nR,1]) ))#insert name
      }
     }
     refT = rbind(refT,newT) #add profiles  
    }
    rownames(refT) = newsn 
  
    assign("refDataTABLE",refT,envir=nnTK) #store table 
    
    sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
    updateTables(type="ref",sort=sortTypes[2]) #updates ref-table again 
  })
  
  tabimportC[1,1] = gWidgets2::gbutton(text= L$compare ,container=tabimportC,handler= #COMPARE FUNCTION
  function(h,...) {
  
  refDataTABLE =  get("refDataTABLE",envir=nnTK)
  if(is.null(refDataTABLE) || nrow(refDataTABLE)==0) return() #return if no refs
  
  #Reset all earlier comparison-results (since method may have changed):
  print("Resetting all previous comparison-results...")
  assign("resCompMAC",NULL,envir=nnTK);assign("resCompLR1",NULL,envir=nnTK);assign("resCompLR2",NULL,envir=nnTK)  
  assign("resCompLR",NULL,envir=nnTK);assign("resMatches",NULL,envir=nnTK);assign("allMixList",NULL,envir=nnTK) 
  suppressWarnings({
   matchL1GUI[] = NULL #Set to zero if anything
   matchL2GUI[] = NULL #Set to zero if anything
  })
  
  #CHECK AND SET FREQUENCY FILE
  ok = setPopFreq(giveMessage=TRUE) #import population frequency from last selected file
  if(ok) {
   locsUse = names(get("popFreq",envir=nnTK)) #use markers in freq table if found
  } else {
   locsUse = colnames(get("mixDataTABLE",envir=nnTK)) #otherwise use loci from evid data
  }
  print(paste0(length(locsUse)," loci used in COMPARISON:"))
  print(paste0(locsUse,collapse="/")) #Print to console which loci are used
  
  #Step 1: Calculate MAC
  res <- getMatchesMAC(locs=locsUse)  #get name of loci to consider
  if(is.null(res)) return(); #return if nothing to compare.
  refreshTabMATRIX() #update tables
  gWidgets2::svalue(nb) <- 2 #go to comparison tab
  
  #CALCULATING LR BASED SCORES:
  if(nrow(res$MatchList)==0) return() #Return if no candidates to calculate
  modtype <- get("setupModel",envir=nnTK)$modeltype #model type selected {1="qual",2="quan",3="both"} #Otherwise only MAC based
  if(!ok) modtype = 0 # return() #show MAC results if LR can't be calculated
  
  sortTypes = get("setupSorting",envir=nnTK) #Obtain sort types
  #Step 2 (optional): Calculate qual based LR
  if( modtype%in%c(1,3) ) {
   getMatchesLR(type="qual") #LRmix
   refreshTabLIST1(sort=sortTypes[3]) #update table with results
   gWidgets2::svalue(nb) <- 3 #go to qual LR result tab when done
  }
  
  #Step 3 (optional): Calculate quan based LR
  if( modtype%in%c(2,3) ) {
   if( get("setupModel",envir=nnTK)$degrad==1 && !casesolver::canPrintEPG(nnTK) ) { #If degradation model chosen but kit not selected
    gWidgets2::gmessage( paste( L$msg.kitspecify, L$settings ,">", L$select ,L$kit) )
    return()
   } 
   getMatchesLR(type="quan") #EFM based
   refreshTabLIST2(sort=sortTypes[4]) #update tables with results
   gWidgets2::svalue(nb) <- 4 #go to quan LR result tab when done
  }
  
  #Step 3: Create matchlist (Final)
  createMatchlist(modtype=modtype) 
  refreshTabLIST(sort=sortTypes[5]) #update tables  
  gWidgets2::svalue(nb) <- 5 #go to overview when done
  
  })
  tabimportC[1,2] = gWidgets2::gbutton(text= paste( L$create , L$report) ,container=tabimportC,handler=function(h,...) {
   casesolver::createReport(nnTK) #creating report (separate R-script)
  })
  
  tabimportC[1,3] = gWidgets2::gbutton(text= L$data.select ,container=tabimportC,handler=f_markprofs)
  tabimportC[1,4] = gWidgets2::gbutton(text= paste( L$calculate , L$IBS ) ,container=tabimportC,handler=f_calcIBS)
  tabimportC[1,5] = gWidgets2::gbutton(text= paste( L$calculate , L$RMP ),container=tabimportC,handler=f_calcRMP)
  tabimportC[1,6] = gWidgets2::gbutton(text= L$data.concordance ,container=tabimportC,handler=f_calcEvidConc)
  tabimportC[1,7] = gWidgets2::gbutton(text= L$data.editrefs ,container=tabimportC,handler=f_addref)
  tabimportC[1,8] = gWidgets2::gbutton(text= L$restart ,container=tabimportC,handler=
  function(h,...) {
   gWidgets2::dispose(mainwin) #shut down window
   gui() #start an empty session (recognized)
  })
  
  
  #INSERT DATA (TABLE-FORMAT) TO GUI: NOTICE the clicktable handler 
  updateTables <- function(sort=1,type="both") { #function to call to update tables (possibly changed order)
   gWidgets2::visible(mainwin) = FALSE
   
   mixTab <- refTab <- "" #empty tables by default
   if(type%in%c("both","mix")) {    #Add mix-table
    mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #assigned in nnTK-environment
    if( !is.null(mixDataTABLE) && nrow(mixDataTABLE)>0 ) { #make sure that there are data in table
     mixDataMATCHSTATUS <- get("mixDataMATCHSTATUS",envir=nnTK) #assign to nnTK-environment 
  
     #SORT TABLE:     
     ord = casesolver::orderTableSort(rownames(mixDataTABLE),mixDataMATCHSTATUS,sort)
     
     newtab <- cbind(mixDataMATCHSTATUS,mixDataTABLE)
     mixTab <- casesolver::addRownameTable(newtab,type=2,L$samplename)
     colnames(mixTab)[1:3] <- c(" ", L$samplename,L$matchstatus) #insert column name for table
     mixTab[mixTab[,3]=="mixture",3] = L$mixture #insert name: BEWARE THAT SAMPLENAMES SHOULD NOT CONTAIN "mixture"
  
     mixTabGUI[] <- mixTab[ord,,drop=FALSE] #update order in GUI table
     if(!layhor0 && !is.null(ncol(mixTab)) && nrow(mixTab)<nLarge ) { #if vertical layout and less than nlarge rows
      colw1 = c(30,150,150)
      colL = 100 #column width for each locus
      gWidgets2::size(mixTabGUI) <- list(column.widths=c(colw1,rep(colL,ncol(mixTab)-length(colw1))))
     }
    }
   } 
   if(type%in%c("both","ref")) {    #Add ref-table
    refDataTABLE <- get("refDataTABLE",envir=nnTK) #assign to nnTK-environment
    if(!is.null(refDataTABLE) && nrow(refDataTABLE)>0 ) { #make sure that there are data in table
      
     #SORT TABLE:     
     ord = casesolver::orderTableSort(rownames(refDataTABLE),sort=sort)
      
     refTab <- casesolver::addRownameTable(refDataTABLE,type=2,L$samplename)
     refTabGUI[] <- refTab[ord,,drop=FALSE]
  
     if(!layhor0 && !is.null(ncol(refTab)) && nrow(refTab)<nLarge ) { #if vertical layout and less than nlarge rows
      colw1 = c(30,150,150)
      colw2 = c(colw1[1],sum(colw1[-1]))
      colL = 100 #column widt for each locus
      gWidgets2::size(refTabGUI) <- list(column.widths=c(colw2,rep(colL,ncol(refTab)-length(colw2))))
     }
    } 
   }
   gWidgets2::visible(mainwin) = TRUE
   setFocus()
  } #end update Table
  updateTables() #update when program starts (use default sorting)
  
  #Add tooltips:
  gWidgets2::tooltip(tabimportA[1,1]) <- L$tip.data.changeview
  gWidgets2::tooltip(tabimportA[1,6]) <- L$tip.data.importref
  gWidgets2::tooltip(tabimportC[1,1]) <- L$tip.data.compare
  gWidgets2::tooltip(tabimportC[1,2]) <- L$tip.data.createreport
  gWidgets2::tooltip(tabimportC[1,3]) <- L$tip.data.selectprofiles
  gWidgets2::tooltip(tabimportC[1,4]) <- L$tip.data.ibs
  gWidgets2::tooltip(tabimportC[1,5]) <- L$tip.data.rmp
  gWidgets2::tooltip(tabimportC[1,6]) <- L$tip.data.conc
  gWidgets2::tooltip(tabimportC[1,7]) <- L$tip.data.editref
  gWidgets2::tooltip(tabimportC[1,8]) <- L$tip.data.restart
  
  ####################################################################################################################
  #######################################Tab 2: Match matrix: #############################################################
  #####################################################################################################################
  
  f_rotateMatchMatrix = function(h,...) {
    refreshTabMATRIX(rotate=TRUE) #update with rotated table  
  }
  
  f_truncatevals = function(h,...) { #removes values below threshold table
    val = get("setupThresh",envir=nnTK)$MACthresh #get options
    matchmat <- get("resCompMAC",envir=nnTK)$MatchMatrix
    if(is.null(matchmat)) return()
    if(gWidgets2::svalue(gridTab2[1,6])== L$truncate )  {
     matchmat[matchmat<val] = ""   #remove values below threshold
    }
    
    gWidgets2::visible(mainwin) = FALSE
    matchMATGUI[] <-  casesolver::addRownameTable(matchmat,type=2,L$samplename)
    if(!is.null(ncol(matchMATGUI[])) && all(dim(matchmat)<nLarge) ) gWidgets2::size(matchMATGUI) <- list(column.widths=c(30,rep(100,ncol(matchMATGUI[])-1)))
    gWidgets2::visible(mainwin) = TRUE
    setFocus()
  }
  
  gridTab2 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabmatchmatrix)) #kit and population selecter
  gridTab2[1,1] <- gWidgets2::gbutton(text= L$export, container=gridTab2,handler=f_exporttable,action="mac")  #Button to create matchcloud
  gWidgets2::tooltip(gridTab2[1,1]) <- L$tip.export 
  
  gridTab2[1,2] <- gWidgets2::gbutton(text= L$rotatematrix, container=gridTab2,handler=f_rotateMatchMatrix)  #Button to rotate table
  gridTab2[1,3] <- gWidgets2::gbutton(text=paste( L$sort.by , L$sort.column),container=gridTab2,handler= function(h,...) { refreshTabMATRIX(sort = 2)}) #Sort by colnames   
  gridTab2[1,4] <- gWidgets2::gbutton(text=paste( L$sort.by , L$sort.row),container=gridTab2,handler= function(h,...) { refreshTabMATRIX(sort = 3)}) #Sort by rownames   
  gridTab2[1,5] <- gWidgets2::gbutton(text=paste( L$sort.by , L$sort.matchval),container=gridTab2,handler= function(h,...) { refreshTabMATRIX(sort = 4)}) #Sort by rownames   
  gridTab2[1,6] <- gWidgets2::gradio(items= c( L$donttruncate , L$truncate ) ,selected=1,container=gridTab2,handler=f_truncatevals)  #Button to create matchcloud
  gWidgets2::tooltip(gridTab2[1,6]) <- L$tip.matchmatrix.truncate
  tabCompMatrix <- gWidgets2::ggroup(container=tabmatchmatrix,expand=TRUE,fill=TRUE)
  matchMATGUI <- gWidgets2::gtable(items="",container=tabCompMatrix) #add to frame
  gWidgets2::add(tabCompMatrix,matchMATGUI,expand=TRUE,fill=TRUE)
  
  #POSSIBLE TO CHANGE ROTATE-LAYOUT (ORDERING OF EVID/REFS) TO SHOW IN REPORT
  refreshTabMATRIX = function(rotate=FALSE,sort=1) { 
    resobj = get("resCompMAC",envir=nnTK)
    if(is.null(resobj$MatchMatrix)) return() #return if no match matrix found
    
    if(sort==2) { #sort by column names
      resobj$MatchMatrix = resobj$MatchMatrix[,order( colnames(resobj$MatchMatrix) ,decreasing=FALSE),drop=FALSE]
    } else if(sort==3) { #sort by row names
      resobj$MatchMatrix = resobj$MatchMatrix[ order( rownames(resobj$MatchMatrix) ,decreasing=FALSE),,drop=FALSE]
    }else if(sort==4) { #sort by match degree (both rows/columns)
      colOrder = order( colSums(resobj$MatchMatrix) ,decreasing=TRUE ) #order based on total MAC sum
      rowOrder = order( rowSums(resobj$MatchMatrix) ,decreasing=TRUE ) #order based on total MAC sum
      resobj$MatchMatrix = resobj$MatchMatrix[ rowOrder,colOrder, drop=FALSE]
    }
    if(rotate) resobj$MatchMatrix = t(resobj$MatchMatrix) #ROTATE MATRIX
    assign("resCompMAC",resobj,envir=nnTK) #store manipulated object (keep same form of table in report)
    
    gWidgets2::visible(mainwin) = FALSE
    matchMATGUI[] <-  casesolver::addRownameTable(resobj$MatchMatrix,type=2,L$samplename) #rotate view
    if(!is.null(ncol(matchMATGUI[])) && all(dim(resobj$MatchMatrix)<nLarge)) gWidgets2::size(matchMATGUI) <- list(column.widths=c(30,rep(100,ncol(matchMATGUI[])-1)))
    gWidgets2::visible(mainwin) = TRUE
    
    #truncate if option is selected
    if(gWidgets2::svalue(gridTab2[1,6])==L$truncate) f_truncatevals(NULL) 
    
    setFocus()
  }
  refreshTabMATRIX() #use default sort
  
  ####################################################################################################################
  #######################################Tab 3: Match List (Qual): ##################################################
  ####################################################################################################################
  sortMathListTablesTxt = paste( L$sort.by , c( L$sort.lr , L$sort.evid , L$sort.ref ))
  
  gridTab3 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabmatchlist1)) #kit and population selecter
  gridTab3[1,1] <- gWidgets2::gbutton(text= L$export ,container=gridTab3,handler=f_exporttable,action="qual")  
  gWidgets2::tooltip(gridTab3[1,1]) <- L$tip.export 
  
  gridTab3[1,2] <- gWidgets2::gbutton(text= paste( L$calculate , L$all, L$quanLR) ,container=gridTab3,handler=f_calcQuanLRall)  
  gWidgets2::tooltip(gridTab3[1,2]) <- L$tip.matchlistQual.calc
  
  gridTab3[1,3] <- gWidgets2::gcombobox(items=sortMathListTablesTxt,selected=1,container=gridTab3,handler=
   function(h,...) {
     sortval = which(gWidgets2::svalue(gridTab3[1,3])==sortMathListTablesTxt)
     resave_Sorting(3,sortval) #resaving to setupSorting object
     refreshTabLIST1(sort=sortval) #change sorted order
   }
  )  
  tabCompLIST1 = gWidgets2::ggroup(container=tabmatchlist1,expand=TRUE,fill=TRUE)
  
  matchL1GUI <- gWidgets2::gtable(items="",container=tabCompLIST1,handler=clickmatchlistQUAL)
  gWidgets2::add(tabCompLIST1,matchL1GUI,expand=TRUE,fill=TRUE)#add to frame
  
  refreshTabLIST1 = function(sort=1) { 
    matchlist <- get("resCompLR1",envir=nnTK)
    if(is.null(matchlist)) return() #return if no list found
    if(nrow(matchlist)==0) return() #return if no candidate found
  
    resTab = casesolver::addRownameTable(matchlist,type=3,L$samplename)
    ord = casesolver::orderTableSort(resTab[,2],resTab[,3],sort) #obtain order of sorting
    
    matchL1GUI[] <-  resTab[ord,,drop=FALSE]
    if( nrow(resTab)<nLarge ) gWidgets2::size(matchL1GUI) <- list(column.widths=c(30,300,300,70,70,100))
  }
  refreshTabLIST1() #use default sort
  
  ####################################################################################################################
  #######################################Tab 4: Match List (Quan): ##################################################
  ####################################################################################################################
  
  gridTab4 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabmatchlist2)) #kit and population selecter
  gridTab4[1,1] <- gWidgets2::gbutton(text= L$export ,container=gridTab4,handler=f_exporttable,action="quan")#
  gWidgets2::tooltip(gridTab4[1,1]) <- L$tip.export 
  
  gridTab4[1,3] <- gWidgets2::gcombobox(items=sortMathListTablesTxt,selected=1,container=gridTab4,handler=
   function(h,...) {
     sortval = which(gWidgets2::svalue(gridTab4[1,3])==sortMathListTablesTxt)
     resave_Sorting(4,sortval) #resaving to setupSorting object
     refreshTabLIST2(sort=sortval) #change sorted order
   }
  )  
  
  tabCompLIST2 = gWidgets2::ggroup(container=tabmatchlist2,expand=TRUE,fill=TRUE)
  
  matchL2GUI <- gWidgets2::gtable(items="",container=tabCompLIST2,handler=clickmatchlistQUAN)
  gWidgets2::add(tabCompLIST2,matchL2GUI,expand=TRUE,fill=TRUE)#add to frame
  
  refreshTabLIST2 = function(selInd=NULL,sort=1) { 
    matchlist <- get("resCompLR2",envir=nnTK)
    if(is.null(matchlist)) return() #return if no list found
    if(nrow(matchlist)==0) return() #return if no candidate found
    resTab = casesolver::addRownameTable(matchlist,type=3,L$samplename)
    ord = casesolver::orderTableSort(resTab[,2],resTab[,3],sort) #obtain order of sorting
    
    matchL2GUI[] <- resTab[ord,,drop=FALSE]
    if(nrow(resTab)<nLarge) gWidgets2::size(matchL2GUI) <- list(column.widths=c(30,300,300,70,70,rep(100,ncol(matchlist)-4)))
    if(!is.null(selInd)) gWidgets2::svalue(matchL2GUI) <- selInd #select row
  }
  refreshTabLIST2() #use default sort
  
  ################################################################################################
  #################################Tab 5: Matches: ##############################################
  ################################################################################################
  sortMixListTablesTxt = paste( L$sort.by , c( "#" , L$sort.evid , L$sort.ref ))
  
  gridTab5 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabmixtures),expand=TRUE,fill=TRUE) #kit and population selecter
  gridTab5[1,1] <- gWidgets2::gbutton(text= paste( L$woeperform ),container=gridTab5,handler=f_doWOEcall) #Perform WOE
  gWidgets2::tooltip(gridTab5[1,1]) <- L$tip.matches.woe
  gridTab5[1,2] <- gWidgets2::gbutton(text= paste( L$deconvolvemixtures ) ,container=gridTab5,handler=f_doDCall)# perform DC for all mixtures (with default settings)
  gWidgets2::tooltip(gridTab5[1,2]) <- L$tip.matches.dcall
  
  gridTab5[1,3] <- gWidgets2::gbutton(text= paste( L$show , L$match , L$network ) ,container=gridTab5,handler=f_showMatchNetwork,action="all")  #Button to create matchcloud
  gridTab5[1,4] <- gWidgets2::gbutton(text= paste( L$showfor , L$mixtures )  ,container=gridTab5,handler=f_showMatchNetwork,action="onlymix")  #Button to create matchcloud
  gridTab5[1,5] <- gWidgets2::gbutton(text= paste( L$showfor , L$singlesources ),container=gridTab5,handler=f_showMatchNetwork,action="onlyss")  #Button to create matchcloud
  gridTab5[1,6] <- gWidgets2::gbutton(text= L$export ,container=gridTab5,handler=f_exporttable,action="final")#  
  gWidgets2::tooltip(gridTab5[1,6]) <- L$tip.export 
  
  gridTab5[1,7] <- gWidgets2::gcombobox(items=sortMixListTablesTxt,selected=1,container=gridTab5,expand=TRUE,handler=
   function(h,...) {
     sortval = which(gWidgets2::svalue(gridTab5[1,7])==sortMixListTablesTxt)
     resave_Sorting(5,sortval) #resaving to setupSorting object
     refreshTabLIST(sort=sortval) #change sorted order
   }
  )  
  
  tabmixLIST = gWidgets2::ggroup(container=tabmixtures,expand=TRUE,fill=TRUE)
  mixlistGUI <- gWidgets2::gtable(items="",container=tabmixLIST,handler=clickmixlist)
  gWidgets2::add(tabmixLIST,mixlistGUI,expand=TRUE,fill=TRUE)#add to frame
  
  refreshTabLIST = function(evidsel=NULL,sort=1) {  #get list with all mixtures with info about matched elements 
    #store match-results from comparison (those with LR>threshold) together with all mixtures
    matchstat <- get("mixDataMATCHSTATUS",envir=nnTK) #get match status mixtures
    indUse <- indMixEmpty <- matchstat == "mixture" | matchstat=="" #important for recognizing mixtures or empty
    if(get("setupAdvanced",envir=nnTK)$compSS=="TRUE")  indUse = indUse | TRUE #SHOW ALL #grepl("Unknown ",matchstat)
    mixName <- names(matchstat)[indUse]  #get sampleName of all mixtures
    if( length(mixName)==0 ) return() #return if no mixtures existing   
    nClow <- sapply(get("mixDataLIST",envir=nnTK)[mixName],function(elem) { ceiling(max(sapply(elem,function(x) length(x$adata)))/2) })
    mixlist <- cbind(mixName,"",nClow)
    colnames(mixlist) <- c( L$evidence , L$references , L$numcontr )
  
    #INSERT SINGLE SOURCE PROFILES IDENTICAL MATCHES (MATCH STATUS)    
    matchstat2 =  matchstat[!indMixEmpty] #matchstatus to use
    indInsert = mixlist[,1]%in%names(matchstat2) 
    mixlist[indInsert,2] = matchstat2[match(names(matchstat2),mixlist[indInsert,1])] #Update: insert single source match information
    
    #UPDATE CANDIDATE MATCHES BASED ON COMPARE RESULTS
    matchlist <- get("resMatches",envir=nnTK) #get match list (based on results after comparisons)
    mixlist[match(matchlist[,1],mixName),2] <- matchlist[,2] #Update contributing References
    mixlist[match(matchlist[,1],mixName),3] <- matchlist[,3] #Update number of contributors
    assign("allMixList",mixlist,envir=nnTK)  #store match-results from comparison (those with LR>threshold) together with all mixtures
  
    resTab = casesolver::addRownameTable(mixlist,type=3,L$samplename)
    ord = casesolver::orderTableSort(resTab[,2],resTab[,3],sort) #obtain order of sorting
    mixlistGUI[] <- resTab[ord,,drop=FALSE]
  
    if(nrow(resTab)<nLarge) gWidgets2::size(mixlistGUI) <- list(column.widths=c(30,300,500,100))
    if(!is.null(evidsel)) gWidgets2::svalue(mixlistGUI) = which(mixlist[,1]==evidsel)  #marking line
  }
  refreshTabLIST() #use default sort
  
  ################################################################################################
  #################################Tab 6: Deconvoluted: ##############################################
  ################################################################################################
  sortDCTableTxt <- paste( L$sort.by , c( "#" , L$sort.evid , L$similarity ))
  
  gridTab6 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabdeconv)) #kit and population selecter
  gridTab6[1,1] <- gWidgets2::gbutton(text=paste( L$export , L$table ),container=gridTab6,handler=f_exporttable,action="dc")# 
  gWidgets2::tooltip(gridTab6[1,1]) <- L$tip.export
  
  #Export selected
  gridTab6[1,2] <- gWidgets2::gbutton(text=paste( L$export , L$selected ),container=gridTab6,handler=function(h,...) {
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #none found
   selID =  as.integer(gsub("#","",gWidgets2::svalue(DClistGUI) ) )
   if(length(selID)==0) return() #none selected
   dupInd = which(duplicated(colnames(DClist))) #get index of duplicated
   
   refn = DClist[selID,1] #get refnames
   refTab = t(DClist[selID,-c(1:4,dupInd)])
   rownames(refTab) = refn 
   saveTable(sampleListToTable(casesolver::tabToListRef(tab=refTab,setEmpty=FALSE))) #save to file
  })# End export selected
  gWidgets2::tooltip(gridTab6[1,2]) <- L$tip.dc.exportsel
  
  #Delete selected
  gridTab6[1,3] <- gWidgets2::gbutton(text=paste( L$delete , L$selected ),container=gridTab6,handler=function(h,...) {
     DClist <- get("DClist",envir=nnTK)
     if(is.null(DClist)) return() #none found
     selID =  as.integer(gsub("#","",gWidgets2::svalue(DClistGUI) ) )
     if(length(selID)==0) return() #none selected
     DClist <- DClist[-selID,,drop=FALSE] #update
     assign("DClist",DClist,envir=nnTK)
     refreshDCLIST()
  })
  gWidgets2::tooltip(gridTab6[1,3]) <- L$tip.dc.delete
  
  #Show ratio-to-next
  gridTab6[1,4] <- gWidgets2::gbutton(text=paste( L$show , L$probRatioToNext ),container=gridTab6,handler=function(h,...) {
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   selID =  as.integer(gsub("#","",gWidgets2::svalue(DClistGUI) ) )
   if(length(selID)==0) return() #none selected
   dupInd = which(duplicated(colnames(DClist))) #get index of duplicated
   tab = cbind(DClist[selID,-c(1:4,dupInd)], DClist[selID,dupInd])
   colnames(tab) = c(L$Genotype, L$probRatioToNext)
   
   showGDFtable(L$probRatioToNext, casesolver::addRownameTable(tab,samplename = L$Marker))
  })
  gWidgets2::tooltip(gridTab6[1,4]) <- L$tip.dc.showprobs
  
  tabDCLIST = gWidgets2::ggroup(container=tabdeconv,expand=TRUE,fill=TRUE)
  DClistGUI <- gWidgets2::gtable(items="",multiple = FALSE,container=tabDCLIST,handler=addDCprofile)
  gWidgets2::add(tabDCLIST,DClistGUI,expand=TRUE,fill=TRUE)#add to frame
  
  #Show PHexp for deconvolved profile (including other fit):
  gridTab6[1,5] <- gWidgets2::gbutton(text=paste( L$expPHplot ),container=gridTab6,handler=function(h) {
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   selID =  as.integer(gsub("#","",gWidgets2::svalue(DClistGUI) ) )
   if(length(selID)==0) return() #none selected
   
   #Prepare profiles:
   candname = DClist[selID,1] #name of candidate
   splitPtrn = "-C" #split pattern to remove last component
   evids =  strsplit(candname,splitPtrn)[[1]] #remove last component
   candname = paste0("C",evids[length(evids)]) #obtain name of component
   evids = paste(evids[1:(length(evids)-1)],collapse=splitPtrn) #keep only evidence name but make it possible for evidence names to have pattern
   evids = strsplit(evids,"/")[[1]] #obtain evidences to consider (possibly replicates)
   conds = strsplit(DClist[selID,2],"/")[[1]] #obtain conditional refs
   NOC = as.integer(DClist[selID,3])
  
   #Prepare profile of candidate:
   dupInd = which(duplicated(colnames(DClist))) #get index of duplicated column names
   candProfile =DClist[selID,-c(1:4,dupInd),drop=FALSE] #get profile of candidate
   rownames(candProfile) = candname
   
   addedProfiles = tabToListRef(candProfile)
   doDC(nC=NOC,evids=evids,refs=conds,showPlot=TRUE,useplotly=TRUE,addedProfiles=addedProfiles)
  })
  
  #add button for adding ref
  gridTab6[1,6] <- gWidgets2::gbutton(text=paste( L$addasref ),container=gridTab6,handler=addDCprofile)
  gWidgets2::tooltip(gridTab6[1,6]) <- L$tip.dc.addasref
  
  #Possible to sort table
  gridTab6[1,7] <- gWidgets2::gcombobox(items=sortDCTableTxt,selected=1,container=gridTab6,expand=TRUE,handler=
                                         function(h,...) {
                                           sortval = which(gWidgets2::svalue(gridTab6[1,7])==sortDCTableTxt)
                                           #resave_Sorting(6,sortval) #DONT store sortvalue to setupSorting object
                                           refreshDCLIST(sort=sortval) #change sorted order
                                         })  
  
  #Function to update DC-table GUI 
  refreshDCLIST = function(sort=1) {  #get list of matched elements
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   dupInd = which(duplicated(colnames(DClist))) #get index of duplicated column names
  #   locnames = colnames(DClist)[dupInd]  #column names
   
   #ORDER ROWS WRT criteria
   ord = 1:nrow(DClist) #default is no sorting
   if(sort==2) { #sort based on evid name
     ord = order(DClist[,1])
  #Order by clustering:
   } else if(nrow(DClist)>=3 && sort==3) { #perform sorting based on similarity (using clustering), must have at least 3 rows
     
     #Idea: Order samples based on hierachical clustering structure (based on distance outcome)
     ord = casesolver::clusterOrder( table=DClist[,-c(1:4,dupInd),drop=FALSE] )
   }# end if clustering
  
   DCtab =  casesolver::addRownameTable(DClist[,-dupInd,drop=FALSE] ,type=3,L$samplename) #add # to first column 
   
   #PREPARE OUT-PUT TABLE (SORTING + TRUNCATE COLUMNS):
   DClistGUI[] <-  DCtab[ord,,drop=FALSE] #DC table to consider 
   gWidgets2::size(DClistGUI) <- list(column.widths=c(30,100,100,30,50,rep(50,ncol(DClistGUI)-5))) 
  }
  refreshDCLIST()
  
  ################################################################################################
  #################################Tab 7: Weight of evidence: ####################################
  ################################################################################################
  gridTab7 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe( L$further ,container=tabWOE)) #kit and population selecter
  gridTab7[1,1] <- gWidgets2::gbutton(text=paste( L$export , L$table ),container=gridTab7,handler=f_exporttable,action="woe")#  
  gridTab7[1,2] <- gWidgets2::gbutton(text=paste( L$show , L$lrpermarker  ),container=gridTab7,handler=showLRperMarker)  #show LR-marker for selected
  
  #SHOW PARAM ESIMTATES
  gridTab7[1,3] <- gWidgets2::gbutton(text=paste( L$show , "param est."  ),container=gridTab7, handler=function(h) {
   resWOEeval <- get("resWOEeval",envir=nnTK) #obtain WOE results
   if(is.null(resWOEeval)) return()
   ids <- as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI)))
   for(id in ids) {
     hypsettxt = paste0("Hyp. set #",id)
     samplename = resWOEeval$resTable[id,1]
       
     #Show mixture proportions in a separe plot:
     if(require(plotrix)) {
       MxHp = sort(resWOEeval[[id]]$MxRefs$hp,decreasing = TRUE) #sort to make same colors
       MxHd = sort(resWOEeval[[id]]$MxRefs$hd,decreasing = TRUE) #sort to make same colors
       s0 = 2
       labelHp = paste0(names(MxHp)," ",round(MxHp*100),"%") #Obtain labels
       labelHd = paste0(names(MxHd)," ",round(MxHd*100),"%") #Obtain labels
       plotfn = paste0("WOEpie",id) #file name
       graphics.off() #remove other plots first
       png(plotfn,height=400,width=1000)
       par(mfrow=c(1,2))
       marg = c(0,5,0,6)
       plotrix::pie3D(MxHp,radius=0.9,labels=labelHp,explode=0.1,main="Hp",mar=marg)
       plotrix::pie3D(MxHd,radius=0.9,labels=labelHd,explode=0.1,main="Hd",mar=marg)
       mtext(samplename,cex=1.5,outer=TRUE,line=-2)
       dev.off()
       
       TopWin <- gWidgets2::gwindow(paste0("Mixture proportions for ",hypsettxt))#,width=800,height=400)
       gWidgets2::gimage(plotfn,container=TopWin)
       gWidgets2::focus(TopWin) = TRUE
       file.remove(plotfn) #remove file after
     }
     
    #show params in table   
     tab = resWOEeval[[id]]$paramTable
     showGDFtable(hypsettxt,casesolver::addRownameTable(tab,samplename = "Param"))
   }
  })
  
  #SHOW MODEL VALIDATION
  gridTab7[1,4] <- gWidgets2::gbutton(text=paste( L$show , "model validation"  ),container=gridTab7, handler=function(h) {
    resWOEeval <- get("resWOEeval",envir=nnTK) #obtain WOE results
    if(is.null(resWOEeval)) return()
    
    ids <- as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI)))
    for(id in ids) {
      hypsettxt = paste0("Hyp. set #",id)
      print(paste0("Calculating for ", hypsettxt))
      kit0 = casesolver::getEnvirKit(nnTK)
      plotfn1 =  paste0("WOEvalidHp",id) #file name
      plotfn2 =  paste0("WOEvalidHd",id) #file name
      graphics.off() #remove other plots first
      size0 = 600
      png(plotfn1,height=size0,width=size0)
      euroformix::validMLEmodel( resWOEeval[[id]]$mleHp, kit=kit0,plottitle = "Hp")
      dev.off() 
      
      png(plotfn2,height=size0,width=size0)
      euroformix::validMLEmodel( resWOEeval[[id]]$mleHd, kit=kit0,plottitle = "Hd" )
      dev.off() 
      
      TopWin <- gWidgets2::gwindow(paste0("Model validation for ",hypsettxt))#,width=800,height=400)
      ggrp <- gWidgets2::ggroup(container=TopWin)
      gWidgets2::gimage(plotfn1,container=ggrp)
      gWidgets2::gimage(plotfn2,container=ggrp)
      gWidgets2::focus(TopWin) = TRUE
      file.remove(plotfn1) #remove file after
      file.remove(plotfn2) #remove file after
    }    
  })  
  
  #SHOW MODEL FITTED PH
  gridTab7[1,5] <- gWidgets2::gbutton(text=paste( L$show , "model fitted PH"  ),container=gridTab7, handler=function(h) {
   resWOEeval <- get("resWOEeval",envir=nnTK) #obtain WOE results
   if(is.null(resWOEeval)) return()
   id <- as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI)))
   if(length(id)>1) {
     gWidgets2::gmessage("Only one hypothesis must be selected for modifying a statement. Please re-select.","Wrong selection", icon="info")
     return()  
   }
   resWOEeval_id = resWOEeval[[id]]
   if(is.null(resWOEeval_id)) return() #couldnt find object
   
   hyp = gWidgets2::gconfirm("Select hypothesis (Yes=Hp, No=Hd)","Choose hypothesis",icon = "question")
   kitname = casesolver::getEnvirKit(nnTK)
   type = casesolver::getSampleType2( resWOEeval_id$mleHp$model$samples ,kitname) #get sample type
   if(hyp) {
     makePlotTop(type,resWOEeval_id$mleHp,NULL,kitname)
   } else {
     makePlotTop(type,resWOEeval_id$mleHd,NULL,kitname)
   }
  })  
  
  #SHOW STATEMENTS
  gridTab7[1,6] <- gWidgets2::gbutton(text=paste( L$show, L$statement ),container=gridTab7, handler=function(h) {
   resWOEeval <- get("resWOEeval",envir=nnTK) #obtain WOE results
   if(is.null(resWOEeval)) return()
   id <- as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI)))
   if(length(id)>1) {
     gWidgets2::gmessage("Only one hypothesis must be selected for modifying a statement. Please re-select.","Wrong selection", icon="info")
     return()  
   }
   state = resWOEeval[[id]]$statement #obtain statement
   env = new.env( parent = emptyenv() ) 
   assign("text",state,envir=env)
   casesolver::textEditor(env) #use inbuilt text editor to modify statement
   resWOEeval[[id]]$statement = get("text",envir=env)#obtain possibly edited text
   assign("resWOEeval",resWOEeval,envir=nnTK) #store text to environment
  }) 
  
  #CONSERVATIVE CALCULATION:
  gridTab7[1,7] <- gWidgets2::gbutton(text=paste( L$calccons ),container=gridTab7,handler=function(h) {
   resWOEeval <- get("resWOEeval",envir=nnTK) 
   if(is.null(resWOEeval)) return()
   ids =  as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI) ) )
   
   #(RE)CULCULATE CONS LR FOR SELECTED HYPS
   for(id in ids) {
     print(paste0("Performing conservative calculations for hypothesis set #",id))
     mcmcOpt = get("setupMCMC",envir=nnTK)
     
     mcmcObjList=NULL
     if(!is.null(resWOEeval[[id]]$mcmc)) mcmcObjList = resWOEeval[[id]]$mcmc$mcmcObjList #obtain previous calcs
     resWOEeval[[id]]$mcmc = euroformix::calcLRmcmc(resWOEeval[[id]]$mleHp,resWOEeval[[id]]$mleHd, mcmcOpt$niter,mcmcOpt$delta,mcmcOpt$quantile,mcmcOpt$seed, verbose=TRUE,traceplot=TRUE, mcmcObjList=mcmcObjList)
     state <- resWOEeval[[id]]$statementNoWOE #obtain statement without woe value
     consLR <- LRuse <- resWOEeval[[id]]$mcmc$log10LRcons
     bayesLR <- resWOEeval[[id]]$mcmc$log10LRbayes
     
     #insert woe strength (also possibly as text)
     LRtxt = signif(10^LRuse,2)
     useVerbalLR =  get("setupReportOpt",envir=nnTK)["verbalLR"]  #If using verbal LR
     if(useVerbalLR) {
       LRverbal = casesolver::number2word(LRuse,L) #obtain verbal LR
       LRtxt = paste0(LRverbal[1]," (",LRverbal[2],")") #add verbal in additon to down-rounded number
     }
     state = gsub("$LRtxt",LRtxt,state,fixed=TRUE)
  
     #INSERT BACK TO OBJECT:
     s0 = 2
     resWOEeval[[id]]$statement <- state
     resWOEeval[[id]]$bayesLR <- resWOEeval$resTable[id,6] <- round(bayesLR,s0) #insert updated values
     resWOEeval[[id]]$consLR  <- resWOEeval$resTable[id,7] <- round(consLR,s0) #insert updated value
   }
   
   assign("resWOEeval",resWOEeval,envir=nnTK) #store object to environment again
   refreshWOELIST() #refresh table
  })
  
  #DELETE SELECTED WOE CALCULATIONS
  gridTab7[1,8] <- gWidgets2::gbutton(text=paste( L$delete , L$selected ),container=gridTab7,handler=function(h) {
   resWOEeval <- get("resWOEeval",envir=nnTK) 
   if(is.null(resWOEeval)) return()
   nHypSets = nrow(resWOEeval$resTable) #number of hypotheses
   if(nHypSets==0) return() #return if none to select (can't happen)
   ids =  as.integer(gsub("#","",gWidgets2::svalue(WOElistGUI) ) ) #selected row
  
   existingIDs = seq_len(nHypSets)
   if(all(existingIDs%in%ids)) { #IF ALL SELECTED
     WOElistGUI[] <- resWOEeval <- NULL      #removing ALL element (set everything as zero)
   } else {
     keepIDs = setdiff(existingIDs,ids) #get which rows to keep
     tmpTable = resWOEeval$resTable[keepIDs,,drop=FALSE] #obtain subset of restable
     resWOEeval= resWOEeval[keepIDs]
     resWOEeval$resTable = tmpTable #insert updated table
   }
   assign("resWOEeval",resWOEeval,envir=nnTK) #store object to environment again
   refreshWOELIST() #refresh table
  })
  
  #INITITATE WOE TABLE
  tabWOELIST = gWidgets2::ggroup(container=tabWOE,expand=TRUE,fill=TRUE)
  WOElistGUI <- gWidgets2::gtable(items="",multiple = TRUE, container=tabWOELIST, handler=showLRperMarker)
  gWidgets2::add(tabWOELIST,WOElistGUI,expand=TRUE,fill=TRUE)#add to frame
  
  refreshWOELIST = function() {  #get list of matched elements
   resTable <- get("resWOEeval",envir=nnTK)$resTable #obtain table
   if(is.null(resTable)) return() #return if no list found
   WOElistGUI[] <- casesolver::addRownameTable( resTable ,type=3) #add DC-list
   gWidgets2::size(WOElistGUI) <- list(column.widths=c(30,150,75,150,rep(50,ncol(WOElistGUI)-4))) 
  }
  refreshWOELIST()
  
  #Add tooltips
  gWidgets2::tooltip(gridTab7[1,1]) <- L$tip.export
  gWidgets2::tooltip(gridTab7[1,2]) <- L$tip.woe.lrpermarker
  gWidgets2::tooltip(gridTab7[1,3]) <- L$tip.woe.paramest
  gWidgets2::tooltip(gridTab7[1,4]) <- L$tip.woe.modelvalid
  gWidgets2::tooltip(gridTab7[1,5]) <- L$tip.woe.fittedPH
  gWidgets2::tooltip(gridTab7[1,6]) <- L$tip.woe.statement
  gWidgets2::tooltip(gridTab7[1,7]) <- L$tip.woe.cons
  gWidgets2::tooltip(gridTab7[1,8]) <- L$tip.woe.delete
  
  ################################################################################################
  #Running through all windows to avoid "bleed through"
  sapply(7:1, function(x) {
   gWidgets2::svalue(nb) <- x
  })
  gWidgets2::visible(mainwin) <- TRUE
  gWidgets2::size(mainwin) = c(mwW,mwH) #set width and hight after window is open (must)
  setFocus()

} #end main function
