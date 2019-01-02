#' @title gui
#' @description gui (GUI) is a GUI wrapper for CaseSolver
#' @details The function starts the CaseSolver software.
#' @param envirfile A file to a saved environment of a project
#' @export


gui = function(envirfile=NULL) {
 LUSsymbol <- "_"
 #size of main window
 mwH <- 500
 mwW <- 1000
 w0 <- 20 #default textwidth (for gWidgets2::gedit)
 nLarge = 10000 #this is large threshold of number of references

 #type of gwidgets-kit
 library(gWidgets2) #install.packages("gWidgets2")
 library(gWidgets2tcltk) #install.packages("gWidgets2tcltk")
 options(guiToolkit="tcltk")

 #models:
 library(euroformix) #requires 2.0.2 
 reqver <- "2.0.2" #required version of euroformix
 if( packageVersion("euroformix")<reqver ) stop(paste0("This version of CaseSolver requires euroformix ",reqver," or newer!"))
 version =  packageVersion("casesolver") #follows same as rpackage
 #sessionInfo()

 #version:
 version =  packageVersion("casesolver") #follows same as rpackage

 #v0.1: Initial version
 #v0.2: Tables are fixed
 #v0.3: Reports are added
 #v0.4: Calculate RMP added	
 #v0.5: Deconvolution added
 #v0.7: Supporting MPS kit (requires >=euroformix_1.11.0)
 #v0.8: Matches either based on QualLR, QuanLR  or both 
 #v1.0: Community Release
 #v1.1: Small changes (gui,model)
 #v1.2: Added QuanLR and DC functionalities

 #software name:
 softname <- paste0("CaseSolver v",version)

 spc <- 10  #Spacing between widgets

 #REPORT ITEMS (may change over versions)
 reportitems <- paste0("Show ", c("Header","References","Single sources (alleles)","Mixtures (alleles)","Consensus (alleles)","Single sources w/PH","Mixtures w/PH","Metadata","Comparison matrix","MatchList (Qual)","MatchList (Quan)","FinalMatchList","MatchNetwork","RMNE for evidence","RMP for references","IBS for references","Settings","EPG figures for single sources","EPG figures for mixtures","MatchStatus","Expected PH plots"))
# reportitems <- paste0("Show ", c("Header","References","Single sources (alleles)","Mixtures (alleles)","Consensus (alleles)","Single sources w/PH","Mixtures w/PH","Metadata","Comparison matrix","MatchList (Qual)","MatchList (Quan)","FinalMatchList","MatchNetwork","RMNE for evidence","RMP for references","IBS for references","Settings","EPG figures for single sources","EPG figures for mixtures","MatchStatus","Estimated profiles"))

#########################################################
################Environment variables####################
#########################################################

#NB: pgkPath and .sep must be changed before compiling R-package!
 pgkPath <- path.package("casesolver", quiet = FALSE) # Get package path.
 #pgkPath <- "config"# Get package path (possible to synchronize)
 .sep <- .Platform$file.sep # Platform dependent path separator. 

 setupFileThresh <- paste(pgkPath,"configThresh",sep=.sep)  
 setupFileModel <- paste(pgkPath,"configModel",sep=.sep)  
 setupFileKit <- paste(pgkPath,"configKit",sep=.sep)  
 setupFilePop <- paste(pgkPath,"configPop",sep=.sep)  
 setupFileCase <- paste(pgkPath,"configCase",sep=.sep)  
 setupFileImport <- paste(pgkPath,"configImport",sep=.sep)  
 setupFileExport <- paste(pgkPath,"mmTK.Rdata",sep=.sep)  #used to run euroformix
 setupFileExport2 <- paste(pgkPath,"mmTK2.Rdata",sep=.sep)  #used to run euroformix
 setupRead = function(file) scan(file,what=character(),quiet=TRUE,sep="\n") #skip line is new element
 setupWrite = function(vals,file) write(vals,file=file,sep="\n")   #save to file in installation folder
 setupFileReport <- paste(pgkPath,"configReport",sep=.sep) #used to set report layout
 setupFileAdvanced <- paste(pgkPath,"configAdvanced",sep=.sep) #used to set report layout
 setupFileView <- paste(pgkPath,"configView",sep=.sep) #used to set GUI layout

 #the files are stored in system settings:
 if(!file.exists(setupFileThresh)) {  #use default values if not existing
   setupThresh =  list(MACthresh=0.9,LRthresh1=10,LRthresh2=1000,minLociSS=7,minIBS=20,ratio=10,probA=0.99)
 } else {
   optF <- setupRead(file=setupFileThresh)
   setupThresh = list(MACthresh=as.numeric(optF[1]),LRthresh1=as.numeric(optF[2]),LRthresh2=as.numeric(optF[3]),minLociSS=as.integer(optF[4]),minIBS=as.integer(optF[5]),ratio=as.numeric(optF[6]),probA=as.numeric(optF[7]))
 }
 if(!file.exists(setupFileModel)) {  #use default values if not existing
   setupModel =  list(threshT=50,dropinC=0.05,dropinL=0.01,degrad="ON",stutt="OFF",modeltype=3) #threshT: default value of detection threshold value,dropinC: Dropin probability parameter in model,dropinL: Dropin peak height distribution parameter Lambda. Modeltype={1,2,3}={qual,quan,both}
 } else {
   optF <- setupRead(file=setupFileModel)
   setupModel= list(threshT=as.numeric(optF[1]),dropinC=as.numeric(optF[2]),dropinL=as.numeric(optF[3]),degrad=optF[4],stutt=optF[5],modeltype=as.numeric(optF[6]))
 }
 if(!file.exists(setupFileKit)) {  #use default values if not existing
   setupKit = list(kitname="") #empty kit by default. Select by short names (from getKit function)
 } else {
   optF <- setupRead(file=setupFileKit)
   setupKit = list(kitname=optF[1])
 }
 if(!file.exists(setupFilePop)) {  #use default values if not existing
   setupPop = list(popfile="",amel="") #empty population by default. Points to a text-file.
 } else {
   optF <- setupRead(file=setupFilePop)
   setupPop = list(popfile=optF[1],amel=optF[2])
 }
 if(!file.exists(setupFileCase)) {  #use default values if not existing
   setupCase = list(casepath="") #empty casepath by default. Points to a folder with CASEID given as folder names. 
 } else {
   optF <- setupRead(file=setupFileCase)
   setupCase = list(casepath=optF[1])
 }
 if(!file.exists(setupFileImport)) {  #use default values if not existing
   setupImport = list(importfile="") #empty importfile by default (must contain the R-function importData. Points to a text-file.
 } else {
   optF <- setupRead(file=setupFileImport)
   setupImport = list(importfile=optF[1])
 }
 if(!file.exists(setupFileReport)) {  #use default values if not existing
   setupReport = list(checked=rep(TRUE,length(reportitems))) #Gives boolean of how layout should be given
 } else {
   optF <- setupRead(file=setupFileReport)
   setupReport = list(checked=optF=="TRUE")
 }
 if(!file.exists(setupFileAdvanced)) {  #use default values if not existing
   setupAdvanced = list(maxC1=6,maxC2=3,nDone=4,pD0=0.1)
 } else {
   optF <- setupRead(file=setupFileAdvanced)
   setupAdvanced = list(maxC1=as.integer(optF[1]),maxC2=as.integer(optF[2]),nDone=as.integer(optF[3]),pD0=as.numeric(optF[4]))
 }
 if(!file.exists(setupFileView)) {  #use default values if not existing
   setupView = list(importHorizontal="FALSE",matchmatrix="FALSE")
 } else {
   optF <- setupRead(file=setupFileView)
   setupView = list(importHorizontal=optF[1],matchmatrix=optF[2])
 }

 if(is.null(envirfile)) {
  nnTK = new.env() #create new envornment object

  #Toolbar options: can be changed any time by using toolbar
  assign("setupThresh",setupThresh,envir=nnTK) 
  assign("setupModel",setupModel,envir=nnTK) 
  assign("setupKit",setupKit,envir=nnTK) 
  assign("setupPop",setupPop,envir=nnTK) 
  assign("setupCase",setupCase,envir=nnTK) 
  assign("setupImport",setupImport,envir=nnTK) 
  assign("setupReport",setupReport,envir=nnTK) 
  assign("setupAdvanced",setupAdvanced,envir=nnTK) 
  assign("setupView",setupView,envir=nnTK) 

  #initializing environment variables
  assign("workdir",NULL,envir=nnTK) #assign working directory to nnTK-environment

  #Considered CaseID:
  assign("caseID",NULL,envir=nnTK) #Used for report

  #imported data:
  assign("popFreq",NULL,envir=nnTK) #assign to nnTK-environment
  assign("mixDataTABLE",NULL,envir=nnTK) #Table of evidence profiles (only alleles)
  assign("refDataTABLE",NULL,envir=nnTK) #Table of ref profiles (only alleles)
  assign("mixDataMATCHSTATUS",NULL,envir=nnTK) #vector with match status for evidence profiles
  assign("mixDataLIST",NULL,envir=nnTK)  #list of evidence profiles (alleles and heights)
  #assign("refDataLIST",NULL,envir=nnTK) #list of evidence profiles (alleles)
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
  assign("DClistReport",NULL,envir=nnTK)  #store list with confirmed deconvoluted reference candidates 
  assign("resRMP",NULL,envir=nnTK)  #store random match prob results 
  assign("resIBS",NULL,envir=nnTK)  #store IBS compare results (not used)
  assign("storedFitHp",NULL,envir=nnTK) #store mle-fitted objects under Hp (for all match candidates)
 } else {
  load(envirfile) #loading environment
 }

#################################################################################
###########################GUI HELPFUNCTIONS#####################################
#################################################################################

 #This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
 mygfile <- function(text,type,encoding="latin1",filter=list(),initf=NULL) {
   file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
   Encoding(file ) <- encoding
   return(file)
 }

 saveTable = function(tab,sep="txt") {
  tabfile  <- mygfile(text="Save table",type="save") #csv is correct format!
  if(!is.na(tabfile)) {
   if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
   if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
   if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=",",row.names=FALSE) 
   print(paste("Table saved in ",tabfile,sep=""))
  }
 } #end file

 setPopFreq = function(change=FALSE) { #helpfunction to read popFreq from file and set to environment
     if(!change && !is.null( get("popFreq",envir=nnTK))) return(TRUE) #return if already set

  	opt <- get("setupPop",envir=nnTK) 
     tryCatch( {
      popFreq <- getFreqs(opt$popfile)
      AMEL <- c(0.75,0.25)
      names(AMEL) <- c("X","Y")
      if( opt$amel=="TRUE" && !any(grepl("AM",names(popFreq))) ) popFreq$AMEL =AMEL
      assign("popFreq",popFreq,envir=nnTK) #assign popFreq to nnTK-environment
      #if(verbose) print(popFreq) #print first time
     }, error = function(e) return(FALSE) )
     if(is.null(get("popFreq",envir=nnTK))) { #if popFreq is missing
       gWidgets2::gmessage("Please select a valid file with population frequencies for further analysis.\nGo to Setup->Population Frequencies") 
       return(FALSE) 
     }
	return(TRUE)
 } 

 canPrintEPG = function() { #Helpfunction to check if can print EPG (kit specified)
   kit0 <- get("setupKit",envir=nnTK)$kitname
   return( !is.na(kit0) && !is.null(kit0) && kit0!="")  #possible to print EPG?
 }

###########################FILE#####################################
 f_setwd = function(h,...) {
  dirfile = mygfile(text="Select folder",type="selectdir")
  if(length(dirfile)>0) {
   setwd(dirfile)
   assign("workdir",dirfile,envir=nnTK) #assign working directory
  }
 }
 f_openproj = function(h,...) {
  projfile = mygfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata"))))
  if(!is.na(projfile)) {
   gWidgets2::dispose(mainwin)
   gui(projfile) #send environment into program
  }
 }
 f_saveproj = function(h,...) {
  projfile = mygfile(text="Save project",type="save")
  if(!is.na(projfile)) {
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata")
   print("Size of stored objects (in MB):") #prints size of each stored object
   print(sapply(nnTK,object.size)/1e6) #prints size of each stored object
   save(nnTK,file=projfile,compress="xz",eval.promises=FALSE,precheck=FALSE,compression_level=2)
   print(paste("Project saved in ",projfile,sep=""))
  }
 }
 f_quitproj = function(h,...) {
  ubool <- gWidgets2::gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(ubool) {
    f_saveproj(h)
  } else { 
   print("Program terminated without saving")
  }
  gWidgets2::dispose(mainwin) #remove window!
 }

###########################REPORT SET#####################################
f_reportlay = function(h,...) { #GUI function to set report layout
   opt <- get("setupReport",envir=nnTK) #get layout settings
   if(is.null(opt)) {
    checkv <- rep(TRUE,length(reportitems))
   } else {
    checkv  <- opt$checked
   }
   setwin <- gWidgets2::gwindow(paste0("Report layout"),visible=FALSE)
   tabval = gWidgets2::ggroup(spacing=0,container=(setwin),horizontal=FALSE) 
   checkbox <- gWidgets2::gcheckboxgroup(reportitems,checked=checkv ,use.table=FALSE,container=tabval)
   gWidgets2::add(tabval,checkbox)#,expand=TRUE,fill=TRUE)
   savebut <- gWidgets2::gbutton("Save",checked=checkv ,use.table=FALSE,container=tabval,handler = function(h, ...) { 
    opt$checked <- reportitems%in%gWidgets2::svalue(checkbox)
    assign("setupReport",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileReport)    #save to file in installation folder
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
}

###########################SETTINGS#####################################
 f_threshsel=  function(h,...) { #GUI function to set threshold values
   opt <- get("setupThresh",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Threshold settings"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="MAC threshold (Comparison)",container=tabval)
   tabval[1,2] <- gWidgets2::gedit(text=opt$MACthresh,width=w0,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Qual. LR threshold (Comparison)",container=tabval)
   tabval[2,2] <- gWidgets2::gedit(text=opt$LRthresh1,width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text="Quan. LR threshold (Comparison)",container=tabval)
   tabval[3,2] <- gWidgets2::gedit(text=opt$LRthresh2,width=w0,container=tabval)
   tabval[4,1] <- gWidgets2::glabel(text="Minimum loci for being SS match (Import)",container=tabval)
   tabval[4,2] <- gWidgets2::gedit(text=opt$minLociSS ,width=w0,container=tabval)
   tabval[5,1] <- gWidgets2::glabel(text="Minimum IBS for being relative candidate (IBS)",container=tabval)
   tabval[5,2] <- gWidgets2::gedit(text=opt$minIBS ,width=w0,container=tabval)
   tabval[6,1] <- gWidgets2::glabel(text="Prob-ratio to next (Deconvolution)",container=tabval)
   tabval[6,2] <- gWidgets2::gedit(text=opt$ratio ,width=w0,container=tabval)
   tabval[7,1] <- gWidgets2::glabel(text="Prob. single allele (Deconvolution)",container=tabval)
   tabval[7,2] <- gWidgets2::gedit(text=opt$probA ,width=w0,container=tabval)
   tabval[8,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 
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
   } )
   gWidgets2::visible(setwin) <- TRUE
 }
 f_modelsel=  function(h,...) {
   radiotxt = c("ON","OFF")
   modtypetxt <- c("Qualitative (LRmix)","Quantitative (EuroForMix)","Both")
   opt <- get("setupModel",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Model settings"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 

   tabval[1,1] <- gWidgets2::glabel(text="Model type(s)",container=tabval)
   tabval[1,2] <- gWidgets2::gcombobox(items=modtypetxt,selected=opt$modeltype,horizontal=TRUE,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Detection threshold (EFM)",container=tabval)
   tabval[2,2] <- gWidgets2::gedit(text=opt$threshT,width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text="Dropin probability",container=tabval)
   tabval[3,2] <- gWidgets2::gedit(text=opt$dropinC,width=w0,container=tabval)
   tabval[4,1] <- gWidgets2::glabel(text="Dropin peak height Lambda (EFM)",container=tabval)
   tabval[4,2] <- gWidgets2::gedit(text=opt$dropinL,width=w0,container=tabval)
   tabval[5,1] <- gWidgets2::glabel(text="Degradation model (EFM)",container=tabval)
   tabval[5,2] <- gWidgets2::gradio(items=radiotxt,selected=which(opt$degrad==radiotxt),horizontal=TRUE,container=tabval)
   tabval[6,1] <- gWidgets2::glabel(text="Stutter model (EFM)",container=tabval)
   tabval[6,2] <- gWidgets2::gradio(items=radiotxt,selected=which(opt$stutt==radiotxt),horizontal=TRUE,container=tabval)
   tabval[7,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 
    opt$modeltype <- which(gWidgets2::svalue(tabval[1,2])==modtypetxt)
    opt$threshT <- as.numeric(gWidgets2::svalue(tabval[2,2]))
    opt$dropinC <- as.numeric(gWidgets2::svalue(tabval[3,2])) 
    opt$dropinL <- as.numeric(gWidgets2::svalue(tabval[4,2]))
    opt$degrad <- gWidgets2::svalue(tabval[5,2])
    opt$stutt <- gWidgets2::svalue(tabval[6,2])
    assign("setupModel",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileModel)    #save to file in installation folder
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
 }

 f_kitsel=  function(h,...) {
   opt <- get("setupKit",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Select Kit"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="Selected kit:",container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$kitname,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Select kit:",container=tabval)
   tabval[2,2] <- gWidgets2::gcombobox(items=euroformix::getKit(), width=100, selected=0, editable = FALSE, container = tabval)
   tabval[3,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 
    opt$kitname <- svalue(tabval[2,2])
    assign("setupKit",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileKit)    #save to file in installation folder
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
 }
 f_popsel=  function(h,...) {
   opt <- get("setupPop",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Select population frequency file"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="Selected frequency file:",container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$popfile,container=tabval)
   tabval[2,1] <- gWidgets2::gcheckbox(text="Include AMEL",checked = ifelse(opt$amel=="TRUE",TRUE,FALSE),container=tabval)
   tabval[3,1] <- gWidgets2::gbutton("Select file", container = tabval,handler = function(h, ...) { 
    ff <- mygfile("Select file",type="open")
    if(length(ff)==0) return()
    opt$popfile <- ff
    opt$amel <- svalue(tabval[2,1])
    assign("setupPop",opt,envir=nnTK)  #assign user-value to opt-list
    ok <- setPopFreq(change=TRUE)	#Assume that another freq file has been selected
    if(ok) {
      setupWrite(unlist(opt),file=setupFilePop)    #save to file in installation folder if successful
      gWidgets2::dispose(setwin) #remove subwindow
      #f_popsel(); #update gui window again after selecting new folder
	}
   })
   tabval[4,1] <- gWidgets2::gbutton("View frequencies", container = tabval,handler = function(h, ...) { 
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
    colnames(tab) <- paste0("A",1:ncol(tab))
    setwin2 <- gWidgets2::gwindow(paste0("Allele frequencies"),visible=FALSE) 
    guitab <- gWidgets2::gdf(items=tab,container = setwin2) 
    gWidgets2::add(setwin2,guitab,expand=TRUE,fill=TRUE) 
    gWidgets2::visible(setwin2) <- TRUE
   })
   gWidgets2::visible(setwin) <- TRUE
 }#end function
 f_casedirsel =  function(h,...) {
   opt <- get("setupCase",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Select Casefolder Path"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="Selected path with CASE folders:",container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$casepath,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Select path:",container=tabval)
   tabval[2,2] <- gWidgets2::gbutton("Set folder", container = tabval,handler = function(h, ...) { 
    opt$casepath <- mygfile("Select folder",type="selectdir")
    assign("setupCase",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileCase)    #save to file in installation folder
    gWidgets2::dispose(setwin) #remove subwindow
    gWidgets2::dispose(mainwin) #remove main window
    gui(); #update gui window again after selecting new folder
   })
  gWidgets2::visible(setwin) <- TRUE
 }
 f_importsel=  function(h,...) {
   opt <- get("setupImport",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Select ImportData function"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="Selected script file:",container=tabval)
   tabval[1,2] <- gWidgets2::glabel(text=opt$importfile,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Select file:",container=tabval)
   tabval[2,2] <- gWidgets2::gbutton("Set file", container = tabval,handler = function(h, ...) { 
    ff <- mygfile("Select file",type="open")
    if(length(ff)==0) return()
    opt$importfile <- ff
    assign("setupImport",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileImport)    #save to file in installation folder
    gWidgets2::dispose(setwin) #remove subwindow
    #f_importsel(); #update gui window again after selecting new folder
   })
   gWidgets2::visible(setwin) <- TRUE
 }

  #The user can change algorithm options (nDone,maxContributors)
  f_options = function(h,...) { 
   radiotxt = c("ON","OFF")
   opt <- get("setupAdvanced",envir=nnTK) 
   setwin <- gWidgets2::gwindow(paste0("Advanced options"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   tabval[1,1] <- gWidgets2::glabel(text="Maximum contributors in QualLR (LRmix)",container=tabval)
   tabval[1,2] <- gWidgets2::gedit(text=opt$maxC1,width=w0,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Maximum contributors in QuanLR (EFM)",container=tabval)
   tabval[2,2] <- gWidgets2::gedit(text=opt$maxC2,width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text="Number of randoms in optimizer",container=tabval)
   tabval[3,2] <- gWidgets2::gedit(text=opt$nDone,width=w0,container=tabval)
   tabval[4,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 
    opt$maxC1 <- as.numeric(gWidgets2::svalue(tabval[1,2]))  #max number of contributors in LRmix model
    opt$maxC2 <- as.numeric(gWidgets2::svalue(tabval[2,2]))  #max number of contributors in EFM model
    opt$nDone <- as.numeric(gWidgets2::svalue(tabval[3,2]))  #iterations in optimizer
    assign("setupAdvanced",opt,envir=nnTK)  #assign user-value to opt-list
    setupWrite(unlist(opt),file=setupFileAdvanced)    #save to file in installation folder
    gWidgets2::dispose(setwin)
   } )
   gWidgets2::visible(setwin) <- TRUE
  }


 ################################################# 
 ##########FUNCTIONALITIES HELPFUNCTIONS########## 
 #################################################

 
 
 f_importData = function(h,...) { #wrapper function which calls other functions: importData and getStructuredData
   #REMOVE PREV. RESULTS WHEN NEW IMPORT:
   if(!is.null(get("DClist",envir=nnTK)) || !is.null(get("resCompMAC",envir=nnTK))) {
    gWidgets2::gmessage("Remember to press RESTART before importing a new case!")
   }

   caseID = gWidgets2::svalue( tabimportA[1,2] ) #get ID from table
   assign("caseID",caseID,envir=nnTK) #store ID in environment
   fn <- list.files(path=paste0(casedir,.sep,caseID), pattern="",full.names=TRUE) #get full names 

   print("-----------------------------")
   print("-------IMPORTING DATA--------")
   #NB: This is a tailored function which must be created such that it provides a list with Evidence and Reference table:
   #Evidences: colnames(data$mix)=("SampleName","Markers","Alleles","Heights"). Separated with "/" in text format.
   #References: colnames(data$ref)=("SampleName","Markers","Alleles"). Separated with "/" in text format.

   #importData
   found <- FALSE
   if(!is.na(get("setupImport",envir=nnTK)$importfile) && get("setupImport",envir=nnTK)$importfile!="") {
    source(get("setupImport",envir=nnTK)$importfile) #Source selected file
    if(exists('importData', mode='function')) found <- TRUE
   }
   if(!found) gWidgets2::gmessage("The user must specify a script file with R-function importData inside.\nGo to Settings->ImportData function")
   data <- list(mix= matrix(nrow=0,ncol=4),ref= matrix(nrow=0,ncol=3)) 
   markers <- numeric()
   metalist <- list() #contains list of table-elements
   consdata <- NULL #default value

   for(ff in fn) { #for each files: 
#ff=fn[3]
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
           metalist[[elem]] <- rbind(metalist[[elem]],tmplist[[elem]]) #add to matrix
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
   print("-------STRUCTURING DATA-------")
   datalist <- getStructuredData(data,ln=toupper(markers),minLoc=get("setupThresh",envir=nnTK)$minLociSS) #get Data in both List-format and Table-format (mixDataTABLE,refDataTABLE,mixDataLIST)

   ##Notice: The alleles in a loci should be ordered: Hence samples with same alleles can be detected   
   #STORE DATA (BOTH TYPES) -> EASY AVAILABLE THROUGH ENVIRONMENT
   assign("mixDataTABLE",datalist$mixDataTABLE,envir=nnTK) #assign to nnTK-environment
   assign("refDataTABLE",datalist$refDataTABLE,envir=nnTK) #assign to nnTK-environment
   assign("mixDataMATCHSTATUS",datalist$mixDataMATCHSTATUS,envir=nnTK) #assign to nnTK-environment
   assign("mixDataLIST",datalist$mixDataLIST,envir=nnTK) #assign to nnTK-environment
  # assign("refDataLIST",datalist$refDataLIST,envir=nnTK) #NOT USED ANYMORE!
   assign("metaDataLIST",metalist,envir=nnTK) #assign to nnTK-environment
   assign("consDataTABLE",consdata,envir=nnTK) #list of imported consensus data

   updateTables() #update datables
   refreshTabLIST() #update mixture-list
 } #end function

 f_runEuroformix = function(h,...) {  #Opens gui of euroformix directly
  last <- get("clicktableLast",envir=nnTK) #assign to nnTK-environment
  mixL <- get("mixDataLIST",envir=nnTK)[last$mix]
  refL <- tabToListRef(tab=get("refDataTABLE",envir=nnTK)[last$ref,,drop=FALSE],forceDi=TRUE)

  load(setupFileExport) #load emtpy object with ESX17 kit popfreq in euroformix 
  mmTK$mixData <- mixL
  mmTK$refData <- refL
  save(mmTK,file=setupFileExport2,compress="xz") #save object 
  euroformix::efm(envirfile=setupFileExport2) #run efm with saved file
 }

 f_saveprof = function(h,...) { #function to store profile
   if(h$action==0) names(Data) <- paste0("stain",sample(100,1))
   if(h$action>0) names(Data) <- paste0("ref",h$action)
   Data = sample_listToTable(Data) #convert from table to list
   saveTable(Data,"csv")
 }

 f_exportData = function(h,...) {  #allows user to export selected samples
  last <- get("clicktableLast",envir=nnTK) #assign to nnTK-environment
  mixL <- get("mixDataLIST",envir=nnTK)[last$mix]
  refL <- tabToListRef(tab=get("refDataTABLE",envir=nnTK)[last$ref,,drop=FALSE],forceDi=FALSE) #NOT FORCING DUP alleles
  evids <- names(mixL)
  refs <- names(refL)

  selwin <- gWidgets2::gwindow("Export data", visible=TRUE)
  tabtmp <- gWidgets2::glayout(container=selwin)

  tabSoft = gWidgets2::glayout(spacing=0,container=(tabtmp[1,1] <- gWidgets2::gframe("Open directly in software",container=tabtmp))) 
  tabSel = gWidgets2::glayout(spacing=0,container=(tabtmp[2,1] <- gWidgets2::gframe("Data selection",container=tabtmp)))  
  tabExp = gWidgets2::glayout(spacing=0,container=(tabtmp[3,1] <- gWidgets2::gframe("Data export",container=tabtmp)))  

  tabSoft[1,1] = gWidgets2::gbutton(text="EuroForMix",container=tabSoft,handler=f_runEuroformix)
  tabSel[1,1] = gWidgets2::glabel(text="Evidence(s)",container=tabSel)
  tabSel[1,2] = gWidgets2::glabel(text="Reference(s)",container=tabSel)
  tabSel[2,1] <- gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
  tabSel[2,2] <- gWidgets2::gcheckboxgroup(items=refs,container=tabSel,checked=TRUE)

  tabExp[2,1] = gWidgets2::gcheckbox(text="Export\npeak heights",container=tabExp,checked=TRUE) #export peak heights?
  tabExp[1,1] = gWidgets2::gbutton(text="Store\nevidence(s)",container=tabExp,handler=function(h,...) { 
   if(length(mixL)>0)  saveTable( sample_listToTable(mixL,PH=gWidgets2::svalue(tabExp[2,1])) )
  })
  tabExp[1,2] = gWidgets2::gbutton(text="Store\nreference(s)",container=tabExp,handler=function(h,...) { 
   if(length(refL)>0) saveTable(sample_listToTable(refL))
  })
 }

 f_calcRMP = function(h,...) {  #Function to calculate RMP for each references and RMNE for all evidence
  sig <- 3 #number of signif values
  refL <- get("refDataTABLE",envir=nnTK) 
  mixL <- get("mixDataLIST",envir=nnTK) #import mixtures
  ok = setPopFreq() #import population frequency from last selected file
  if(!ok) return() 
  popL <- get("popFreq",envir=nnTK) #import popFreq
  sn <- names(mixL) #get evid names
  rn <- rownames(refL) #get ref names
  locs <- names(popL) #loci to consider
  locs  <- setdiff(locs,"AMEL") #remove AMEL
  print( paste0("Considered loci from freq file: ",paste(locs,collapse="/")) )
  minF = 0.001 #minimium freq inserted when missing alleles
 
  calcRMNE = function(X) { #import a list X[[loc]] with allele vector and compare incl alleles in popFreqs
   #notice newly observed alleles are assumed to have allele freq=0
   rmne = 1 #default rmne value
   for(loc in locs) {
     if( is.null(X[[loc]]) || length(X[[loc]])==0) next; #skip if empty or not existing
     av <- names(popL[[loc]]) #get allele names
     freqs <- popL[[loc]][match( X[[loc]],av )] #[av%in%X[[loc]]] #get allele freqs
     freqs[is.na(freqs)] <- minF #.001 #insert minimum allele equal 0.001 if prev. not seen, alternative=min(unlist(popL))
     rmne <- rmne*sum(freqs)^2 #apply formula
   }
   return(rmne) #return calculated
  }

  
  #CALCULATE RMNE FOR EVIDENCE
  evidList <- matrix("",ncol=3,nrow=length(sn))
  colnames(evidList) <- c("ID","SampleName","RMNE")
  for(ss in sn) { #for each evidence
    ind <- which(sn==ss)
    rmne <- calcRMNE(X=lapply(mixL[[ss]][locs],function(x) x$adata)) #return only adata: DONT CONSIDER THE DETECTION THRESHOLD HERE
    evidList[ind,] <- c(paste0("#",ind),ss,format(rmne,digits=sig))
  }

  #CALCULATE RMP FOR REFERENCES
  popL = lapply(popL,function(x) x/sum(x)) #require sum to 1
  Glist = euroformix::getGlist(popL)
  rmp = rep(1,nrow(refL)) #default rmp value
  for(loc in locs) {
     sub = refL[,toupper(colnames(refL))==toupper(loc)] 
     G = paste0(Glist[[loc]]$G[,1],"/",Glist[[loc]]$G[,2])
     indUse = match(sub,G)
     isna = is.na(indUse) #get those not found
     rmp[!isna] = rmp[!isna]*Glist[[loc]]$Gprob[indUse[!isna]]
     rmp[isna] = rmp[isna]*minF #IN SITUATION OF MISSING 
  }
  refList = cbind(paste0("#",1:length(rn)),rn,format(rmp,digits=sig))
  colnames(refList) <- c("ID","SampleName","RMP")
  resRMP <- list(evid=evidList,ref=refList) #store results
  assign("resRMP",resRMP,envir=nnTK)  #store random match prob results (to be shown in report)

  setwin <- gWidgets2::gwindow(paste0("Random match probability calculations"),visible=FALSE) 
  tabval <- gWidgets2::ggroup(container=gWidgets2::gframe( paste0("Evidence(s)",paste0(rep("\t",4),collapse=""),"Reference(s)"),container=setwin,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
  guitab1 <- gWidgets2::gtable(items=evidList,container = tabval) 
  guitab2 <- gWidgets2::gtable(items=refList,container = tabval) 
  gWidgets2::add(tabval,guitab1,expand=TRUE,fill=TRUE) 
  gWidgets2::add(tabval,guitab2,expand=TRUE,fill=TRUE) 
  gWidgets2::visible(setwin) <- TRUE 
 } #end function


 f_calcIBS = function(h,...) {  #Function to calculate IBD between references
 # return()
  DBref <- get("refDataTABLE",envir=nnTK) #consider lists
  if(nrow(DBref)>=nLarge) {
    cat(paste0("The number of references were massive (>",nLarge,").\n Only references given in MatchStatus are considered!\n"))
    cand <- get("mixDataMATCHSTATUS",envir=nnTK) #get matchstatus lists
    cand = unique(cand[cand!="mixture"]) #consider only refs which fits to single sources 
    DBref = DBref[rownames(DBref)%in%cand,]
  }
  allrefsn <- rownames( DBref )
  nR <- length(allrefsn)
  if(nR<2) return() #return if less than 2
  locs <- colnames(DBref)
  ok = setPopFreq() #import population frequency from last selected file
  if(ok) locs = locs[ toupper(locs)%in%toupper(names(get("popFreq",envir=nnTK))) ] #Use ONLY THOSE IN POPFREQ (AVOIDS Y-markers) 
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
     gWidgets2::gmessage("Could not find R-package bigmemory. Please install it, or increase the memory and try again!")   
    )
  } 
  #describe(mac)
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
   tab <- cbind(mac[tabind],paste0(allrefsn[tabind[,1]]," - ",allrefsn[tabind[,2]]))
   ord <- order(as.numeric(tab[,1]),decreasing=TRUE)
   out <- tab[ord,,drop=FALSE]
   colnames(out) <- c("IBS","Comparison")
   assign("resIBS",out,envir=nnTK)  #store random match prob results (to be shown in report)

   #show candidates in GUI:
   setwin <- gWidgets2::gwindow(paste0("Similarity of references"),visible=FALSE) 
   tabval <- gWidgets2::ggroup(container=gWidgets2::gframe("Comparisons",container=setwin,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
   guitab1 <- gWidgets2::gtable(items=out,container = tabval) 
   gWidgets2::add(tabval,guitab1,expand=TRUE,fill=TRUE) 
   gWidgets2::visible(setwin) <- TRUE 
  } else {
   gWidgets2::gmessage("No candidates found")
   assign("resIBS",matrix(nrow=0,ncol=1),envir=nnTK)  #store empty IBS 
  }
 }

 f_calcIBSdist = function(h,...) {  #Function to calculate IBD between references
  ok = setPopFreq() #import population frequency from last selected file
  if(!ok) return() 
  dist <- getIBSdistr(get("popFreq",envir=nnTK)) #simulate random match probas
  tab <- 1-cumsum(dist[,2]) #consider cumulative probs
  names(tab) <- dist[,1]
  barplot(tab,main="Random allele sharing under unrelatedness",xlab="x=number of shared alleles",ylab="Prob. sharing>=x") 
  #abline(h=0.05,lty=2)
  print(tab)
 }

 #Function executed when double clicked one of the data tables
 clicktable = function(h,...) {
  last <- get("clicktableLast",envir=nnTK) #assign to nnTK-environment
  mixL <- get("mixDataLIST",envir=nnTK)
  isLUS = all(unlist(sapply(mixL,function(x)  sapply(x,function(y) all(grepl(LUSsymbol,y$adata),na.rm=TRUE)) )))  #ADDED: check if alleles are given on LUS format 
  kit0 = get("setupKit",envir=nnTK)$kitname
  if(is.null(kit0)) {
    gWidgets2::gmessage("Please specify kit in order to show EPG!")
    return()
  }
  if(h$action=="mix") {
   id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))

   for(id2 in id) { #print out evidence profiles 
       print(names(mixL)[id2]) 
       out = rbind(sapply(mixL[[id2]],function(x) paste0(x$adata,collapse="/")),sapply(mixL[[id2]],function(x) paste0(x$hdata,collapse="/")))
       rownames(out) = c("Alleles","Heights")
       print(out) 
   }    
   last <- list(mix=id) #store id, refs deleted
   assign("clicktableLast",last,envir=nnTK) #

   dev.new(width=25, height=10)
    if(isLUS) {
     if(length(id)>1) gWidgets2::gmessage("Not able to show multiple replicates for LUS variant data")
     euroformix::plotLUS(mixData=mixL[[id]],sn=names(mixL[id]),LUSsymbol=LUSsymbol) 
    } else {
     euroformix::plotEPG(mixL[id],kitname=kit0)
#Data=mixL[id];kitname=kit0;threshT=0;refcond=NULL;showPH=FALSE
    }
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    par(op)
    focus(mainwin)
   
  } 
  if(h$action=="ref" && !is.null(last$mix)) {
   mixL <- get("mixDataLIST",envir=nnTK)
   id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
   allid <- c(last$ref,id) #remember the other refs also
   last <- list(mix=last$mix,ref=allid) #store ids
   assign("clicktableLast",last,envir=nnTK) #
   refL <- tabToListRef(tab=get("refDataTABLE",envir=nnTK)[allid,,drop=FALSE],forceDi=FALSE) #NOT FORCING DUP alleles

   dev.new(width=25, height=10)
    if(isLUS) {
     if(length(id)>1) gWidgets2::gmessage("Not able to show multiple replicates for LUS variant data")
     euroformix::plotLUS(mixData=mixL[[id]],sn=names(mixL[id]),refData=refL,LUSsymbol=LUSsymbol) 
    } else {
     if( canPrintEPG() ) {
      euroformix::plotEPG(mixL[last$mix],refcond=refL,kitname=get("setupKit",envir=nnTK)$kitname)
     } else {
      gWidgets2::gmessage("Specify a valid kit to show EPG!")
     }
    }
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    par(op)
    focus(mainwin)   
  }
 }

 f_calcLRefm = function(h,...) { #function to run LR for all comparisons having LR>LRthresh(Qual)
   getMatchesLR(type="quan") #run EFM 
   refreshTabLIST2() #update QUAN LR table with results

   #Create matchlist (Final step)
   createMatchlist(modtype=2)
   refreshTabLIST() #update tables  
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
     lrval2 <- format( as.numeric( LRtab0[ LRtab0[,1]==evid & LRtab0[,2]%in%condRefs,4] ),digits=4) #get LR for ref
   }
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
  suppressWarnings({
   dc <- euroformix::deconvolve(contFit,maxlist=1) #get top candidate profiles
  })
  euroformix::plotTopEPG(MLEobj=contFit,DCobj=dc,kitname=get("setupKit",envir=nnTK)$kitname,threshT=get("setupModel",envir=nnTK)$threshT)
 }

 #Helpfunction to get refList for specific given refs
 getRefL = function(refs,forceDi=FALSE) { #return list with same order as for refs
  reftab=get("refDataTABLE",envir=nnTK)
  refL <- tabToListRef(tab=reftab[match(refs,rownames(reftab)),,drop=FALSE],forceDi=TRUE) #FORCING DUP alleles
  rm(reftab);gc()
  return(refL)
 }
 
 #FUNCTION WHICH PERFORMS DC (uses settings in GUI)
 doDC = function(nC,evids=NULL,refs=NULL,showPlot=TRUE) {
   refData = NULL
   if(!is.null(refs) && length(refs)>0) {
     refData = getRefL(refs,forceDi=FALSE) # get list data
   }
   evidData <- get("mixDataLIST",envir=nnTK)[evids] #evidence to consider
   evid = paste0(evids,collapse="/") #collapse multiple evidence names
   condrefs = paste0(refs,collapse="/")
    suppressWarnings({ #calculate both MLE fit and DC:
	 kit0 = NULL #no kit specified by default
	 xi0 = 0 #no stutters assumed by default
      if( get("setupModel",envir=nnTK)$degrad=="ON" ) kit0 = get("setupKit",envir=nnTK)$kitname #assign kit
      if( get("setupModel",envir=nnTK)$stutt=="ON" ) xi0 = NULL
      mlefit <- calcHp(evidData,refData,nC=nC,popFreq=get("popFreq",envir=nnTK),kit=kit0,pC=get("setupModel",envir=nnTK)$dropinC,lambda=get("setupModel",envir=nnTK)$dropinL,threshT=get("setupModel",envir=nnTK)$threshT, xi=xi0 , nDone=get("setupAdvanced",envir=nnTK)$nDone) #get fitted object
      dc <- euroformix::deconvolve(mlefit,maxlist=1) #get top candidate profiles
      if(showPlot) euroformix::plotTopEPG(MLEobj=mlefit,DCobj=dc,kitname=get("setupKit",envir=nnTK)$kitname,threshT=get("setupModel",envir=nnTK)$threshT)
    })

    #INSERTING CANDIDATE DECONVOLED PROFILES:
    ratio <- get("setupThresh",envir=nnTK)$ratio #get ratio-threshold
    probA <- get("setupThresh",envir=nnTK)$probA #get probability of allele - threshold

    locs <- names(dc$toprankGi) #get loci
    contrs = colnames(dc$toprankGi[[1]]) #get contributors
    candtab <- matrix(nrow=0,ncol=length(locs)+4) 
    colnames(candtab) <- c("Component","Conditional","nC","MixProp",locs)
    nR = length(refData) #number of conditional refs
    for(cind in 1:length(contrs)) { #for each contributors
      compn <-  contrs[cind] #paste0(evid,"_C",uind) #component name
      addRef = FALSE #Should the profile be added? (Must have at least 1 deduced allele for an unknown component)
      mxhat <- mlefit$fit$thetahat2[cind] #get mixture proportion
      newrow <- rep(NA,length(locs))
      for(loc in locs) { #for each locus 
        insind <- which(locs==loc) #insert index
        cand <- dc$toprankGi[[loc]][,cind] #get candidate  
        if(!is.na(cand[3]) &&  as.numeric(cand[3])<ratio) { #if not a likely genotype
         ind <- which(dc$table4[,1]==compn & dc$table4[,2]==loc)[1] #find top ranked single allele
         if( as.numeric(dc$table4[ind,4])>= probA )  newrow[insind] <- dc$table4[ind,3] #insert allele if prob>probA
        } else {
         newrow[insind] <- cand[1]
        }
        if( cind <= nR && length(refData[[cind]][[loc]]$adata)==0 && !is.na(newrow[insind]) ) addRef  = TRUE #indicate that ref prof. should be added
      } #end for each loci
      if(all(is.na(newrow))) next #skip to next
      if( cind <= nR && !addRef) next #skip if ref should not be added
      if(cind <= nR ) compn = names(refData)[cind] #use ref name instead if conditioned on
      newrow <- c(paste0(evid,"-",compn),condrefs,nC,signif(mxhat,2),newrow)
      candtab <- rbind(candtab, newrow)
    } 
    return(candtab)
 } #end doDC

 getAlleles = function(x) {
   if(length(x)==1) return(as.character()) #dont consider if only 1 allele
   return(x)
 }

 #Helpfunction to calc LR (fit quan model for both hyps)
 fitEFMHYPs = function(nC,evids,ref,condref=NULL) {
    ok = setPopFreq()
    if(!ok) return()
  
    #Prepare params:
    pC=get("setupModel",envir=nnTK)$dropinC;
    lambda=get("setupModel",envir=nnTK)$dropinL;
    threshT=get("setupModel",envir=nnTK)$threshT;
    kit = NULL #no kit specified by default
    xi = 0 #no stutters assumed by default
    if( get("setupModel",envir=nnTK)$degrad=="ON" ) kit = get("setupKit",envir=nnTK)$kitname #assign kit
    if( get("setupModel",envir=nnTK)$stutt=="ON" ) xi = NULL
    if( !is.null(kit) && !canPrintEPG() ) {
     gWidgets2::gmessage("Please specify a kit for further comparisons.\nGo to Settings->Kit selection")
     return() 
    } 

    #Prepare data:
    popFreq=get("popFreq",envir=nnTK);
    samples <- get("mixDataLIST",envir=nnTK)[evids] #consider lists
    uselocs <-  names(popFreq)[toupper(names(popFreq))%in%toupper(names(samples[[1]]))] #loci to consider
    if(is.null(xi)) uselocs <- setdiff(uselocs,"AMEL") #EXCLUDE AMEL IF STUTTER CONSIDERD (QUICK SOLUTION) 
    samples = lapply(samples,function(x) x[uselocs])
    allrefs = c(ref,condref) #get all refs. POI is first ref
    refD = getRefL(allrefs,forceDi=FALSE) # get list data
    refD2 <- list() #reference list must have other structure than considered here...
    for(loc in uselocs) refD2[[loc]] <- lapply(refD ,function(x) getAlleles(x[[loc]]$adata)) #extract reference
    data <- euroformix::Qassignate(samples, popFreq[uselocs],refD2,incS=FALSE,incR=FALSE) #don't include stutters, use all loci

    condhp <- 1:length(allrefs) #conditional order must be increasing order
    condhd <- condhp - 1 #don't condition on POI under Hd

    #RUN CALCULATIONS:
    nDone=get("setupAdvanced",envir=nnTK)$nDone
    fitMLE = function(condOrder) {
     nU = nC-sum(condOrder>0)
     if(nU>2) fit <- euroformix::contLikMLEpara(nC,data$samples,data$popFreq,data$refData,condOrder,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=TRUE)
     if(nU<=2) fit <- euroformix::contLikMLE(nC,data$samples,data$popFreq,data$refData,condOrder,xi=xi,prC=pC,lambda=lambda,nDone=nDone,threshT=threshT,kit=kit,verbose=TRUE)
     return(fit)
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
  selwin <- gWidgets2::gwindow(paste0("Quan LR calculation (EFM)"), visible=TRUE)
  tabSel <- gWidgets2::glayout(container=selwin)
  tabSel[1,1] = gWidgets2::glabel(text="Evidence(s):",container=tabSel)
  tabSel[1,2] = gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
  tabSel[2,1] = gWidgets2::glabel(text="POI:",container=tabSel)
  tabSel[2,2] = gWidgets2::glabel(text=ref,container=tabSel)
  tabSel[3,1] = gWidgets2::glabel(text="Number of contributors=",container=tabSel)
  tabSel[3,2] = gWidgets2::gcombobox(items=Krange,selected=selind,editable=TRUE,container=tabSel) 
# gWidgets2::gedit(text=nC,container=tabSel)
  tabSel[4,1] = gWidgets2::glabel(text="Condition on:",container=tabSel)
  tabSel[4,2] = gWidgets2::gcheckboxgroup(items=condRefs,container=tabSel,checked=FALSE) #don't select by default
  if(!is.null(lrvals)) { #if lr values given
   tabSel[1,3] = gWidgets2::glabel(text="log10 LR",container=tabSel)
   tabSel[2,3] = gWidgets2::glabel(lrvals,container=tabSel)
  }
  if(!is.null(lrvals2)) { #if lr values given
   tabSel[4,3] = gWidgets2::glabel(lrvals2,container=tabSel)
  }
  tabSel[5,1] = gWidgets2::gbutton(text="Calculate",container=tabSel,handler=function(h,...) {
    nC2 <- as.integer(svalue(tabSel[3,2])) #get selected number of contributors
    if(nC2>get("setupAdvanced",envir=nnTK)$maxC2) gWidgets2::gmessage("Notice: You have selected more contributors than given in Advanced options.",title="Warning")

    condsel <- NULL #conditional references
    if(!is.null(condRefs) && length(condRefs)>0) {
     condsel <- svalue(tabSel[4,2]) #get selected ref to condition on
    }
    mixsel <-  svalue(tabSel[1,2]) #get selected samples
    if(length(mixsel)==0) {
      gWidgets2::gmessage("Please select at least one evidence")
      return()
    }
    gWidgets2::dispose(selwin) #remove window 
    ret = fitEFMHYPs(nC=nC2,evids=mixsel,ref=ref,condref=condsel) 

    #STORE MATCH HERE (not in function directly)  
    matchlist0 = get("resCompLR",envir=nnTK)  #get stored results from Qual/Quan comparison (unsorted, but truncated)
    matchlist1 = get("resCompLR1",envir=nnTK)  #get stored qual based results from comparison (sorted wrt QualLR)
    matchlist2 = get("resCompLR2",envir=nnTK)  #get stored quan based results from comparison (sorted wrt QuanLR)
    hpfitlist = get("storedFitHp",envir=nnTK)  #get stored model results under Hp

    condREF = FALSE #bool if cond. ref is included
    hasCondREF = FALSE #whether condtional ref columns is in table
    if(!is.null(condsel) && length(condsel)>0) condREF = TRUE #add cond.refs to table?
    if(!is.null(matchlist2) && "condRef"%in%colnames(matchlist2)) hasCondREF = TRUE 
    if(!is.null(matchlist2) && condREF && !hasCondREF) { #need to extend the table
      matchlist2 = cbind(matchlist2,"") #add an extra column
      colnames(matchlist2)[ncol(matchlist2)] = "condRef" #add name
    } 
     
    #NEED TO REMOVE ALREADY EXISTING COMPARISON AND EXCHANGE WITH NEW RESULTS:
    ind = which(matchlist0[,1]==mixsel & matchlist0[,2]==ref)  #get index in matchlist0 (hp stored in this index)
    if(length(ind)==0) { #NB: this variant was not computet before
      ind = which(matchlist1[,1]==mixsel & matchlist1[,2]==ref)  #get index in matchlist0 (hp NOT stored in this index)
      newrow = matchlist1[ind,]
    } else { #Found
      newrow = matchlist0[ind,]
    }
    if(condREF || hasCondREF) newrow = c(newrow,"") #add column if condRef considered OR before
    newrow[4] = signif(ret$mleLR,digits=4) #round
    newrow[5] = nC2 #insert number of contributors
    if(condREF) newrow[6] = paste0(condsel,collapse="/") #show conditionalRefs in table
    if(is.null(matchlist2) && condREF) names(newrow)[6] = "condRef" #need to add name
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

    #STORE HPFIT RESULTS
    hpfitlist = get("storedFitHp",envir=nnTK)  #get already stored objects
    if(is.null(hpfitlist)) hpfitlist = replicate(nrow(matchlist1),list()) #init. list if first time
    hpfitlist[[ind]] <- ret$fithp  #insert on right index
    assign("storedFitHp",hpfitlist,envir=nnTK)  #store object

    #UPDATE TABLES:
    refreshTabLIST2(which(ord==insInd)) #update tables with results, with marked on selected one
    createMatchlist(modtype=3) #final step is to update MIXTURES table
    refreshTabLIST(evidsel=mixsel) #update match-tables: mark on considered evidence  
    gWidgets2::svalue(nb) <- 4 #go to quan LR result tab when done
  }) #end button
 }

 #FUNCTION TO SPECIFY HYPOTHESIS AND PERFORMS A SINGLE DC
 createHypDCWindow = function(evids,refs,nC,lrvals=NULL) { 
  Krange = 1:4 #range of number of contributors
  selind = which(nC==Krange)
  if(length(selind)==0) selind = max(Krange)
  selwin <- gWidgets2::gwindow(paste0("Deconvolution/Show expected peak heights"), visible=TRUE)
  tabSel <- gWidgets2::glayout(container=selwin)
  tabSel[1,1] = gWidgets2::glabel(text="Evidence(s):",container=tabSel)
  tabSel[1,2] = gWidgets2::gcheckboxgroup(items=evids,container=tabSel,checked=TRUE)
  tabSel[2,1] = gWidgets2::glabel(text="Number of contributors=",container=tabSel)
  tabSel[2,2] = gWidgets2::gcombobox(items=Krange,selected=selind,editable=TRUE,container=tabSel) 
# gWidgets2::gedit(text=nC,container=tabSel)
  tabSel[3,1] = gWidgets2::glabel(text="Condition on:",container=tabSel)
  tabSel[3,2] = gWidgets2::gcheckboxgroup(items=refs,container=tabSel,checked=TRUE)
  if(!is.null(lrvals)) { #if lr values given
   tabSel[2,3] = gWidgets2::glabel(text="log10 LR",container=tabSel)
   tabSel[3,3] = gWidgets2::glabel(lrvals,container=tabSel)
  }
  tabSel[4,1] = gWidgets2::gbutton(text="Calculate",container=tabSel,handler=function(h,...) {
    nC2 <- as.integer(svalue(tabSel[2,2])) #get selected references
    refsel <- NULL #conditional references
    if(length(refs)>0) {
     refsel <- svalue(tabSel[3,2]) #get selected refs
    }
    mixsel <-  svalue(tabSel[1,2]) #get selected samples
    if(length(mixsel)==0) {
      gWidgets2::gmessage("Please select at least one evidence")
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
 clickDClist = function(h,...) {
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   id <- as.integer(gsub("#","",gWidgets2::svalue(h$obj)))
   answ <- gWidgets2::gconfirm(paste0("Do you want to add the deconvoluted profile\n",DClist[id,1],"\nto the references?"))
   if(answ) f_addref(h=list(action=DClist[id,])) #open edit window of references
   #assign("DClistReport",NULL,envir=nnTK)  #store confirmed deconvoluted reference candidates to own list (used for reporting)
 }

 #function which takes all matches (with LR>threshold) and create a list to double click on (showing confirming under all conded)
 createMatchlist = function(modtype) { #directly after calculations are done
  #modtype: 1=All Qual LR, 2=All Quan LR, 3=Original Qual LR, but some updated with Quan LR
  threshLR <- get("setupThresh",envir=nnTK)$LRthresh1 #QualLR used by default
  tab <- get("resCompLR1",envir=nnTK) #QualLR used by default
  if(modtype==2)  {
    threshLR <- get("setupThresh",envir=nnTK)$LRthresh2 #QuanLR used if type 2
    tab <- get("resCompLR2",envir=nnTK)
  }
  score <- as.numeric(tab[,4])
  tab <- tab[score>=log10(threshLR),,drop=FALSE] #combinations to consider (all above LR treshold)
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
   resMatches <- numeric() #create table
   for(evid in unEvid) { #for each evidence we condition on all evidence 
    evidind <- tab[,1]==evid
    refs <- tab[evidind,2] #get references
    nC <- as.numeric(tab[evidind,5][1]) #get number of contributors (equal for all)
    resMatches <- rbind(resMatches, c(evid,paste0(refs,collapse="/"),nC) )
   }
  } else {
   resMatches =  matrix(nrow=0,ncol=3)
  }
  colnames(resMatches) <- c("Evidence","Reference(s)","Num Contr.")
  assign("resMatches",resMatches,envir=nnTK) 
 }

 #Show matches in a graph:
 showMatchNetwork = function(h,...) {
   require(igraph)
   tab <- get("resCompLR",envir=nnTK) #get last stored LR comparison values
   if( is.null(tab)) return(FALSE)

   threshLR <- get("setupThresh",envir=nnTK)$LRthresh2 #can be dynamically changed..
   if( is.null(get("resCompLR2",envir=nnTK))) {  #IF ONLY QUAL BASED LR IS CONSIDERED
      threshLR <- get("setupThresh",envir=nnTK)$LRthresh1 #can be dynamically changed..
   }
   score <- as.numeric(tab[,4])
   rem <- duplicated(tab[,1:2]) #indices to remove
   if(nrow(tab)==1) rem <- FALSE #duplicated doesnt work for case having only 1 row
   rem <- rem | score<log10(threshLR) #remove matches below LR threshold
   tab2 <-  data.frame(from=tab[!rem,2],to=tab[!rem,1],weight=sqrt(score[!rem]),nC=as.integer(tab[!rem,5]))
   gg <- igraph::graph.data.frame(tab2,directed=FALSE)
   nods <- igraph::V(gg)

   #Colorcode circles wrt number of contributors
   cols <- rep("green",length( nods  ))
   cols[grepl("Unknown ",names(nods))] <- "cyan" #unknowns 
   cols[names(nods)%in%tab2$to[tab2$nC==2]] <- "orange" 
   cols[names(nods)%in%tab2$to[tab2$nC>2]] <- "red" 
   plot(gg,edge.width=E(gg)$weight,vertex.color=cols,vertex.size=10,vertex.label.color="black",vertex.label.cex=0.8,main=paste0("Matches for ",get("caseID",envir=nnTK)) )
   return(TRUE)
 }


####################################
#STEP 1) CALCULATE MATCHING ALLELES#
####################################
 #Function executed when clicking "Comparison"
 getMatchesMAC = function(locs=NULL) { #Compare only alleles (given selection of loci
  DBmix <- get("mixDataTABLE",envir=nnTK) #consider lists
  DBref <- get("refDataTABLE",envir=nnTK) #consider lists
  DBmixmatch <- get("mixDataMATCHSTATUS",envir=nnTK) #consider lists

  #Add evidence profiles as "mixture" considered
  DBmix <- DBmix[DBmixmatch=="mixture",,drop=FALSE] #All evidences must be mixtures

  if(!is.null(locs)) {
   keep = toupper(colnames(DBmix))%in%toupper(locs)
   DBmix <- DBmix[,keep,drop=FALSE] #use relevant loci
   DBref <- DBref[,keep,drop=FALSE] #use relevant loci
  }

  #Perform the Matching allele counting algorithm-> Output is score for all comparisons
  matchMAC <- calcMACcomparison(DBmix=DBmix,DBref=DBref,threshMAC=get("setupThresh",envir=nnTK)$MACthresh) #calculating the score=normalized number of allele match counting for all combinations
  #print(matchMAC$MatchMatrix)
  assign("resCompMAC",matchMAC,envir=nnTK)  #store match-matrix in environment 
 } #end MAC comparison

############################################
#STEP 2) CALCULATE LR FOR REMAINING MATCHES#
############################################
 getMatchesLR = function(type="quan") { #Calculate LR for individuals in MatchList 
  #type={"qual","quan"} 
  Clist <- get("resCompMAC",envir=nnTK)$MatchList  #list to consider for calculating LR

  if(type=="quan" && !is.null(get("resCompLR1",envir=nnTK)) ) {
   Clist <- get("resCompLR1",envir=nnTK)  #list to consider for calculating LR (based on qualLR)
   Clist <- Clist[as.numeric(Clist[,4])>=log10(get("setupThresh",envir=nnTK)$LRthresh1),1:3,drop=FALSE] #keep only variants above thrshold AND COLUMNS in MAC
  }
  DBmix <- get("mixDataLIST",envir=nnTK)[unique(Clist[,1])] #get relevant evidence
  DBref <- get("refDataTABLE",envir=nnTK)
  DBref = DBref[rownames(DBref)%in%unique(Clist[,2]),,drop=FALSE] #get only relevant references

  #matchlist=Clist;popFreq=get("popFreq",envir=nnTK);kit=get("setupKit",envir=nnTK)$kitname;pC=get("setupModel",envir=nnTK)$dropinC;lambda=get("setupModel",envir=nnTK)$dropinL;threshT=get("setupModel",envir=nnTK)$threshT;maxC=6;p0=0.1;#xi=xi0
  suppressWarnings({
   if(type=="qual")  matchLRres <- calcQualLRcomparison(DBmix,DBref,matchlist=Clist,popFreq=get("popFreq",envir=nnTK),pC=get("setupModel",envir=nnTK)$dropinC,pD0=get("setupAdvanced",envir=nnTK)$pD0,maxC=get("setupAdvanced",envir=nnTK)$maxC1)
   if(type=="quan") {
      
      kit0 = NULL #no kit specified by default
	 xi0 = 0 #no stutters assumed by default
      if( get("setupModel",envir=nnTK)$degrad=="ON" ) kit0 = get("setupKit",envir=nnTK)$kitname #assign kit
      if( get("setupModel",envir=nnTK)$stutt=="ON" ) xi0 = NULL
      matchLRres <- calcQuanLRcomparison(DBmix,DBref,matchlist=Clist,popFreq=get("popFreq",envir=nnTK),kit=kit0,xi=xi0,pC=get("setupModel",envir=nnTK)$dropinC,lambda=get("setupModel",envir=nnTK)$dropinL,threshT=get("setupModel",envir=nnTK)$threshT,nDone=get("setupAdvanced",envir=nnTK)$nDone,maxC=get("setupAdvanced",envir=nnTK)$maxC2) 
   }
  })
  matchlist <- matchLRres$MatchList

  if(type=="qual") {
   assign("resCompLR",matchlist,envir=nnTK)  #store results from comparison

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
   assign("resCompLR",matchlist,envir=nnTK)  #store results from comparison

   #sort the matchlist with respect to the LRs:  
   LRcol <- which(colnames(matchlist)=="log10LR") #get column where LR is
   LRval <- as.numeric(matchlist[,LRcol])
   ord <- order( LRval,decreasing=TRUE)
   matchlist[,LRcol] <- round(LRval,2) #round to 2 dec
   #matchlist[ord,] #sort list by LR
   assign("resCompLR2", matchlist[ord,,drop=FALSE],envir=nnTK)  #store sorted matchlist 
  }
 } #end LR comparison


 #Function giving window for editing alleles for new references
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
    ret = gWidgets2::ginput("Define index of reference profiles to show (use comma, AND/OR minus-)", text = "1-10", title = "User input", icon = c("info"))
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
  defN <- "NEW" #default name
  refN <- rownames(refT) #refnames
  newTab <- cbind(refN,refT)
  colnames(newTab)[1] <- "SampleName"

  newline <- c(defN,rep(def,length(locs)))
  if(!is.null(h$action)) { #opened with deconvolution
    tmp <- h$action #get 
    newline[1] <- tmp[1] #get sample name (component)
    dat <- tmp[5:length(tmp)] #get suggested deconvolved elem
    newline[-1] <- dat[match(locs,names(dat))] #
    newline[is.na(newline)] <- ""
  }
  newTab <- rbind(newTab,newline) #add empty/new row

  setwin <- gWidgets2::gwindow(paste0("Edit or add new reference profiles"),visible=FALSE) 
  tabval = gWidgets2::ggroup(spacing=0,container=(setwin),horizontal=FALSE)  
  guitab <- gWidgets2::gdf(items=newTab,container = tabval) 
  gWidgets2::add(tabval,guitab,expand=TRUE,fill=TRUE) 
  savbut = gWidgets2::gbutton(text="Save profile",spacing=0,container=(tabval),handler = function(h, ...) { 
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
     answ <- gWidgets2::gconfirm(paste0("Are you sure you want to delete:\n",paste0(refN[delindsOld],collapse="/"),"?"))
     if(answ) { #if agree then refs are deleted
      refT <- refT[!delindsOld,,drop=FALSE] #update ref-table
      newref = newref[!delindsOld,] #update ref-table
      refN = refN[!delindsOld]
     }
    } 

    #CHECK FOR CHANGES (not added ref)
    if(nrow(refT)>0) { #must have at least one ref
     isChanged = which(refT != newref[,-1],arr.ind=TRUE)
     if(length(isChanged)>0) {
      text = paste0(refN[isChanged[,1]]," (",colnames(refT)[isChanged[,2]],")")
      answ <- gWidgets2::gconfirm(paste0("Are you sure you want to apply the changes in:\n",paste0(text,collapse="/"),"?"))
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
      gWidgets2::gmessage("Reference name already taken. Please consider a new name!")
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
    updateTables() #updates tables again 
  }) #end button
  gWidgets2::add(tabval,savbut,expand=TRUE,fill=TRUE) 
  gWidgets2::visible(setwin) <- TRUE 
 }  #end add ref

 #FUNCTION TO STORE TABLE:
 f_exporttable = function(h,...) { 
   tab <- NULL
   if(h$action=="mac")    tab <- addRownameTable(get("resCompMAC",envir=nnTK)$MatchMatrix,type=1)
   if(h$action=="qual")    tab <- get("resCompLR1",envir=nnTK)
   if(h$action=="quan")    tab <- get("resCompLR2",envir=nnTK)
   if(h$action=="final")    tab <- get("allMixList",envir=nnTK)
   if(h$action=="dc")    tab <- get("DClist",envir=nnTK)
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
    candtab <- doDC(nC,evid,refs,showPlot=FALSE)
    candtabs = rbind(candtabs,candtab) #add candidate table
    print(paste0(round(ss/nS* 100), "% DC calculation complete...")) 
  }
  if(length(candtabs)>0) { #if more than one new 
   dclist <- get("DClist",envir=nnTK) #get stored DC-list
   dclist <- rbind(dclist,candtabs) #add candidates
   assign("DClist",dclist,envir=nnTK) #get stored DC-list  
   refreshDCLIST() #refresh DC-list
   gWidgets2::svalue(nb) <- 6 #go to DC-tab after calculation
  }
} #end f_doDCall


####################
#######REPORT#######	
####################

  f_createreport = function(h,...) { #this function loads data and put them in a HTML script
    optReport <- get("setupReport",envir=nnTK)  #get report options
    reportname <- gWidgets2::ginput("Name report file", text="report", title="User input", icon="question")
    if(reportname=="") return() #report not created

    library(R2HTML) 
    print("Creating report...")
    graphics.off() #close all plots before running..
    brdt1 <- 1 #inner bord type
    brdt2 <- NULL #outer bord type
    caseID  <- get("caseID",envir=nnTK) #get ID fromenvironment
    if(is.null(caseID)) return()
    path <- paste0(casedir,.sep,caseID)
    path2 <- paste0(path,.sep,"report")
    dir.create(path2, showWarnings = FALSE) #create folder if not existing
    
    #GENERATE HTML CODE:
    #cssfile <- "http://www.stat.ucl.ac.be/R2HTML/Pastel.css"
    cssfile <- system.file("samples", "R2HTML.css", package="R2HTML")
    htmlf <- HTMLInitFile(path2,filename=reportname,CSSFile=cssfile,Title=paste0("Case ",get("caseID",envir=nnTK)))

    #Header:    
    HTML(as.title(paste0("Report for Case ",get("caseID",envir=nnTK))),file = htmlf,HR=1)

  if(optReport$checked[1]=="TRUE") { #show header?
    HTML( paste0("CaseSolver version ",version," (euroformix_",packageVersion("euroformix"),").") ,file = htmlf)#,HR=0)
    HTML( R.version.string ,file = htmlf)#,HR=0)
    HTML( paste0("User: ",Sys.getenv("USERNAME")) ,file = htmlf)#,HR=0)
    HTML( paste0("Created: ",Sys.time()) ,file = htmlf)#,HR=0)
  }

    #Provide data tables:
    mixLIST <- get("mixDataLIST",envir=nnTK) #get mix list
    mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #get data table from nnTK-environment
    refDataTABLE <- get("refDataTABLE",envir=nnTK) #get data table from nnTK-environment
    mixDataMATCHSTATUS <- get("mixDataMATCHSTATUS",envir=nnTK) #assign to nnTK-environment
    metaDataList <- get("metaDataLIST",envir=nnTK) #list of imported metadata
    consDataTABLE <- get("consDataTABLE",envir=nnTK) #list of imported metadata
    locs <- colnames(mixDataTABLE) #assume at least 1 evidence sample

    mixTab <- ssTab <- refTab <- NULL #empty
    if(!is.null(refDataTABLE)) {    #Add mix-table
     refTab <- addRownameTable(refDataTABLE,type=4)
    }
    if(!is.null(mixDataTABLE)) {    #Add evid-tables
     ssDataTABLE <-  cbind(mixDataMATCHSTATUS,mixDataTABLE)[mixDataMATCHSTATUS!="mixture",,drop=FALSE]
     mixDataTABLE <-  mixDataTABLE[mixDataMATCHSTATUS=="mixture",,drop=FALSE]
     colnames(ssDataTABLE)[1] <- "MatchStatus"
     if(!optReport$checked[20]) ssDataTABLE = ssDataTABLE[,-1,drop=FALSE] #drop MatchStatus column
     ssTab <-  addRownameTable(ssDataTABLE,type=4)
     mixTab <- addRownameTable(mixDataTABLE,type=4)
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

   # HTML(as.title(paste0("Data")),file = htmlf,HR=1)
    reftab <-  addRownameTable(get("refDataTABLE",envir=nnTK),type=4)

###############
###SHOW DATA###
###############
  if(optReport$checked[2]=="TRUE") { #show references?
    HTML(as.title("References"),file = htmlf,HR=2)
    if(!is.null(reftab)) {
      if(nrow(reftab)>=nLarge) print("THIS MAY TAKE A WHILE!")
      HTML(as.data.frame(reftab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
      HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }

  if(optReport$checked[3]=="TRUE") { #show SS profiles?
    HTML(as.title("Single source profiles"),file = htmlf,HR=2)
    if(!is.null(ssTab)) {
     HTML(as.data.frame(ssTab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }

  if(optReport$checked[4]=="TRUE") { #show Mix profiles?
    HTML(as.title("Mixture profiles"),file = htmlf,HR=2)
    if(!is.null(mixTab)) {
     HTML(as.data.frame(mixTab), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }

  if(optReport$checked[5]=="TRUE") { #Show consensus profiles?
    HTML(as.title("Consensus profiles"),file = htmlf,HR=2)
    if(!is.null(consDataTABLE)) {
     sn <- unique(consDataTABLE[,1]) #get sample names
     consDataOUT <- matrix(ncol=length(locs)+1,nrow=length(sn))
     consDataOUT[,1] <- sn
     for(ss in sn) { #for each samples
       for(loc in locs) { #for each locus
         consDataOUT[which(sn==ss),which(locs==loc)+1] <- consDataTABLE[consDataTABLE[,1]==ss &consDataTABLE[,2]==loc,3]
       }
     }
     colnames(consDataOUT) <- c("SampleName",locs)
     HTML(as.data.frame(consDataOUT), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }

###############
###SHOW W/PH###
###############
  if(optReport$checked[6]=="TRUE") { #should SS with peak heights be included?
    HTML(as.title("Single sources w/peak heights"),file = htmlf,HR=2)
    if(!is.null(mixLIST) && !is.null(ssTab)) {
     for(ss in ssTab[,1] ) { #for each single sources
       HTML(as.title( paste0(ss) ),file = htmlf,HR=3)
       if(nrow(allTabLIST[[ss]])==0) { #check if empty
         HTML(paste0("Empty."),file = htmlf)
       } else {
         HTML(as.data.frame(addRownameTable(allTabLIST[[ss]],type=0)), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
       }  
     }
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }

  if(optReport$checked[7]=="TRUE") { #should Mixture with peak heights be included?
    HTML(as.title("Mixtures w/peak heights"),file = htmlf,HR=2)
    if(!is.null(mixLIST) && !is.null(mixTab)) {
     for(ss in mixTab[,1] ) { #for each mixtures
       HTML(as.title( paste0(ss) ),file = htmlf,HR=3)
       if(nrow(allTabLIST[[ss]])==0) { #check if empty
         HTML(paste0("Empty."),file = htmlf)
       } else {
         HTML(as.data.frame(addRownameTable(allTabLIST[[ss]],type=0)), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
       }  
     }
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
  }


 if(optReport$checked[8]=="TRUE") { #Show metadata?
    HTML(as.title("Metadata"),file = htmlf,HR=2)
    if(length(metaDataList)>0) { #is empty list?
     for(elem in names(metaDataList)) {
       if(length(metaDataList[[elem]])==0) next #skip if no info
       HTML(as.title( paste0(elem) ),file = htmlf,HR=3)
       HTML(metaDataList[[elem]], file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)       
     }
    } else {
     HTML(paste0("None."),file = htmlf)#,HR=10)
    }
 }

    #Helpfunction
    easyAdd = function(X,type=2,inclRN=TRUE) {
     if(!is.null(X)) {
      HTML(as.data.frame(addRownameTable(X,type=type)), file= htmlf,row.names=inclRN,align="left",innerBorder = brdt1,Border=brdt2)
     } else {
      HTML(paste0("Not completed."),file = htmlf)#,HR=10)
     }
    }

#################
###COMPARISONS###
#################
    if(any(optReport$checked[9:12])) { #Any comparisons to show?
     HTML(as.title("Comparisons"),file = htmlf,HR=1)
    }
    #Provide comparison matrix:
    tab <- get("resCompMAC",envir=nnTK)$MatchMatrix
  if(optReport$checked[9]=="TRUE" && !is.null(tab) ) { #Show match list1?
    HTML(as.title("Comparison matrix"),file = htmlf,HR=2)
    if(nrow(tab)>0) {
     ind <- as.numeric(tab)<get("setupThresh",envir=nnTK)$MAC #get smaller indices
     tab[ind] <- "" #show only greater than threshold
     easyAdd(tab,type=4)
    }
  }

    #Provide match list 1 (qual based):
    tab <- get("resCompLR1",envir=nnTK)
  if(optReport$checked[10]=="TRUE" && !is.null(tab) ) { #Show match list1?
    HTML(as.title("Match list (Qual LR)"),file = htmlf,HR=2)
    if(nrow(tab)>0) {
     tab <- tab[as.numeric(tab[,4])>=log10(get("setupThresh",envir=nnTK)$LRthresh1),,drop=FALSE]
     easyAdd(tab,type=0)
    }
  }
    #Provide match list 2 (quan based):
    tab <- get("resCompLR2",envir=nnTK)
 if(optReport$checked[11]=="TRUE" && !is.null(tab) ) { #Show match list1?
    HTML(as.title("Match list (Quan LR)"),file = htmlf,HR=2)
    if(nrow(tab)>0) {
     tab <- tab[as.numeric(tab[,4])>=log10(get("setupThresh",envir=nnTK)$LRthresh2),,drop=FALSE]
     easyAdd(tab,type=0)
    }
 }
    #Provide match list 2:
 if(optReport$checked[12]=="TRUE") { #Show match list2?
    HTML(as.title("Final match list (w/all mixtures)"),file = htmlf,HR=2)
    easyAdd(get("allMixList",envir=nnTK),type=0) #show all mixtures
 }

 if(optReport$checked[13]=="TRUE") { #Show match network?
     HTML(as.title("Match network"),file = htmlf,HR=2)
    okplot <- !is.null(get("resCompLR",envir=nnTK))
    if(okplot) {  #only add if OK
     netf <- file.path(path2,"matchnetwork.png") #file of picture 
     png(netf,width = 2000, height = 2000,res=200)
     showMatchNetwork()
     dev.off()
     HTMLInsertGraph(netf ,file=htmlf, Align = "left", WidthHTML = 1000)
    } else { 
     HTML(paste0("Not completed."),file = htmlf)#,HR=10)
    }
  }

#########
###RMP###
#########

  #Provide random match probabilities
  resRMP <- get("resRMP",envir=nnTK)  #store random match prob results (to be shown in report)
  if(any(optReport$checked[14:15])) { #Any comparisons to show?
    HTML(as.title("Random match probabilities"),file = htmlf,HR=1)
  }
  if(optReport$checked[14]=="TRUE") { #Show RMNE probs?
    if(!is.null(resRMP)) { #if random match probabilities are calculated
     if(nrow(resRMP$evid)>0) {
      HTML(as.title("RMNE for evidence profiles"),file = htmlf,HR=2)
      HTML(as.data.frame(resRMP$evid), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
     }
    } else { 
     HTML(paste0("Not completed."),file = htmlf)#,HR=10)
    }
  }
  if(optReport$checked[15]=="TRUE") { #Show RMP probs?
    if(!is.null(resRMP)) { #if random match probabilities are calculated
     if(nrow(resRMP$ref)>0) {
      HTML(as.title("RMP for reference profiles"),file = htmlf,HR=2)
      HTML(as.data.frame(resRMP$ref), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
     }
    } else { 
     HTML(paste0("Not completed."),file = htmlf)#,HR=10)
    }
  }
  resIBS <- get("resIBS",envir=nnTK)  #store random match prob results (to be shown in report)
  if(optReport$checked[16]=="TRUE") { #Show IBS table?
    HTML(as.title("Related reference candidates"),file = htmlf,HR=2)
    if(!is.null(resIBS)) { #if random match probabilities are calculated
     if(nrow(resIBS)>0) {
      HTML(as.data.frame(resIBS), file= htmlf,row.names=TRUE,align="left",innerBorder = brdt1,Border=brdt2)
     }
    } else { 
     HTML(paste0("Not completed."),file = htmlf)#,HR=10)
    }
  }


##############
###Settings###
##############
  kit0 <- get("setupKit",envir=nnTK)$kitname
  printEPG <- canPrintEPG()  && !is.null(mixLIST) #possible to print EPG?
  if(optReport$checked[17]=="TRUE") { #should settings be shown?
    HTML(as.title(paste0("Settings:")),file = htmlf,HR=1)
    HTML(as.title(paste0("Thresholds:")),file = htmlf,HR=2)
    HTML( paste0("MAC threshold=",get("setupThresh",envir=nnTK)$MACthresh) ,file = htmlf)#,HR=0)
    HTML( paste0("Qual. LR threshold=",get("setupThresh",envir=nnTK)$LRthresh1) ,file = htmlf)#,HR=0)
    HTML( paste0("Quan. LR threshold=",get("setupThresh",envir=nnTK)$LRthresh2) ,file = htmlf)#,HR=0)
    HTML( paste0("Minimum loci=",get("setupThresh",envir=nnTK)$minLociSS) ,file = htmlf)#,HR=0)
    HTML( paste0("Deconv.ratio=",get("setupThresh",envir=nnTK)$ratio) ,file = htmlf)#,HR=0)
    HTML( paste0("Deconv.alleleProb=",get("setupThresh",envir=nnTK)$probA) ,file = htmlf)#,HR=0)

    HTML(as.title(paste0("Model:")),file = htmlf,HR=2)
    HTML( paste0("Frequency file=",get("setupPop",envir=nnTK)$popfile) ,file = htmlf)#,HR=0)

    modeltype=get("setupModel",envir=nnTK)$modeltype
    if(modeltype%in%c(1,3)) {
     HTML(as.title(paste0("Qualitative (LRmix):")),file = htmlf,HR=3)
     HTML( paste0("Drop-in prob=",get("setupModel",envir=nnTK)$dropinC) ,file = htmlf)#,HR=0)
    } 
    if(modeltype%in%c(2,3)) {
     HTML(as.title(paste0("Quantitative (EuroForMix):")),file = htmlf,HR=3)
     HTML( paste0("Detection threshold=",get("setupModel",envir=nnTK)$threshT) ,file = htmlf)#,HR=0)
     HTML( paste0("Kit=",kit0) ,file = htmlf)#,HR=0)
     HTML( paste0("Degradation model=",get("setupModel",envir=nnTK)$degrad) ,file = htmlf)#,HR=0)
     HTML( paste0("Stutter model=",get("setupModel",envir=nnTK)$stutt) ,file = htmlf)#,HR=0)
     HTML( paste0("Drop-in prob=",get("setupModel",envir=nnTK)$dropinC) ,file = htmlf)#,HR=0)
     HTML( paste0("Drop-in Lambda=",get("setupModel",envir=nnTK)$dropinL) ,file = htmlf)#,HR=0)
    } 
  }

################
###ATTACHMENT###
################
  if(any(optReport$checked[18:19])) { #Any comparisons to show?
     HTML(as.title("Attachments"),file = htmlf,HR=1)
  }
  whsize <- c(1920*2,1080*2,120*2) #number of pixels and resolution

  #Plot EPG figures for single sources
  if(printEPG && optReport$checked[18]=="TRUE") { 
      HTML(as.title(paste0("Single source profiles")),file = htmlf,HR=2)
      if(nrow(ssDataTABLE)>0) {
       unREF = unique(ssDataTABLE[,1]) #get unique refs
       refL = getRefL(unREF,forceDi=FALSE) # get relevant references

       for(i in 1:nrow(ssDataTABLE)) { #for each single source profiles
        evid <- rownames(ssDataTABLE)[i]
        ref <- ssDataTABLE[i,1]

        if(evid=="empty") {
         HTML("Empty.",file = htmlf)
        } else {
         condref <- refL[ref] #extract reference
         if(length(condref)==0) condref = NULL

         epgf <- file.path(path2,paste0("epg_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
         png(epgf ,width =whsize[1], height = whsize[2],res=whsize[3])
         tryCatch( euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE) ,error = function(e) euroformix::plotLUS(mixLIST[[evid]],sn=evid,condref))
         dev.off()
         HTMLInsertGraph(epgf,file=htmlf, Align = "left", WidthHTML =  whsize[2]*0.8)
        } #end if not empty
         HTML(as.title(paste0("#",i," - ",evid)),file = htmlf,HR=3)
       } #end for each samples
      } else {
        HTML(paste0("None."),file = htmlf)#,HR=10)
      }#end if
  } #end if plot EPG

  #Plot EPG figures for mixtures?
  if(printEPG && optReport$checked[19]=="TRUE") { 
      HTML(as.title(paste0("Mixture profiles")),file = htmlf,HR=2)
      if(nrow(mixDataTABLE)>0) {
       matches <- get("resMatches",envir=nnTK)  #get match list of mixtures
       if(!is.null(matches) && nrow(matches)>0) { #require match table
        unREF = unique(unlist(strsplit(matches[,2],"/"))) #get unique refs
        refL = getRefL(unREF,forceDi=FALSE) # get relevant references
       }

       for(i in 1:nrow(mixDataTABLE)) { #for each single source profiles
        evid <- rownames(mixDataTABLE)[i]
        condref = NULL
        if(!is.null(matches) && nrow(matches)>0) { #require match table
         ind <- matches[,1]%in%evid
         if(any(ind)) {
           refs <- unlist(strsplit(matches[ind,2],"/")) #get refs
           condref <- refL[refs]
           if(length(condref)==0) condref = NULL
         }
        }
        epgf <- file.path(path2,paste0("epg_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
        png(epgf ,width =whsize[1], height = whsize[2],res=whsize[3])
        tryCatch( euroformix::plotEPG(mixLIST[evid],refcond=condref,kitname=kit0, showPH = TRUE) ,error = function(e) euroformix::plotLUS(mixLIST[[evid]],sn=evid,condref))
        dev.off()
        HTMLInsertGraph(epgf,file=htmlf, Align = "left", WidthHTML = whsize[2]*0.8)
        HTML(as.title(paste0("#",i," - ",evid)),file = htmlf,HR=3)
       } #end for each samples
      } else {
       HTML(paste0("None."),file = htmlf)#,HR=10)
      }#end if
  } #end if plot EPG

 #Plot expected PH for matching references (+DC):
  if(printEPG && optReport$checked[21]=="TRUE") { 
      mixList= get("allMixList",envir=nnTK) #get mixture list (MIXTURES)
      HTML(as.title(paste0("Expected PH for match candidates")),file = htmlf,HR=2)
      if(is.null(mixList) || all(mixList[,2]=="")) {
        HTML(paste0("None."),file = htmlf)#,HR=10)
      } else {
       	#Go through each evidence with a match and show
       	nS = nrow(mixList) #get number of mixtures 
		if(nS==0) return()
		candtabs = numeric()
    		for(ss in 1:nS) {
			refs = unlist(strsplit(mixList[ss,2],"/"))
			if(length(refs)==0) next 
			evid = mixList[ss,1]
			nC = as.integer(mixList[ss,3])

		  tryCatch(	{ 
   			epgf <- file.path(path2,paste0("exp_",gsub(.Platform$file.sep,"_",evid),".png")) #file of picture
			png(epgf ,width =whsize[1], height = whsize[2],res=whsize[3])
			candtab <- doDC(nC,evid,refs,showPlot=TRUE)  
			dev.off()
			HTMLInsertGraph(epgf,file=htmlf, Align = "left", WidthHTML = whsize[2]*0.8)
			HTML(as.title(paste0("#",ss," - ",evid)),file = htmlf,HR=3)
			candtabs = rbind(candtabs,candtab) #add candidate table
		  } )
		  print(paste0(round(ss/nS* 100), "% calculation complete...")) 
		} #end for each mixtures
		if(length(candtabs)>0) { #if more than one new after DC 
  			dclist <- rbind(get("DClist",envir=nnTK),candtabs) #add candidates
               dclist  <- unique(dclist) #consider only the unique variants
			assign("DClist",dclist,envir=nnTK) #get stored DC-list  
			refreshDCLIST() #refresh DC-list
   			gWidgets2::svalue(nb) <- 6 #go to DC-tab after calculation
		} #end if any candidate DC (added to list)
	} #end if else
  } #end plot EPG
  print("REPORT SAVED:")
  print(htmlf)
  browseURL(htmlf) #look directly in browser after creating report
} #end create report


##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 ###############
 #start layout:#
 ###############
 mblst = list( #NOTICE THE NEW CODE IN gWidgets2
  File=list(  
   gWidgets2::gaction('Set directory',handler=f_setwd),
   gWidgets2::gaction('Open project',handler=f_openproj),
   gWidgets2::gaction('Save project',handler=f_saveproj),
   gWidgets2::gaction('Quit',handler=f_quitproj,icon="close")
  ),
  Setup=list( #The values of these are stored only ones
   gWidgets2::gaction('Threshold settings',handler=f_threshsel),
   gWidgets2::gaction('Model settings',handler=f_modelsel),
   gWidgets2::gaction('Kit selection',handler=f_kitsel),
   gWidgets2::gaction('Population frequencies',handler=f_popsel),
   gWidgets2::gaction('Set importData function',handler=f_importsel),
   gWidgets2::gaction('Set case directory',handler=f_casedirsel)
  ),
  Report=list( #The values of these are stored only ones
   gWidgets2::gaction('Report layout',handler=f_reportlay)
  ),
  Advanced=list( #The values of these are stored only ones
   gWidgets2::gaction('Advanced options',handler=f_options),
   gWidgets2::gaction('Random IBS',handler=f_calcIBSdist) #calculates random IBS
  )
 )

 #change working directory to the one stored in nnTK-environment
 wd=get("workdir",envir=nnTK) #assign working directory to nnTK-environment
 if(!is.null(wd)) setwd(wd)
 
 #Main window: 
 mainwin <- gWidgets2::gwindow(softname, visible=TRUE, width=mwW,height=mwH)
 gWidgets2::addHandlerDestroy( mainwin, handler = f_quitproj  ) #call quit function

 gWidgets2::gmenu(mblst,container=mainwin)
 nb = gWidgets2::gnotebook(container=mainwin)
 tabimport = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Data") #tab2: (imports all files)
 tabmatchmatrix = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Match matrix") #comparison results
 tabmatchlist1 = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Match list (Qual LR)") #results from Qualitative model (LRmix)
 tabmatchlist2 = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Match list (Quan LR)") #results from Quantitative model (EFM)
 tabmixtures = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Mixtures") #comparison results (collapsed)
 tabdeconv = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Deconvoluted") #comparison results


####################################################
###############Tab 1: Import Data:##################
####################################################

 #TAB layout
 layhor0 = as.logical(get("setupView",envir=nnTK)$importHorizontal) #get configured layout of tables
 tabimportA = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe("Select Case ID",container=tabimport)) #kit and population selecter
 tabimportB = gWidgets2::ggroup(horizontal=layhor0,container=gWidgets2::gframe( paste0("Evidence profile(s) and Reference profile(s)"),container=tabimport,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
 tabimportC = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe("Functionalities",container=tabimport)) #Tasks button

 #Choose box and import button
 casedir <-  get("setupCase",envir=nnTK)$casepath
 if(casedir=="") gWidgets2::gmessage("The user must specify the directory of the Case folders.\nGo to Settings->Case directory")
 casefolds <- list.dirs(casedir,recursive=FALSE)
 casefolds2 <- sapply(casefolds, function(x) tail(unlist(strsplit(x,.sep)),1))
 
 caseid <- get("caseID",envir=nnTK) 
 if(is.null(caseid)) {
  caseid <- 0
 } else {
  caseid <- which(casefolds2==caseid) #must be index
 }
 tabimportA[1,1] = gWidgets2::gbutton(text="Change view",container=tabimportA,handler=function(h,...) {
  opt = get("setupView",envir=nnTK) #get options
  if(layhor0) {
    opt$importHorizontal = "FALSE"
  } else {
    opt$importHorizontal = "TRUE"
  }
  setupWrite(unlist(opt),file=setupFileView)    #save to file in installation folder
  gWidgets2::dispose(mainwin) #shut down window
  gc() #empty garbage memory
  gui() #start an empty session
 })
 tabimportA[1,2] <- gWidgets2::gcombobox(items=casefolds2, width=100, selected =caseid  , editable = TRUE, container = tabimportA)

 tabimportA[1,3] = gWidgets2::gbutton(text="Import",container=tabimportA,handler=f_importData)
 mixTabGUI <- gWidgets2::gtable(items="",multiple = TRUE,container = tabimportB, handler=clicktable,action="mix")
 refTabGUI <- gWidgets2::gtable(items="",multiple = TRUE,container = tabimportB, handler=clicktable,action="ref")
 gWidgets2::add(tabimportB,mixTabGUI,expand=TRUE,fill=TRUE)
 gWidgets2::add(tabimportB,refTabGUI,expand=TRUE,fill=TRUE)

 sortstring1 <- c("Sort Evids by #","Sort Evids by name")
 tabimportA[1,4] = gWidgets2::gradio(items=sortstring1,selected=1,container=tabimportA,horizontal = TRUE,expand=TRUE,fill=TRUE,handler=
  function(h,...) { 
   updateTables(sortmix=which(sortstring1==svalue(tabimportA[1,4])),sortref=which(sortstring2==svalue(tabimportA[1,5]))) 
 })
 sortstring2 <- c("Sort Refs by #","Sort Refs by name")
 tabimportA[1,5] = gWidgets2::gradio(items=sortstring2,selected=1,container=tabimportA,horizontal = TRUE,expand=TRUE,fill=TRUE,handler=
  function(h,...) { 
   updateTables(sortmix=which(sortstring1==svalue(tabimportA[1,4])),sortref=which(sortstring2==svalue(tabimportA[1,5]))) 
 })


 tabimportC[1,1] = gWidgets2::gbutton(text="Compare",container=tabimportC,handler=
  function(h,...) {

  if(nrow(get("refDataTABLE",envir=nnTK))==0) return() #return if no refs

  #CHECK AND SET FREQUENCY FILE
  ok = setPopFreq() #import population frequency from last selected file
  if(!ok) return() 

  #Step 1: Calculate MAC
  res <- getMatchesMAC(locs=names(get("popFreq",envir=nnTK)))  #get name of loci to consider
  if(is.null(res)) return(); #return if nothing to compare.
  refreshTabMATRIX() #update tables
  gWidgets2::svalue(nb) <- 2 #go to comparison tab
  
  #CALCULATING LR BASED SCORES:
  if(nrow(res$MatchList)==0) return() #Return if no candidates to calculate
  modtype <- get("setupModel",envir=nnTK)$modeltype #model type selected {1="qual",2="quan",3="both"}
  #Step 2 (optional): Calculate qual based LR
  if( modtype%in%c(1,3) ) {
   getMatchesLR(type="qual") #LRmix
   refreshTabLIST1() #update table with results
   gWidgets2::svalue(nb) <- 3 #go to qual LR result tab when done
  }

  #Step 3 (optional): Calculate quan based LR
  if( modtype%in%c(2,3) ) {
   if( get("setupModel",envir=nnTK)$degrad=="ON" && !canPrintEPG() ) {
    gWidgets2::gmessage("Please specify a kit for further comparisons.\nGo to Settings->Kit selection")
    return()
   } 
   getMatchesLR(type="quan") #EFM 
   refreshTabLIST2() #update tables with results
   gWidgets2::svalue(nb) <- 4 #go to quan LR result tab when done
  }

  #Step 3: Create matchlist (Final)
  createMatchlist(modtype=get("setupModel",envir=nnTK)$modeltype)
  refreshTabLIST() #update tables  
  gWidgets2::svalue(nb) <- 5 #go to overview when done

 })
 tabimportC[1,2] = gWidgets2::gbutton(text="Create Report",container=tabimportC,handler=f_createreport)

 tabimportC[1,3] = gWidgets2::gbutton(text="Add/Edit references",container=tabimportC,handler=f_addref)

 tabimportC[1,4] = gWidgets2::gbutton(text="Export profile(s)",container=tabimportC,handler=f_exportData)

 tabimportC[1,5] = gWidgets2::gbutton(text="Calculate RMP",container=tabimportC,handler=f_calcRMP)

 tabimportC[1,6] = gWidgets2::gbutton(text="Calculate IBS",container=tabimportC,handler=f_calcIBS)

 tabimportC[1,7] = gWidgets2::gbutton(text="Deconvolve",container=tabimportC,handler=function(h,...){
  last <- get("clicktableLast",envir=nnTK) #assign to nnTK-environment
  evids <- rownames(get("mixDataTABLE",envir=nnTK))[last$mix] #get name
  refs <- rownames(get("refDataTABLE",envir=nnTK))[last$ref] #get name
  if(length(evids)==0) return() 
  createHypDCWindow(evids,refs,nC=1) #create window for specify model for DC  
 }) #end button 


 tabimportC[1,8] = gWidgets2::gbutton(text="RESTART",container=tabimportC,handler=
  function(h,...) {
   gWidgets2::dispose(mainwin) #shut down window
   gc() #empty garbage memory
   gui() #start an empty session
 })

 #INSERT DATA (TABLE-FORMAT) TO GUI: NOTICE the clicktable handler 
 updateTables <- function(sortmix=1,sortref=1) { #function to call to update tables  
   mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #assign to nnTK-environment
   refDataTABLE <- get("refDataTABLE",envir=nnTK) #assign to nnTK-environment
   mixDataMATCHSTATUS <- get("mixDataMATCHSTATUS",envir=nnTK) #assign to nnTK-environment
  
   mixTab <- refTab <- "" #empty
   if(!is.null(mixDataTABLE)) {    #Add mix-table
    newtab <- cbind(mixDataMATCHSTATUS,mixDataTABLE)
    colnames(newtab)[1] <- "MatchStatus"
    mixTab <- addRownameTable(newtab,type=2)
    if(sortmix==2)  mixTab <- mixTab[order(mixTab[,2],decreasing=FALSE),,drop=FALSE]
   } 
   if(!is.null(refDataTABLE)) {    #Add mix-table
    refTab <- addRownameTable(refDataTABLE,type=2)
    if(sortref==2)  refTab <- refTab[order(refTab[,2],decreasing=FALSE),,drop=FALSE]
   }
   mixTabGUI[] <- mixTab
   refTabGUI[] <- refTab

   #ADJUST COLUMN WIDTH:
   if(!layhor0 && !is.null(ncol(mixTab)) && !is.null(ncol(refTab)) ) { #if vertical layout:
    colw1 = c(30,150,150)
    colw2 = c(colw1[1],sum(colw1[-1]))
    colL = 100 #column widt for each locus

    gWidgets2::size(mixTabGUI) <- list(column.widths=c(colw1,rep(colL,ncol(mixTab)-length(colw1))))
    gWidgets2::size(refTabGUI) <- list(column.widths=c(colw2,rep(colL,ncol(refTab)-length(colw2))))
   }
 }
 updateTables() #update when program starts

####################################################################################################################
#######################################Tab 2: Match matrix: #############################################################
#####################################################################################################################

 f_rotateMatchMatrix = function(h,...) {
  opt = get("setupView",envir=nnTK) #get options
  if(opt$matchmatrix == "TRUE") {
    opt$matchmatrix = "FALSE"
  } else {
    opt$matchmatrix = "TRUE"
  }
  setupWrite(unlist(opt),file=setupFileView)    #save to file in installation folder
  assign("setupView",opt,envir=nnTK) #store again
  refreshTabMATRIX(rotate=TRUE) #update with rotated table  
 }

 f_truncatevals = function(h,...) { #removes values below threshold table
  val = get("setupThresh",envir=nnTK)$MACthresh #get options
  matchmat <- get("resCompMAC",envir=nnTK)$MatchMatrix
  if(svalue(tabCompA0[1,3])=="Truncate")  {
   matchmat[matchmat<val] = ""   #remove values below threshold
  }
  matchMATGUI[] <-  addRownameTable(matchmat,type=2)
  if(!is.null(ncol(matchMATGUI[]))) gWidgets2::size(matchMATGUI) <- list(column.widths=c(30,rep(100,ncol(matchMATGUI[])-1)))
 }

 tabCompA0 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe("Further",container=tabmatchmatrix)) #kit and population selecter
 tabCompA0[1,1] <- gWidgets2::gbutton(text="Export",container=tabCompA0,handler=f_exporttable,action="mac")  #Button to create matchcloud
 tabCompA0[1,2] <- gWidgets2::gbutton(text="Change view",container=tabCompA0,handler=f_rotateMatchMatrix)  #Button to rotate table
 tabCompA0[1,3] <- gWidgets2::gradio(items=c("No truncation","Truncate"),selected=1,container=tabCompA0,handler=f_truncatevals)  #Button to create matchcloud
 tabCompMatrix = gWidgets2::ggroup(container=tabmatchmatrix,expand=TRUE,fill=TRUE)
 matchMATGUI <- gWidgets2::gtable(items="",container=tabCompMatrix) #add to frame
 gWidgets2::add(tabCompMatrix,matchMATGUI,expand=TRUE,fill=TRUE)
 refreshTabMATRIX = function(rotate=FALSE) { 
    resobj = get("resCompMAC",envir=nnTK)
    if(is.null(resobj$MatchMatrix)) return() #return if no match matrix found
    if(rotate) {
     resobj$MatchMatrix = t(resobj$MatchMatrix)
     assign("resCompMAC",resobj,envir=nnTK) #store rotated object (will be rotated in report)
    }
    matchMATGUI[] <-  addRownameTable(resobj$MatchMatrix,type=2) #rotate view
    if(!is.null(ncol(matchMATGUI[]))) gWidgets2::size(matchMATGUI) <- list(column.widths=c(30,rep(100,ncol(matchMATGUI[])-1)))
 }
 refreshTabMATRIX()

####################################################################################################################
#######################################Tab 3a: Match List (Qual): ##################################################
####################################################################################################################

 tabCompA1 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe("Further",container=tabmatchlist1)) #kit and population selecter
 tabCompA1[1,1] <- gWidgets2::gbutton(text="Export",container=tabCompA1,handler=f_exporttable,action="qual")  
 tabCompA1[1,2] <- gWidgets2::gbutton(text="Calculate Quan LRs",container=tabCompA1,handler=f_calcLRefm)  
 tabCompLIST1 = gWidgets2::ggroup(container=tabmatchlist1,expand=TRUE,fill=TRUE)

 matchL1GUI <- gWidgets2::gtable(items="",container=tabCompLIST1,handler=clickmatchlistQUAL)
 gWidgets2::add(tabCompLIST1,matchL1GUI,expand=TRUE,fill=TRUE)#add to frame

 refreshTabLIST1 = function() { 
    matchlist <- get("resCompLR1",envir=nnTK)
    if(is.null(matchlist)) return() #return if no list found
    if(nrow(matchlist)==0) return() #return if no candidate found
    matchL1GUI[] <-  addRownameTable(matchlist,type=3)
    gWidgets2::size(matchL1GUI) <- list(column.widths=c(30,300,300,70,70,100))
 }
 refreshTabLIST1()


####################################################################################################################
#######################################Tab 3b: Match List (Quan): ##################################################
####################################################################################################################

 tabCompA2 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe("Further",container=tabmatchlist2)) #kit and population selecter
 tabCompA2[1,1] <- gWidgets2::gbutton(text="Export",container=tabCompA2,handler=f_exporttable,action="quan")#  
 tabCompLIST2 = gWidgets2::ggroup(container=tabmatchlist2,expand=TRUE,fill=TRUE)

 matchL2GUI <- gWidgets2::gtable(items="",container=tabCompLIST2,handler=clickmatchlistQUAN)
 gWidgets2::add(tabCompLIST2,matchL2GUI,expand=TRUE,fill=TRUE)#add to frame

 refreshTabLIST2 = function(selInd=NULL) { 
    matchlist <- get("resCompLR2",envir=nnTK)
    if(is.null(matchlist)) return() #return if no list found
    if(nrow(matchlist)==0) return() #return if no candidate found
    matchL2GUI[] <-  addRownameTable(matchlist,type=3)
    gWidgets2::size(matchL2GUI) <- list(column.widths=c(30,300,300,70,70,rep(100,ncol(matchlist)-4)))
    if(!is.null(selInd)) gWidgets2::svalue(matchL2GUI) <- selInd #select row
 }
 refreshTabLIST2()

################################################################################################
#################################Tab 4: Mixtures: ##############################################
################################################################################################
 tabCompA3 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe("Further",container=tabmixtures)) #kit and population selecter
 tabCompA3[1,1] <- gWidgets2::gbutton(text="Export",container=tabCompA3,handler=f_exporttable,action="final")#  
 tabCompA3[1,2] <- gWidgets2::gbutton(text="Show match network",container=tabCompA3,handler=showMatchNetwork)  #Button to create matchcloud
 tabCompA3[1,3] <- gWidgets2::gbutton(text="Deconvolve all",container=tabCompA3,handler=f_doDCall)# perform DC for all mixtures (with default settings)
 tabmixLIST = gWidgets2::ggroup(container=tabmixtures,expand=TRUE,fill=TRUE)
 mixlistGUI <- gWidgets2::gtable(items="",container=tabmixLIST,handler=clickmixlist)
 gWidgets2::add(tabmixLIST,mixlistGUI,expand=TRUE,fill=TRUE)#add to frame

 refreshTabLIST = function(evidsel=NULL) {  #get list with all mixtures with info about matched elements 
    #store match-results from comparison (those with LR>threshold) together with all mixtures
    mixN <- get("mixDataMATCHSTATUS",envir=nnTK)
    mixN <- names(mixN)[mixN=="mixture"]  #get sampleName of all mixtures
    if( length(mixN)==0 ) return() #return if no mixtures existing   
    nClow <- sapply(get("mixDataLIST",envir=nnTK)[mixN],function(elem) { ceiling(max(sapply(elem,function(x) length(x$adata)))/2) })
    mixlist <- cbind(mixN,"",nClow)
    colnames(mixlist) <- c("Evidence","Reference(s)","Num Contr.")
   
    #import matchlist and insert reference candidates
    matchlist <- get("resMatches",envir=nnTK) #get match list (from function)
    mixlist[match(matchlist[,1],mixN),2] <- matchlist[,2] #Update contributing References
    mixlist[match(matchlist[,1],mixN),3] <- matchlist[,3] #Update number of contributors
    assign("allMixList",mixlist,envir=nnTK)  #store match-results from comparison (those with LR>threshold) together with all mixtures

    mixlistGUI[] <- addRownameTable(mixlist,type=3)
    gWidgets2::size(mixlistGUI) <- list(column.widths=c(30,300,500,100))
    if(!is.null(evidsel)) gWidgets2::svalue(mixlistGUI) = which(mixlist[,1]==evidsel)  #marking line
  }
  refreshTabLIST()

################################################################################################
#################################Tab 5: Deconvoluted: ##############################################
################################################################################################
 tabCompA4 = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=gWidgets2::gframe("Further",container=tabdeconv)) #kit and population selecter
 tabCompA4[1,1] <- gWidgets2::gbutton(text="Export",container=tabCompA4,handler=f_exporttable,action="dc")#  

 tabDCLIST = gWidgets2::ggroup(container=tabdeconv,expand=TRUE,fill=TRUE)
 DClistGUI <- gWidgets2::gtable(items="",container=tabDCLIST,handler=clickDClist)
 gWidgets2::add(tabDCLIST,DClistGUI,expand=TRUE,fill=TRUE)#add to frame

 refreshDCLIST = function() {  #get list of matched elements
   DClist <- get("DClist",envir=nnTK)
   if(is.null(DClist)) return() #return if no list found
   DClistGUI[] <- addRownameTable( DClist ,type=3) #add DC-list
   gWidgets2::size(DClistGUI) <- list(column.widths=c(30,100,100,30,50,rep(50,ncol(DClistGUI)-5))) 
 }
 refreshDCLIST()

#####################################################################################################

 gWidgets2::visible(mainwin) <- TRUE
 gWidgets2::focus(mainwin)
 gWidgets2::svalue(nb) <- 1 #initial start at second tab

} #end funcions
