#' @title setMarkerSettings
#' @description Defining marker settings
#' @details The user defines marker settings by inserting values to table
#' @param nnTK an environment object from stored CaseSolver object
#' @export

#library(casesolver);load("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/CaseSolverDev/DNAxs_DNAStatistX Test data NFI/CSmarkersettings.Rdata")

#creates a table with different settings for each markers (based on popfreq)
setMarkerSettings = function(nnTK) { 
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  
  #Step 1: Obtain and structure info from environment (loci names etc)
  #Used environment object: setupMarkers
  setupMarkers = get("setupMarkers",envir=nnTK) #get already stored marker setup
  model <- get("setupModel",envir=nnTK)
  paramDefault = c(model$threshT, model$dropinC, model$dropinL,model$fst ) #obtain default global values
  nParams = length(paramDefault) #number of parameter types
  paramNames = c(L$Marker, L$analyticalthreshold, L$dropinprob, L$dropinpeakheightlambda , L$fstsetting)
  
  popFreq = get("popFreq",envir=nnTK) #get population freqs
  kit0 = getEnvirKit(nnTK) #get kit
  
  #NEED TO HAVE POPULATION FREQUENCIES SPECIFIED TO PROCEED:
  if(is.null(popFreq)) { 
    gWidgets2::gmessage(L$msg.setPopFreq)
    return()
  }
  
  popLocs = names(popFreq) #locus names from population
  
  #Create new marker setting setup: Insert global param by default 
  setupMarkers2 = list()
  setupMarkers2[[1]] <- popLocs #insert 
  for(cc in 1:length(paramDefault)) setupMarkers2[[cc+1]] <- rep(paramDefault[cc],length(popLocs)) #insert defualt
  
  #If earlier specification (INSERT VALUES ALREADY SPECIFIED, SUBSET ON SPECIFIED MARKERS)
  if(!is.null(setupMarkers)) { 
    commonLocMatchInd = match(setupMarkers2[[1]],setupMarkers[[1]]) #obtain index of common loci
    indVec1 = !is.na(commonLocMatchInd) #index vector to insert
    indVec2 =  na.omit(commonLocMatchInd) #index vector to extract
    #all(setupMarkers2[[1]][indVec1]==setupMarkers[[1]][indVec2])
    for(cc in 2:length(setupMarkers)) setupMarkers2[[cc]][indVec1] <- setupMarkers[[cc]][indVec2] #subset
  }
  setupMarkers = setupMarkers2 #update
  
  #Obtaining table with marker vs color overview (kitdyes)
  kitdyes = NULL #default is no colors (dyes)
  if(!is.null(kit0)) { #If kit specified: Subset on markers existing there
    kitMarkers = toupper(euroformix::getKit(kit0,"MARKER")) #Set upper name
    locUseInd = setupMarkers[[1]]%in%kitMarkers #obtain those to use
    for(cc in 1:length(setupMarkers)) setupMarkers[[cc]] <- setupMarkers[[cc]][locUseInd] #subset

    #Create DYE TABLE
    kitdyes = euroformix::getKit(kit0,"COLOR") #get kitinfo from selected kit
    if(!is.na(kitdyes) && length(kitdyes)>1) { #if kitinformation found
      kitdyes$Marker = toupper(kitdyes$Marker) #ensure upper case
      kitdyes = kitdyes[kitdyes$Marker%in%setupMarkers[[1]],,drop=FALSE  ]
      
      #Change order of markers
      ordMatch <- match(kitdyes$Marker,setupMarkers[[1]])#obtain order
      for(cc in 1:length(setupMarkers)) setupMarkers[[cc]] <- setupMarkers[[cc]][ordMatch] #change order for each vector
    }
  }
  locNames = setupMarkers[[1]] #obtain locus names
  nLocs = length(locNames)
    
  setToDefault = function(init=FALSE,usecommon=FALSE) { #helpfunctions to set values back to default (init is variable declaration)
    for(m in 1:nLocs) { #for each markers
      for(c in 1:length(paramDefault)) { #for each param type
        insval = paramDefault[c] #default value to insert 
        if(!usecommon && !is.null(setupMarkers)) { #if not using common values (from settings) restoring data 
          insval = setupMarkers[[c+1]][m] #obtain value
        } 
        if(init) {
          tabval[m+1,c+1] <- gWidgets2::gedit(text=insval,width=w0,container=tabval) #insert value
        } else {
          gWidgets2::svalue(tabval[m+1,c+1]) <- insval #insert default value
        }
      } #end for each objects       
    } #end for each loci
  }
  
  helpShowWindow = function() {
    gWidgets2::focus(setwin) = TRUE #set top
    gWidgets2::visible(setwin) = TRUE #show window
  }
  
  #CREATE GUI WINDOW:
  flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
  setwin <- gWidgets2::gwindow(paste0("Marker specific settings"),visible=FALSE, width=750, height=500)
  gWidgets2::addHandlerUnrealize( setwin, handler = function(h,...) { #
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    return(NULL) #return selected profiles when window is exit 
  }  ) #call quit function
  
  grouplay <- gWidgets2::ggroup(spacing=2,container=(setwin),horizontal = FALSE, use.scrollwindow = TRUE)  #set group layout

  tabsel = gWidgets2::glayout(spacing=10,container=(grouplay),horizontal = TRUE)  #set grid  (will contain buttons)
  tabsel[1,1] = gWidgets2::gbutton(text="Restore",container=tabsel, handler = function(h,...) {
    gWidgets2::visible(setwin) = FALSE #hide window
    setToDefault(FALSE) #dont init 
    helpShowWindow()
  }) 
  tabsel[1,2] = gWidgets2::gbutton(text="Set to default values",container=tabsel, handler = function(h,...) {
    gWidgets2::visible(setwin) = FALSE #hide window
    setToDefault(FALSE,TRUE) #dont init 
    helpShowWindow()
  }) 
  tabsel[1,3] = gWidgets2::gbutton(text="Fill out dye info",container=tabsel, handler = function(h,...) {
    gWidgets2::visible(setwin) = FALSE #hide window) 
    #filling in values over all dyes (using those specified in  top dye )
    dyes = kitdyes$Color
    unDyes = unique(dyes) #obtain unique colors
    for(c in 1:nParams) { #for each param type
      for(dye in unDyes) {  #for each dye
        rowind = which(dyes==dye) #obtain row indices of selected dyes
        vals = sapply(rowind,function(x) gWidgets2::svalue(tabval[x+1,c+1])) #obtain values
        useval = vals[vals!=""] #obtain value in cells
        if(length(useval)>0) {
          for(ii in rowind) gWidgets2::svalue(tabval[ii+1,c+1]) <- useval[1] #fill with first value only
        }
      }
    }
    helpShowWindow()
  })
  if(is.null(kitdyes)) gWidgets2::enabled(tabsel[1,3]) = FALSE #deactivate if no dye information found
  
  tabsel[1,4] = gWidgets2::gbutton(text="Empty all",container=tabsel, handler = function(h,...) {
    gWidgets2::visible(setwin) = FALSE #hide window
    for(m in 1:nLocs) {
      for(c in 1:nParams) {
        gWidgets2::svalue(tabval[m+1,c+1]) = "" #insert emtpy cell
      }
    }
    helpShowWindow()
  } ) 
  
  #extracting values from table
  tabsel[1,5] = gWidgets2::gbutton(text="Save settings",container=tabsel, handler = function(h,...) {
    for(c in 1:nParams) {
      for(m in 1:nLocs) {
        val = as.numeric( gWidgets2::svalue(tabval[m+1,c+1]) )  #convert from text to number
        if(is.na(val) || val<0) {
          gWidgets2::gmessage("Invalid input found.")
          return()
        }
        setupMarkers[[c+1]][m] = val #insert value
      }
    }  
    assign("setupMarkers",setupMarkers,envir=nnTK)  #store
    gWidgets2::dispose(setwin) #close window when successful
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
  }) 
  
  tabval = gWidgets2::glayout(spacing=0,container=(grouplay))  #create grid layout
  w0 <- 15 #width of textbox
  
  for(cc in 1:length(paramNames)) tabval[1,cc] <- gWidgets2::glabel(text=paramNames[cc],container=tabval) 
  if(!is.null(kitdyes)) tabval[1,length(paramNames)+1] <- gWidgets2::glabel(text="Dye (color)",container=tabval) 
  
  #Insert marker names and possibly dyes
  for(cc in 1:nLocs) {
    tabval[cc+1,1] <- gWidgets2::glabel(text=locNames[cc],container=tabval) #init markers
    if(!is.null(kitdyes))  tabval[cc+1,length(paramNames)+1] <- gWidgets2::glabel(text=kitdyes[cc,2],container=tabval) #insert dye
  }

  setToDefault(init=TRUE) #set to default values
  gWidgets2::visible(setwin) <- TRUE
  tcltk::tkwait.variable(flag) #important to not quit window before broken
} #end gui-function # f_markerSettings()
