#' @title profileSelectorGUI
#' @description A function for letting user selecting profiles
#' @details The function returns the names of the selected profiles (when window is closed)
#' @param env an environment with elements evids and refs
#' @export
#' @examples
#' \dontrun{
#' env = new.env( parent = emptyenv() ) 
#' assign("selected",list(evids=1:10,refs=1:100),envir=env)
#' profileSelectorGUI(env)
#' }

profileSelectorGUI = function(env=NULL) {
  selList = get("selected",envir=env)
  refs = selList$refs
  evids = selList$evids
  
  flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
  inwin = gWidgets2::gwindow("Select specific data",visible=FALSE,expand=TRUE,fill=TRUE)
  gWidgets2::addHandlerUnrealize( inwin, handler = function(h,...) { #
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    return(NULL) #return selected profiles when window is exit 
  }  ) #call quit function
  
  gg = gWidgets2::ggroup(container=inwin,horizontal = FALSE) #set up window layout
  butlay <- gWidgets2::glayout(container=gg,expand=FALSE,fill=FALSE) #set button layout
  butlay[1,1] = gWidgets2::gbutton("Use selected data",container=butlay,handler=function(x) {
    evidSel = gWidgets2::svalue(evidtab) #get selected evid profiles
    refSel = gWidgets2::svalue(reftab)  #get selected ref profiles
    
    bool = gWidgets2::gconfirm("Are you sure you want to continue?")
    if(!bool) return() #
    assign("selected",list(evids=evidSel,refs=refSel),envir=env)
    gWidgets2::dispose(inwin) #close window
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
  })
  
  tablay <- gWidgets2::ggroup(container=gWidgets2::gframe( paste0("Data selection") ,container=gg,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE,horizontal = TRUE) #evidence,ref dataframe
  
  tablayEVID <- gWidgets2::ggroup(space=0,container=tablay,expand=TRUE,fill=TRUE,horizontal = FALSE) #vertical groups
  butlayEVID <- gWidgets2::glayout(container=tablayEVID,expand=FALSE,fill=FALSE) #set buttons
  butlayEVID[1,1] = gWidgets2::gbutton("Deselect all",container=butlayEVID,handler=function(x) {
    gWidgets2::svalue(evidtab) = FALSE  #deselect all
  })
  butlayEVID[1,2] = gWidgets2::gbutton("Select all",container=butlayEVID,handler=function(x) {
    gWidgets2::svalue(evidtab) = evids  #select all
  })
  evidtab = gWidgets2::gcheckboxgroup(evids,checked=TRUE ,container=tablayEVID, use.table = TRUE,expand=TRUE,fill=TRUE)
  
  tablayREF <- gWidgets2::ggroup(container=tablay,expand=TRUE,fill=TRUE,horizontal = FALSE) #vertical groups
  butlayREF <- gWidgets2::glayout(container=tablayREF,expand=FALSE,fill=FALSE) #set buttons
  butlayREF[1,1] = gWidgets2::gbutton("Deselect all",container=butlayREF,handler=function(x) {
    gWidgets2::svalue(reftab) = FALSE  #deselect all
  })
  butlayREF[1,2] = gWidgets2::gbutton("Select all",container=butlayREF,handler=function(x) {
    gWidgets2::svalue(reftab) = refs  #select all
  })
  reftab = gWidgets2::gcheckboxgroup(refs,checked=TRUE ,container=tablayREF, use.table = TRUE,expand=TRUE,fill=TRUE)
  
  gWidgets2::visible(inwin) = TRUE
  tcltk::tkwait.variable(flag) #important to not quit window before broken
} #end function