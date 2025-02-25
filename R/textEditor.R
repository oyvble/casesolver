#' @title textEditor
#' @description User select items between two tables 
#' @details The function returns the items selected in the two lists (when window is closed)
#' @param env an environment with elements items1 and items2
#' @export
#' @examples
#' \dontrun{
#' env = new.env( parent = emptyenv() ) 
#' assign("text","this is my text",envir=env)
#' textEditor(env)
#' }

textEditor = function(env) {
  objStore = "text"
  text0 = get(objStore,envir=env)

  #event when clicking OK
  f_clickbut = function(h) {
    bool = TRUE #gWidgets2::gconfirm("Are you sure you want to continue?")
    if(!bool) return() #
    
    
    if(h$action=="restore") {
      gWidgets2::svalue(textGUI) <- text0
    }
    if(h$action=="save") {
      assign(objStore, gWidgets2::svalue(textGUI) ,envir=env)
      gWidgets2::dispose(inwin) #close window
      tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    }
    if(h$action=="quit") {
      gWidgets2::dispose(inwin) #close window
      tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    }
  }
  
  flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
  
  #inwin = gWidgets2::gwindow("Text editor",visible=FALSE,expand=TRUE,fill=TRUE,width=500,height = 200)
  inwin = gWidgets2::gwindow("Text editor",visible=FALSE,width=500,height = 200)
  gWidgets2::addHandlerUnrealize( inwin, handler = function(h,...) { #
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    return(NULL) #return selected profiles when window is exit 
  }  ) #call quit function
  
  grp = gWidgets2::ggroup(horizontal = FALSE,container=inwin,expand=TRUE,fill=TRUE)
  butLay = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=grp) #kit and population selecter
  butLay[1,1] = gWidgets2::gbutton("Restore",container=butLay, action="restore", handler=f_clickbut)
  butLay[1,2] = gWidgets2::gbutton("Save",container=butLay, action="save", handler=f_clickbut)
  butLay[1,3] = gWidgets2::gbutton("Quit",container=butLay, action="quit", handler=f_clickbut)
  textGUI <- gWidgets2::gtext(text0, container=grp,width = 300,height=150,wrap = TRUE,fill=TRUE,expand=TRUE) #editor
  gWidgets2::visible(inwin) = TRUE #NOTE HAVING ISSUE TO SEE THE EDITOR

  tcltk::tkwait.variable(flag) #important to not quit window before broken
} #end function



