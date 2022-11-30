#' @title profileSwapperGUI
#' @description User select items between two tables 
#' @details The function returns the items selected in the two lists (when window is closed)
#' @param env an environment with elements items1 and items2
#' @export
#' @examples
#' \dontrun{
#' env = new.env( parent = emptyenv() ) 
#' assign("itemList",list(items1=paste("item",1:10),items2=paste("item",21:100)),envir=env)
#' profileSwapperGUI(env)
#' }

profileSwapperGUI = function(env) {
  objStore = "itemList"
  itemList = get(objStore,envir=env)
  width0=150 #width of table
  
  #event when clicking OK
  f_clickOK = function(h) {
    #bool = gWidgets2::gconfirm("Are you sure you want to continue?")
    #if(!bool) return() #
    assign(objStore,list(items1=tab1GUI[],items2=tab2GUI[]),envir=env)
    gWidgets2::dispose(inwin) #close window
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
  }
  
 #event when clicking swap button (right=1 or left=2)
  f_swap = function(h) {
    items1=tab1GUI[]
    items2=tab2GUI[]
    
    #right --->  
    if(h$action==1) {  
      selItem = gWidgets2::svalue(tab1GUI)
      items2 <- union(items2,selItem)
      items1 <- setdiff(items1,selItem)
    }
    
    #left <----  
    if(h$action==2) {  
      selItem = gWidgets2::svalue(tab2GUI)
      items1 <- union(items1,selItem)
      items2 <- setdiff(items2,selItem)
    }
    tab1GUI[] <- items1
    tab2GUI[] <- items2
    gWidgets2::size(tab1GUI) <- gWidgets2::size(tab2GUI) <- list(column.widths=width0) 
  }
  
  flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
  inwin = gWidgets2::gwindow("Profile selection",visible=FALSE,expand=TRUE,fill=TRUE,width=500)
  gWidgets2::addHandlerUnrealize( inwin, handler = function(h,...) { #
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    return(NULL) #return selected profiles when window is exit 
  }  ) #call quit function

  #header <- ggroup(horizontal = TRUE, container=inwin,expand=TRUE,fill=TRUE)
  #glabel("Included", container=header)  
  #glabel("Swap", container=header)  
  #glabel("Excluded", container=header)  

  frame <- gWidgets2::ggroup(horizontal = TRUE, container=inwin,expand=TRUE,fill=TRUE)
  suppressWarnings({
    tab1GUI <- gWidgets2::gtable(items=itemList$items1,multiple = TRUE,container = frame,expand=TRUE,fill=TRUE, handler=f_swap,action=1)
    gridBut <- gWidgets2::ggroup(space=5, container=frame,horizontal = FALSE)
    gWidgets2::gbutton("-->",container=gridBut, handler=f_swap, action=1)
    gWidgets2::gbutton("<--",container=gridBut, handler=f_swap, action=2)
    gWidgets2::gbutton("OK",container=gridBut, handler=f_clickOK)
    tab2GUI <- gWidgets2::gtable(items=itemList$items2,multiple = TRUE,container = frame,expand=TRUE,fill=TRUE, handler=f_swap,action=2)
    
    gWidgets2::size(tab1GUI) <- gWidgets2::size(tab2GUI) <- list(column.widths=width0) 
  })
  gWidgets2::visible(inwin) = TRUE
  
  
  tcltk::tkwait.variable(flag) #important to not quit window before broken
} #end function