#' @title addRownameTable
#' @description Helpfunction for CaseSolver
#' @export

 addRownameTable = function(tab,type=1) {
  NAtoSign <- function(x) {
   x[is.na(x)] <- "-" #NB: New version of gtable does not accept NA values
   return(x)
  }
  if( nrow(tab)==0 ) return(NULL)
  tmp <- colnames(tab)
  if(type==0) { #dont add rowname but put index as rowname
   rownames(tab) <- 1:nrow(tab)
  } 
  if(type==1) {
   tab <- cbind(rownames(tab),tab)
   colnames(tab) <- c("SampleName",tmp)
  } 
  if(type==2) {
   tab <- cbind(paste0("#",1:nrow(tab)),rownames(tab),tab)
   colnames(tab) <- c("ID","SampleName",tmp)
  }
  if(type==3) {
   tab <- cbind(paste0("#",1:nrow(tab)),tab)
   colnames(tab) <- c(".",tmp)
  }
  if(type==4) {
   tab <- cbind(rownames(tab),tab)
   colnames(tab) <- c("SampleName",tmp)
   rownames(tab) <- 1:nrow(tab)
  } 
  return(NAtoSign(tab))
 }
