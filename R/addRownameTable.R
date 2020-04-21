#' @title addRownameTable
#' @description Helpfunction for CaseSolver to set rownames as first column in table (modification)
#' @param tab The table to modify
#' @param type The type of modification: {0,3}=First column is indices, 1=First column is rownames, 2=1st column is indices AND 2nd column is rownames, 4=1st column is rownames AND collesponding column named "samplename"
#' @param samplename Stringname for "SampleName" (may be in another language)
#' @return A modified table (wrt rownames)
#' @export

 addRownameTable = function(tab,type=1,samplename="SampleName") {
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
   colnames(tab) <- c(samplename,tmp)
  } 
  if(type==2) {
   tab <- cbind(paste0("#",1:nrow(tab)),rownames(tab),tab)
   colnames(tab) <- c(".",samplename,tmp)
  }
  if(type==3) {
   tab <- cbind(paste0("#",1:nrow(tab)),tab)
   colnames(tab) <- c(".",tmp)
  }
  if(type==4) {
   tab <- cbind(rownames(tab),tab)
   colnames(tab) <- c(samplename,tmp)
   rownames(tab) <- 1:nrow(tab)
  } 
  return(NAtoSign(tab))
 }
