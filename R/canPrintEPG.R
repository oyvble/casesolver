#' @title canPrintEPG
#' @description A helpfunction to check if EPGs can be printed (kit specified)
#' @param nnTK an environment object from stored CaseSolver object
#' @export

canPrintEPG = function(nnTK) { 
  return( !is.null(casesolver::getEnvirKit(nnTK)) )  #possible to print EPG?
}