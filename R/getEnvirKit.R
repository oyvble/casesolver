#' @title getEnvirKit
#' @description Obtain defined kit from environment
#' @param nnTK an environment object from stored CaseSolver object
#' @return Name of kit if found, otherwise NULL is returned.
#' @export

getEnvirKit = function(nnTK) { 
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  kitname = get("setupKit",envir=nnTK)$kitname 
  if(is.null(kitname) || is.na(kitname[1]) || kitname==L$none || kitname=="" ) kitname = NULL
  return(kitname)
}

