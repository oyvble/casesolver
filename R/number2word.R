#' @title number2word
#' @description Convert number to verbal
#' @param x Number (on log10 scale)
#' @param L Language object (returned from getLanguage)
#' @param prefix Defined key-prefix for obtaining word equivalents
#' @param suffix Defined key for end of word when value>1 (plural) 
#' @return String
#' @export
#' 
#
number2word = function(x, L,prefix = "numword_",suffix = "numwordPlural" ) { 

  val = floor(x) #value to use (rounded down)
  indUse = grep(prefix,names(L))
  vals = as.integer(gsub(prefix,"",names(L)[indUse]))
  if( val < min(vals) ) return( "") #don't return number if less than registered (allowed)
  itemuse = max(which(val>=vals)) #row to use
  word = L[[indUse[itemuse]]] 
  
  valUse = vals[itemuse] #value to use (basis)
  
  valRem = val-valUse #calc remainder:
  prefixVal = 1 #default
  if(valRem>0) {
    prefixVal = 10^valRem
    ending = L[[suffix]] #obtain plural ending
    word = paste0(word,ending) #add ending 
  }
  word = paste(prefixVal,word) #add number in front of word 
  return( word )
  
}