#' @title number2word
#' @description Convert number to verbal
#' @param x LR number (on log10 scale)
#' @param L Language object (returned from getLanguage)
#' @param prefix Defined key-prefix for obtaining word equivalents
#' @param suffix Defined key for end of word when value>1 (plural) 
#' @return A vector with two elements (verbal and numeric expression) 
#' @export
#' 
#
number2word = function(x, L, prefix = "numword_",suffix = "numwordPlural" ) { 

  val = floor(x) #log10 value to use (rounded down)
  indUse = grep(prefix,names(L))
  indLast = indUse[grep("more",names(L)[indUse])]
  atleast = L[[indLast]] #index of "at least"
  indUse = setdiff(indUse,indLast) #get only numbers
  vals = as.integer(gsub(prefix,"",names(L)[indUse]))
  if( val < min(vals) ) return( "") #don't return number if less than registered (allowed)
  itemuse = max(which(val>=vals)) #row to use
  word = L[[indUse[itemuse]]] 
  
  valUse = vals[itemuse] #value to use (basis)
  valRem = val-valUse #calc remainder:
  
  prefixVal = 1 #default value
  outval = val #value to report
  if(valRem>0) { #if number in front is greater than 1
    #Separate situation of being in last verbal scale or not
    if(itemuse==length(vals)) {
      prefixVal = paste(atleast,prefixVal) #include "at least" in front
      outval = valUse #use restricted value 
    } else {
      prefixVal = 10^valRem
      ending = L[[suffix]] #obtain plural ending
      word = paste0(word,ending) #add ending 
    }
  }
  word = paste(prefixVal,word) #add number in front of word 
  outval = paste0("1e",outval) #use scientific notation 
  return( c(word,outval) )
  
}