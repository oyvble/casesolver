#' @title orderTableSort
#' @description helpfunction to obtain order of a table
#' @param var1 Variable to use for sorting (sort=2)
#' @param var2 Additional variable to sort wrt (sort=3)
#' @param sort Type of sorting
#' @export

orderTableSort = function(var1,var2=NULL,sort=1) {
  ord = 1:length(var1)  #default is no sorting
  if(sort==2) {
    if(is.null(var2)) {
      ord = order(var1,decreasing=FALSE)  #sort by name (increasing name)
    } else{
      ord = order(var1,var2,decreasing=FALSE)  #sort by name (increasing name)
    }
  }
  if(sort==3 && !is.null(var2)) {
    ord = order(var2,var1,decreasing=FALSE)
  }
  return(ord)
}
