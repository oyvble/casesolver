#' @title sampleListToTable
#' @description  A function for converting profiles from list to table.
#' @description Converting profiles from list to table format (helpfunction)
#' @param outL A list on form outL[[samplename]][[locusname]]$adata,outL[[samplename]][[locusname]]$hdata (Evid or refs)
#' @param PH Boolean of whether Peak heights should be exported (used for LRmix studio)
#' @return table 
#' @export

 sampleListToTable = function(outL,PH=TRUE) {
   #Includes peak heights only if PH=TRUE
   sn = names(outL) #outL is a list on form outL[[samplename]][[locusname]]$adata,outL[[samplename]][[locusname]]$hdata
   aM = 0   #count number of max allele data:
   hM = 0   #count number of max allele heights:
   for(ss in sn) { #for each sample
    aM = max(unlist( lapply(outL[[ss]],function(x) length(x$adata)) ),aM)
    hM = max(unlist( lapply(outL[[ss]],function(x) length(x$hdata)) ),hM)
   }
   #create tables:
   X=numeric()
   for(ss in sn) { #for each sample
    newsample=numeric() #for allele
    ln = names(outL[[ss]])
    for(loc in ln) {
     newrow = outL[[ss]][[loc]]$adata
     newsample = rbind(newsample, c(newrow,rep("",aM-length(newrow))))
    }
    newsample2=numeric() #for heights
    if(hM>0 && PH) {
     for(loc in ln) {
      newrow = outL[[ss]][[loc]]$hdata
      newsample2 = rbind(newsample2, c(newrow,rep("",hM-length(newrow))))
     }      
    }
    X = rbind(X,cbind(ss,ln,newsample,newsample2))
   }
   cn = c("SampleName","Marker", paste("Allele",1:aM,sep=""))
   if(hM>0 && PH) cn = c(cn,paste("Height",1:hM,sep=""))
   colnames(X)  = cn
   return(X)
 } #end of functions
