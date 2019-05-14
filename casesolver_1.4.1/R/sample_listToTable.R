#' @title sample_listToTable
#' @description  A function for converting profiles from list to table.
#' @export

 sample_listToTable = function(Y,PH=TRUE) {
   #Includes peak heights only if PH=TRUE
   sn = names(Y) #Y is a list on form Y[[samplename]][[locusname]]$adata,Y[[samplename]][[locusname]]$hdata
   aM = 0   #count number of max allele data:
   hM = 0   #count number of max allele heights:
   for(ss in sn) { #for each sample
    aM = max(unlist( lapply(Y[[ss]],function(x) length(x$adata)) ),aM)
    hM = max(unlist( lapply(Y[[ss]],function(x) length(x$hdata)) ),hM)
   }
   #create tables:
   X=numeric()
   for(ss in sn) { #for each sample
    newsample=numeric() #for allele
    ln = names(Y[[ss]])
    for(loc in ln) {
     newrow = Y[[ss]][[loc]]$adata
     newsample = rbind(newsample, c(newrow,rep("",aM-length(newrow))))
    }
    newsample2=numeric() #for heights
    if(hM>0 && PH) {
     for(loc in ln) {
      newrow = Y[[ss]][[loc]]$hdata
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
