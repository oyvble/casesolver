#' @title getFreqs
#' @description get allelefreq-info from file
#' @param filename A string giving the file with allelefrequencies
#' @export

 getFreqs = function(filename) {
    tab=euroformix::tableReader(filename)
   Anames = tab[,1] #first column is allele frequeneies
   tab = tab[,-1] 
   freqlist = list()
   for(j in 1:ncol(tab)) { #for each locus
     tmp = tab[,j]
     tmp2 = tmp[!is.na(tmp)]
     names(tmp2) = Anames[!is.na(tmp)]
     freqlist[[j]] = tmp2
   }
   names(freqlist) = toupper(colnames(tab)) #LOCUS-names are assigned as Upper-case! This is important to do!
   return(freqlist)
 }
 