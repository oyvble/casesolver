#' @title getLanguage
#' @description A function to extract list of words for selected language (files found in installation folder)
#' @param langSel Selected language (NULL returns list of available languages)
#' @param encoding Encoding used to read text file (assuming english only)
#' @param langFile Name of language files (without extensions)
#' @param defaultLanguage Default language
#' @export

getLanguage = function(langSel=NULL,encoding="unknown",langFile="Language",defaultLanguage="English") { 
  bool = require(readxl)
  if(!bool) {
    install.packages("readxl") #request installing package
    if(!require(readxl)) {
      print("Could not install readxl..")
      return()
    }
  }
  
  pgkPath <- path.package("casesolver", quiet = FALSE) # Get package path.
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  
  languagefile1 = paste(pgkPath,paste0(langFile,".xlsx"),sep=.sep) #Contains multiple langues
  
  langTab = NULL
  tryCatch( { langTab = readxl::read_xlsx(languagefile1,sheet=1)  }, #try read from language file (in installation folder)
            error=function(e) { #if not possible to read from file:
              cat(paste0("Could not open excel language file at: ",pgkPath,"\nPlease make sure this is readable. ")) 
  })
  if(is.null(langTab)) stop("Could not open language file. Program stops!") 
  
  if(is.null(langSel)) return(colnames(langTab)[-1]) #return only language names if language not selected
  nwords = length(langTab$Key) #get total number of words
  if( nwords!=sum(!is.na(langTab[[langSel]])) ) {
    langSel = defaultLanguage #set default language
    print(paste0("Selected language did not contain all required words. Language was set back to default (",defaultLanguage,")"))
  }
  wordList = list() #insert phrases into word list
  for(i in 1:nrow(langTab)) wordList[[langTab$Key[i]]] = langTab[[langSel]][i]
  return(wordList)  #get language table
}
