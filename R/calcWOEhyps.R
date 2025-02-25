#' @title calcWOEhyps
#' @description A function for designing WoE hypothesis and calculate these, also statements are generated
#' @details The WOE window must be closed in order to return from function (wait flag)
#' @param nnTK an environment object from stored CaseSolver object
#' @param verbose Whether progress should be printed
#' @export


calcWOEhyps = function(nnTK,verbose=TRUE) {
  library(gWidgets2)
  setupAdvanced = get("setupAdvanced",envir=nnTK) #obtain match candidates to base calculation on
  useEFMex = FALSE #wheter to use EFMex calculations
  if(!is.null(setupAdvanced$useEFMex)) useEFMex = setupAdvanced$useEFMex=="TRUE"
  matchTable = get("resMatches",envir= nnTK) #obtain match candidates to base calculation on
  tableSortType = get("setupSorting",envir= nnTK)[5]  #obtain sorting order for matches table
  ord = casesolver::orderTableSort(matchTable[,1],matchTable[,2],tableSortType) #obtain selected sorted order for evids
  matchTable <- matchTable[ord,,drop=FALSE] #rearrange table wrt selected sorting

  #Store a copy of (already) existing WOE calculations
  resWOEevalBackup <- get("resWOEeval",envir=nnTK) #obtain mle fitted results
  
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language

  #Explaining CONS
  CONStxt = "CONS means CaseSolver will calculate a 'conservative LR'. \n
             The x% quantile can be selected by the user under 'Advanced->MCMC options' (5% as default). \n
             This analysis is equivalent with the 'LR sensitivity analysis' carried out in EuroForMix."

  mixLIST <- get("mixDataLIST",envir=nnTK) #get mix list
  #mixDataTABLE <- get("mixDataTABLE",envir=nnTK) #get data table from nnTK-environment
  refDataTABLE <- get("refDataTABLE",envir=nnTK) #get data table from nnTK-environment
  if(is.null(mixLIST) || is.null(refDataTABLE)) {
    print("No data to perform analysis on.")
    return()#no evid/ref data to compare
  }

  #extract profile names (used in tables):
  evidNames = names(mixLIST)
  refNames = rownames(refDataTABLE)
  
  #HELPFUNCTIONS:
  
  #helpfunction when extracting reference data
  getRefL = function(refs) { #return list with same order as for refs
    refL <- casesolver::tabToListRef(tab=refDataTABLE[ match(refs,rownames(refDataTABLE)),,drop=FALSE],setEmpty=TRUE) #FORCING 1-alleles to be empty
    return(refL)
  }
  
  
  #Helpfunction to get hypothesis text
  getHypTxt = function(mlefit,getNamesOnly=FALSE) {
    model = mlefit$model
    condOrder = model$condOrder
    nU = model$nC - sum(condOrder>0)
    locNames =  names(mlefit$model$popFreq) #locus names
    refNames = names(model$refData) #obtain reference names
    
    #If reference names was any of the locus names we have another format:
    if( any(toupper(refNames)%in%toupper(locNames)) ) refNames = names(model$refData[[1]])
    
    condNames <- refNames[which(condOrder>0)]
    hyptxt = paste0(condNames,collapse="/")
    if(nU>0) {
      hyptxt = paste(hyptxt,L$and,nU)
      if(nU==1)  hyptxt = paste(hyptxt,L$unknown)
      if(nU>1)  hyptxt = paste(hyptxt,L$unknowns)
      hyptxt = paste0(hyptxt,"*")
      condNames = c(condNames,paste0("U",1:nU)) #update with unknowns
    }
    if(getNamesOnly) return(condNames) #return only vector of names
    return(hyptxt)
  }
  
  #Evaluating hypotheses based on table information
  f_eval = function(h,...) {
    if(verbose) print("Performing calculations...")
    
    #TRAVERSE ALL HYPOTHESIS SETS and store result in nnTK environment
    nHyps = length(evalList)
    if(nHyps==0) {
      if(verbose) print("No hypotheses to evaluate!")
      return()
    }
    
    evalList2 = list() #new list to insert results

    #IMPORTANT: APPEND PREVIOUS CALCULATIONS (resWOEevalBackup) WITH NEW CALCULATIONS:
    if(!is.null(resWOEevalBackup)) { #if previous calculations already conducted
      evalList2 = resWOEevalBackup #copy already existing
      evalList2$resTable = NULL #reset table as it will be updated later
    }   
    
    for(i in 1:nHyps) {
#      i=2
      hyptxt = paste0("Hypothesis #",i)
      if( !gWidgets2::svalue(grid[i+1,1]) ) next #skip if not selected

      NOC = as.integer(gWidgets2::svalue(grid[i+1,5])) #obtain number of contributors
      calcCons = gWidgets2::svalue(grid[i+1,6])  #whether to calculate conservative LR
      evids = evalList[[i]]$evid
      poi = evalList[[i]]$poi
      conds = evalList[[i]]$cond
      
      nCond = length(conds) #number of conditional
      nUhd = NOC - nCond #number of unknowns under Hd
      nUhp = nUhd - 1 #number of unknowns under Hp (assuming fixed possibility)
      
      #Skip hypothesis sets which can't be evaluated
      if(length(evids)==0 || length(poi)==0) {
        if(verbose) print("No evidence or poi to evaluate. Ignoring hypothesis set...")
        next
      } else if( nCond >= NOC ) { 
        if(verbose) print("Number of conditionals was equal or greater than number of contributors. Ignoring hypothesis set...")
        next 
      }  
      
      #CONSTRUCT  TXT BASED
      evidtxt = paste0(evids,collapse="/")
      condtxt = paste0( conds ,collapse="/")

      if(verbose) {
        print(paste0("---------",hyptxt,"--------------"))
        print(paste0("Evidence(s): ",evidtxt ))
        print(paste0("Person of Interest(s) (POI): ", paste0(poi,collapse="/") ))
        print(paste0("Conditional references: ",condtxt))
        print(paste0("Number of contributors: ",NOC))
      }

      #Obtain and prepare data:
      nPOI = length(poi)
      refs = c(poi,conds) #Note: poi always first index
      evidData = mixLIST[evids]  #obtain evidence data in list
      refData = getRefL(refs) #obtain reference data in list
      #all(names(refData)==refs)

      if(useEFMex) { #perform calculation based on EFMex
        mod <- casesolver::getModelSettings(nnTK) #Obtain model settings
        condRefNames = NULL
        if(nCond>0) condRefNames = conds
        calc = EFMex::calculateExhaustive(evidData,refData, mod$popFreq, NOC, mod$kit,  condRefNames, #data settings
                                          modelSetting = list(AT=mod$threshT,fst=mod$fst, pC=mod$pC, lambda=mod$lambda), #model settings
                                          modelConfig = list(degrade=!is.null(mod$kit), stutterBW=is.null(mod$xiBW),stutterFW=is.null(mod$xiFW),priorBWS = NULL,priorFWS = NULL), #model config
                                          optimConfig = list(seed=NULL,nDone=mod$nDone, steptol=1e-4, minF=mod$minFreq, normalize=mod$normalize,adjQbp=FALSE),
                                          storeModelFit = TRUE, verbose = TRUE, LRthresh=1)
      }
      
      #TRAVERSING THROUGH EACH POI (for both approaches)    
      for(j in seq_len(nPOI)) { 
#          j=1
        poiSel = poi[j] #select specific POI (also when not using EFMex)
             
        if(useEFMex) { #special handling results from EFMex
          jointlr = calc$LRtable[poiSel,"GLR"] #obtain LR
          POIinHyp = calc$hypCalcs[,poiSel]>0 #obtain whether POI is considered
          POIinHypIdx = which(POIinHyp)
          POIninHypIdx = which(!POIinHyp)
          logliks = calc$hypCalcs[,nPOI+1] #this is logLik columns
          bestFitHpIdx = POIinHypIdx[which.max(logliks[POIinHyp])]
          bestFitHdIdx = POIninHypIdx[which.max(logliks[!POIinHyp])]
          fithp = calc$mleFit[[bestFitHpIdx]]
          fithd = calc$mleFit[[bestFitHdIdx]]
          POIcontrIdx = calc$hypCalcs[bestFitHpIdx,poiSel]
          MxPOI = fithp$fit$thetahat2[POIcontrIdx]
          
        } else { 
          #Evaluate LR in ordinary way:
            
          condhp <- 1:length(refs) #conditional order must be increasing order
          condhd <- condhp - 1 #don't condition on POI under Hd
          
          #RUN CALCULATIONS:
          fitMLE = function(condOrder,knownRef=NULL) {
            mlefit <- casesolver::calcQuanMLE(evidData,refData,condOrder,NOC,nnTK,isWOE=TRUE,knownRef=knownRef)
            return(mlefit)
          }
          if(verbose) print("calculating MLE under Hp...")
          fithp <- fitMLE(condhp)
          if(verbose) print(paste0("logLik (Hp): ",fithp$fit$loglik))
          
          if(verbose) print("calculating MLE under Hd...")
          fithd <- fitMLE(condhd,knownRef=1) #put index of known non-contributor Ref as 1
          if(verbose) print(paste0("logLik (Hd): ",fithd$fit$loglik))
          
          #Extract certain parameters:
          MxPOI = fithp$fit$thetahat2[1]
        } #end if not using EFMex
    
        #Obtain hypothesis references for obtained MLE fit objects
        hptxt = getHypTxt(fithp) #only text
        hdtxt = getHypTxt(fithd) #only text
        condhpNames = getHypTxt(fithp,getNamesOnly=TRUE) #also as vector
        condhdNames = getHypTxt(fithd,getNamesOnly=TRUE) #also as vector

        #CONTINUE (SAME FOR BOTH APPROACHS)
        #Evaluate LR per marker and validation:
        logHp <- euroformix::logLiki(fithp) #log Pr(Data|Hp)
        logHd <- euroformix::logLiki(fithd) #log Pr(Data|Hd)
        
        mleLR <- (fithp$fit$loglik - fithd$fit$loglik)/log(10) #get estimated log10LR (Quan based)
        mleLRi <- exp(logHp-logHd) #LR per marker (ordinary scale) /log(10)
        
        kit = casesolver::getEnvirKit(nnTK) #obtain kit
        if(verbose) print("calculating valid under Hp...")
        validHp <- euroformix::validMLEmodel(fithp,kit=kit, createplot=FALSE,alpha=0.01,verbose=FALSE)
        if(verbose) print("calculating valid under Hd...")
        validHd <- euroformix::validMLEmodel(fithd,kit=kit, createplot=FALSE,alpha=0.01,verbose=FALSE)
        nFailedHp=sum(validHp$Significant)
        nFailedHd=sum(validHd$Significant)
        if(verbose) print(paste0("log10LR (mle): ",round(mleLR,2)))
        
        #store table with fitted parameters (Hp vs Hd). Also include adjusted logLik
        paramHp = fithp$fit$thetahat2
        paramHd = fithd$fit$thetahat2
        tab = cbind(Hp=paramHp,Hd=paramHd)
        paramTable = round(tab,3) #round table
        paramTable[NOC+1,] = round(paramTable[NOC+1,]) #don't round PHexp
        
        #Also add adjusted logLIk:
        getAdjLogLik = function(mle) mle$fit$loglik - length(mle$fit$thetahat)
        adjLogLik = round( c(getAdjLogLik(fithp),getAdjLogLik(fithd)),2)
        paramTable = rbind(paramTable,adjLogLik) #add adj.loglik
        
        #obtain names with mixture proportions for conditional references:
        MxRefs = list(hp=setNames(paramHp[1:NOC],condhpNames),hd=setNames(paramHd[1:NOC],condhdNames))
  
        #If calculating conservative LR
        consLR <- bayesLR <- NA #none by default
        mcmc = NULL #default mcmc object (NULL if not run)
        if(calcCons) {
          mcmcOpt = get("setupMCMC",envir=nnTK)
          if(verbose) print("Performing MCMC to estimate conservative LR")
          
          #verbose=TRUE
          mcmc = euroformix::calcLRmcmc(fithp,fithd, mcmcOpt$niter,mcmcOpt$delta,mcmcOpt$quantile,mcmcOpt$seed, verbose=TRUE,traceplot=FALSE)
          #cons = calcCONS(fithp,fithd, mcmcOpt$niter,mcmcOpt$delta,mcmcOpt$quantile,mcmcOpt$seed, verbose=TRUE)
          consLR = mcmc$log10LRcons #extract conservative LR
          bayesLR = mcmc$log10LRbayes # extract  bayesian based LR
          if(verbose) {
            print(paste0("log10LR (CONS): ",round(consLR,2)))
            print(paste0("log10LR (BAYES): ",round(bayesLR,2)))
          }
        } #end if conservative
        
        #Construct statement
        state = L$statementScheme  #obtain statement (includes $ to indicate text)
        state = gsub("$evidtxt",evidtxt,state,fixed=TRUE)
        state = gsub("$hptxt",hptxt,state,fixed=TRUE)
        state = gsub("$hdtxt",hdtxt,state,fixed=TRUE)
        state2 <- state #store a copy before inserting LRvals (makes it possible to insert it afterwards)
  
        #insert woe strength      
        LRuse = mleLR
        if(calcCons) LRuse = consLR
        LRtxt = signif(10^LRuse,2)
        
        #Obtain value as text (if )
        useVerbalLR =  get("setupReportOpt",envir=nnTK)["verbalLR"]  #If using verbal LR
        if(useVerbalLR) {
          LRverbal = number2word(LRuse,L) #obtain verbal LR
          LRtxt = paste0(LRverbal[1]," (",LRverbal[2],")") #add verbal in additon to down-rounded number
        }
        state = gsub("$LRtxt",LRtxt,state,fixed=TRUE)
        # if(length(LRtxt)>1) print(LRtxt)
        
        #statement = paste0("The evidence (",evidtxt,") is ",LRtxt," times more likely if the DNA came from ",hptxt," than if it came from ",hdtxt)
        evalList2[[length(evalList2)+1]] = list(evid=evids,poi=poiSel,cond=conds,NOC=NOC,
                                     mleHp=fithp,mleHd=fithd, MxPOI=MxPOI,
                                     nFailedHp=nFailedHp,nFailedHd=nFailedHd,
                                     mleLR=mleLR,mleLRi=mleLRi,consLR=consLR,bayesLR=bayesLR,
                                     statement=state,statementNoWOE=state2,
                                     paramTable=paramTable,MxRefs=MxRefs,
                                     mcmc = mcmc) #also store MCMC run (this is new)
      } #end for each POI
    } #end for each hypothesis
    assign("resWOEeval",evalList2,envir=nnTK)  #store match-results from comparison (those with LR>threshold)
    gWidgets2::dispose(win) #close main window
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
  } #end evaluation function
  
  #Edit profiles in hyp-set (dynamic update of GUI table)
  f_profileEdit = function(h) {
    objStore = "itemList"
    type = h$action$typ #obtain type
    hypID  = as.integer(h$action$hyp) #Obtain hypothesis ID
    
    items1 = evalList[[hypID]][[type]] #obtained stored elements for specific type
    #items1 = gsub(stringMultiple,"",items1) #remove ++ from name
    if(type=="evid") {
      items2 = setdiff(evidNames,items1)
    } else {
      items2 = setdiff(refNames,items1)
    }
    #create new environment for using profileSwapperGUI:
    env = new.env( parent = emptyenv() ) 
    assign(objStore,list(items1=items1,items2=items2),envir=env)
    profileSwapperGUI(env) #USER MAY CHANGE ITEMS
    
    #OBTAIN SELECTED ITEM(S) FROM USER:
    items1 = get(objStore,envir=env)$items1 
    #print(item1)
    tooltipItems = items1
    if(length(items1)==0) tooltipItems = "" #set as empty if none selected

    if(type=="poi") {
      if(length(items1)==0) {
        gWidgets2::gmessage("A profile must be selected.")
      } else if(!useEFMex && length(items1)>1) {
        gWidgets2::gmessage("Exactly one profile must be selected.")
      }
    }     
    
    evalList[[hypID]][[type]] <<- items1 #insert selected items (must be before modifying it)
    
    if( length(items1)>1 ) { #if multiple elements selected:
      items1 = paste0(items1[1],stringMultiple) #Select only first as reference #"multiple"
    }  else if(length(items1)==0) { #if no elements selected
      items1 = L$none 
    }
    
    gWidgets2::tooltip(grid[h$action$ij[1],h$action$ij[2]]) <- tooltipItems #always insert as tooltip
    gWidgets2::svalue(grid[h$action$ij[1],h$action$ij[2]]) <- items1 #insert button value
  }
  
  #helpfunction for selecting all (or none) hyp IDs
  f_selAllHyp = function(h) {
    if(IDcount>0) {
      for(id in 1:IDcount) gWidgets2::svalue(grid[id+1,1]) <- h$action #traverse all hyps and select all or none
    }
  }
  
  #helpfunction for selecting all (or none) CONS
  f_selAllCONS = function(h) {
    if(IDcount>0) {
      for(id in 1:IDcount) gWidgets2::svalue(grid[id+1,6]) <- h$action #traverse all hyps and select all or none
    }
  }
  
  #helpfunction for unconditioning all
  f_selUNCOND = function(h) {
    if(IDcount>0) {
      for(id in 1:IDcount) { #traverse all hyps 
        gWidgets2::svalue(grid[id+1,4]) <- stringNone #set all as "none"
        evalList[[id]]$cond <<- emptyObj #set as not conditional 
      } 
    }
  }
  
  #Setup:
  IDcount = 0 #counter for hypothesis set
  stringNone = L$none #"none"
  stringMultiple = "++" #"multiple"
  NOCrange = 1:5 #number can be edited in textbox
  emptyObj =   as.character()
  evalList = list() #Construct object (table) with hypothesis 
  
  #Obtain window from table:
  flag <- tcltk::tclVar("") #init flag to avoid quitting GUI
  win <- gWidgets2::gwindow( "Hypothesis window for weight of evidence" ,visible=FALSE, width=1000, height=500)
  gWidgets2::addHandlerUnrealize( win, handler = function(h,...) { #
    tcltk::tclvalue(flag) <- "destroy" #Destroy wait flag
    return(NULL) #return selected profiles when window is exit 
  }  ) #call quit function
  
  frame <- gWidgets2::ggroup(horizontal = FALSE, container=win,use.scrollwindow = T, fill=T,expand=T)
  
  #insert buttons:
  eval = gWidgets2::glayout(spacing=3,container=frame) 
  eval[1,1] = gWidgets2::gbutton("Add hyp. set",container=eval, handler=function(h) { insGridRow(NULL) })#Add another hypothesis set 
  eval[1,2] = gWidgets2::gbutton("Select all",container=eval, handler=f_selAllHyp, action=TRUE)
  eval[1,3] = gWidgets2::gbutton("Unselect all",container=eval, handler=f_selAllHyp, action=FALSE)
  eval[1,4] = gWidgets2::gbutton("Evaluate",container=eval, handler=f_eval)
  
  #Add choices for Conservative
  eval[2,1] = gWidgets2::glabel("CONS:",container=eval)
  eval[2,2] = gWidgets2::gbutton("Select all",container=eval, handler=f_selAllCONS, action=TRUE)
  eval[2,3] = gWidgets2::gbutton("Unselect all",container=eval, handler=f_selAllCONS, action=FALSE)
  
  #Possible to "uncondition for all hypotheses" (NEW IN v1.9)
  eval[2,4] = gWidgets2::gbutton("Uncondition all",container=eval, handler=f_selUNCOND, action=FALSE)
  
  gWidgets2::gseparator(frame,horizontal = T)
  grid = gWidgets2::glayout(spacing=3,container=frame) 
  grid[1,1] = gWidgets2::glabel("Set",container=grid)
  grid[1,2] = gWidgets2::glabel("Evid(s)",container=grid)
  grid[1,3] = gWidgets2::glabel("POI",container=grid)
  grid[1,4] = gWidgets2::glabel("Cond(s)",container=grid)
  grid[1,5] = gWidgets2::glabel("NOC",container=grid)
  grid[1,6] = gWidgets2::glabel("CONS",container=grid)
  
  #Adding tooltip about CONS label:
  gWidgets2::tooltip(eval[2,1]) <- gWidgets2::tooltip(grid[1,6]) <- CONStxt
  
  #Function to insert hypothesis sets into GUI
  insGridRow = function(i) {
    #val=opt$priority[i] #obtain value
    if(is.null(i)) {
      IDcount <<- IDcount + 1  #update counter
      grid[IDcount+1,1] <<- gWidgets2::gcheckbox( paste0("#",IDcount) ,checked=TRUE,container=grid)
      grid[IDcount+1,2] <<- gWidgets2::gbutton(stringNone, container=grid, handler=f_profileEdit, action=list(typ="evid", hyp=IDcount,ij=c(IDcount+1,2) ))
      grid[IDcount+1,3] <<- gWidgets2::gbutton(stringNone, container=grid, handler=f_profileEdit, action=list(typ="poi", hyp=IDcount,ij=c(IDcount+1,3)))
      grid[IDcount+1,4] <<- gWidgets2::gbutton(stringNone, container=grid, handler=f_profileEdit, action=list(typ="cond", hyp=IDcount,ij=c(IDcount+1,4)))
      grid[IDcount+1,5] <<- gWidgets2::gcombobox(items=NOCrange,selected=2,editable=TRUE,container=grid)
      gWidgets2::size( grid[IDcount+1,5] ) = 2
      grid[IDcount+1,6] <<- gWidgets2::gcheckbox( "",checked=FALSE,container=grid)
      evalList[[IDcount]] <<-  list(evid=emptyObj,poi=emptyObj,cond=emptyObj) #add empty object
    }  else {
      refs = unlist(strsplit(matchTable[i,2],"/"))
      nRefs = length(refs)
      if(nRefs==0) return(NULL) #return if none
      evid = matchTable[i,1]
      NOC =  matchTable[i,3]
      
      #BEFORE INSERTING HYPOTHESIS: CHECK IF ALREADY CALCULATED (then not include):
      if(!is.null(resWOEevalBackup) && nrow(resWOEevalBackup$resTable)>0) {
        bool1 = resWOEevalBackup$resTable[,1]%in%evid
        bool2 = resWOEevalBackup$resTable[,2]%in%refs
        #bool3 = resWOEevalBackup$resTable[,3]%in%condcheck #removed from v1.8.1
        #bool4 = resWOEevalBackup$resTable[,4]%in%NOC
        if( any(bool1 & bool2) ) return(NULL) #Skip suggested WOE if already calculated
      }
      
      
      if(useEFMex) { #run EFMex if 
        poi = refs #allow multiple POIs

        IDcount <<- IDcount + 1  #update counter (pre)
        grid[IDcount+1,1] <<- gWidgets2::gcheckbox( paste0("#",IDcount) ,checked=TRUE,container=grid)
        
        #EVID BUTTON
        grid[IDcount+1,2] <<- gWidgets2::gbutton( evid , container=grid, handler=f_profileEdit, action=list(typ="evid", hyp=IDcount, ij=c(IDcount+1,2) ))
        
        #POI BUTTON
        if(nRefs>1)  poi = paste0(poi[1],stringMultiple) #stringMultiple
        grid[IDcount+1,3] <<- gWidgets2::gbutton( poi , container=grid, handler=f_profileEdit, action=list(typ="poi", hyp=IDcount, ij=c(IDcount+1,3)))
        gWidgets2::tooltip( grid[IDcount+1,3]) = paste0(refs,collapse="\n")
        
        #No conditional by default        
        grid[IDcount+1,4] <<- gWidgets2::gbutton( stringNone, container=grid, handler=f_profileEdit, action=list(typ="cond", hyp=IDcount, ij=c(IDcount+1,4)))

        #Define NOC        
        grid[IDcount+1,5] <<- gWidgets2::gcombobox(items=NOCrange,selected=NOC,editable=TRUE,container=grid)
        gWidgets2::size( grid[IDcount+1,5] ) = 2
        
        #Last column is whether to consider CONSERVATIVE CALCULATIONS
        grid[IDcount+1,6] <<- gWidgets2::gcheckbox( "",checked=FALSE,container=grid)
        #gWidgets2::enabled(grid[IDcount+1,6]) <<- FALSE #deactivating box
        
        #Insert to list:
        evalList[[IDcount]] <<- list(evid=matchTable[i,1],poi=refs,cond=emptyObj)
        
      } else {
        for(j in 1:nRefs) {
          poi = refs[j] #obtain j'th ref (set as POI)
          conds = setdiff(refs,refs[j]) #obtain other refs
          nConds = length(conds) #number of conditionals
          
          IDcount <<- IDcount + 1  #update counter (pre)
          grid[IDcount+1,1] <<- gWidgets2::gcheckbox( paste0("#",IDcount) ,checked=TRUE,container=grid)
          
          #EVID BUTTON
          grid[IDcount+1,2] <<- gWidgets2::gbutton( evid , container=grid, handler=f_profileEdit, action=list(typ="evid", hyp=IDcount, ij=c(IDcount+1,2) ))
          
          #POI BUTTON
          grid[IDcount+1,3] <<- gWidgets2::gbutton( poi , container=grid, handler=f_profileEdit, action=list(typ="poi", hyp=IDcount, ij=c(IDcount+1,3)))
          
          if(nConds==0) {
            cond = stringNone
          } else if(nConds==1) { #if only one
            cond = conds
          } else {
            cond = paste0(conds[1],stringMultiple) #stringMultiple
          }
          
          #COND BUTTON
          grid[IDcount+1,4] <<- gWidgets2::gbutton( cond, container=grid, handler=f_profileEdit, action=list(typ="cond", hyp=IDcount, ij=c(IDcount+1,4)))
          if( nConds>0 ) {
            gWidgets2::tooltip( grid[IDcount+1,4]) = paste0(conds,collapse="\n")
          }
          
          grid[IDcount+1,5] <<- gWidgets2::gcombobox(items=NOCrange,selected=NOC,editable=TRUE,container=grid)
          gWidgets2::size( grid[IDcount+1,5] ) = 2
          
          #Last column is whether to consider CONSERVATIVE CALCULATIONS
          grid[IDcount+1,6] <<- gWidgets2::gcheckbox( "",checked=FALSE,container=grid)
          
          #Insert to list:
          evalList[[IDcount]] <<- list(evid=matchTable[i,1],poi=refs[j],cond=conds)
        } #end for each ref
      }
    }
  }
  
  #INSERT GRID TABLE  
  if(!is.null(matchTable) && nrow(matchTable)>0) {
    for(i in 1:nrow(matchTable)) {
      insGridRow(i)
    }
  }
  gWidgets2::visible(win) = TRUE
  gWidgets2::focus(win) = TRUE
  
  tcltk::tkwait.variable(flag) #important to not quit window before broken
  
}

