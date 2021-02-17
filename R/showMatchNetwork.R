#' @title showMatchNetwork
#' @description Calculate identical by state statistics between all references
#' @param nnTK an environment object from stored CaseSolver object (includes data)
#' @param action The type of plot to be produced (all,onlymix,onlyss)
#' @param createInteractive A boolean of whether to create an interactive matchnetwork (shown in browser)
#' @param selList A list (evids,refs) to only show the matchnetwork for. Default is not considered 
#' @return A list with canidate results
#' @export

showMatchNetwork = function(nnTK,action,createInteractive=FALSE,selList=NULL) {
  L = casesolver::getLanguage( get("setupLanguage",envir=nnTK)$language ) #, get("setupLanguage",envir=nnTK)$encoding ) #get list of words from selected language
  
  matchList = get("allMixList",envir=nnTK) #get list with all matches (Matches panel)
  if(nrow(matchList)==0) return() #return if no matches
  
  defaultScore=20 #default score for log10 value (assuming no missmatch)
  tmp = strsplit(matchList[,2],"/") #get all comparison matches
  nMatches = sapply(tmp,length) #get number of matches per evid profiles
  nContrMatch = rep(matchList[,3],nMatches) #number of contributors
  EvidMatch = rep(matchList[,1],nMatches) #evidence profiles
  RefMatch = unlist(tmp) #reference profiles
  ShowMatch = cbind(EvidMatch,RefMatch,nContrMatch,defaultScore)
  
  #Filter matches based on thresholds (after doing compare)
  tab <- get("resCompLR",envir=nnTK) #get last stored LR comparison values 
  scoreColumnInd = 4 #default column with score (LR)
  if( is.null(tab) ) { #no LR-results found: Checking if MAC was calculated..
    tab = get("resCompMAC",envir=nnTK)$MatchList #extract MAC based results
    if( !is.null(tab) ) scoreColumnInd = 3 #column with score (MAC)
  }
  
  #Update scores in "ShowMatch" based on scores from "compared match results"  (MAC, qualLR,quanLR) :  
  if( !is.null(tab) ) { #if any results to filter on (MAC, qualLR,quanLR) 
    key1 = paste0(ShowMatch[,1],ShowMatch[,2])
    key2 = paste0(tab[,1],tab[,2])
    updateScores  =  key2%in%key1 #index to update score for
    if( any(updateScores) ) { #should any scores be update
      tab = tab[updateScores,,drop=FALSE] #update tab with scores
      #all(key1[ match(key2[updateScores],key1)  ] == key2[updateScores]) #all must be true
      ShowMatch[ match(key2[updateScores],key1) ,4] = tab[,scoreColumnInd] #update score
    }
  }
  
  if(createInteractive) { #IF BUTTON CLICKED: 
    if(action=="onlymix") ShowMatch = ShowMatch[ShowMatch[,3]!="1",,drop=FALSE]  # Check if comparing only against Mixture
    if(action=="onlyss") ShowMatch = ShowMatch[ShowMatch[,3]=="1",,drop=FALSE]  # Check if comparing only against Mixture
  }
  
  #FILTER ON DATA SELECTION:
  if(!is.null(selList)) {
    keep  =  ShowMatch[,1]%in%selList$evids & ShowMatch[,2]%in%selList$refs #NOTE THE OREDER
    ShowMatch = ShowMatch[keep,,drop=FALSE]
  }
  
  if(nrow(ShowMatch)==0) {
    if(createInteractive) gWidgets2::gmessage("No candidate matches to show!") #if no candidates found
    return(FALSE)
  }
    
  getWeight = function(x) {
    x[x<1] = 1 #set 1 as minimum score
    return(sqrt(x))
  }
  #rem <- duplicated(tab[,1:2]) #indices to remove
  #if(nrow(tab)==1) rem <- FALSE #duplicated doesnt work for case having only 1 row
  tab2 <-  data.frame(from=ShowMatch[,2],to=ShowMatch[,1],weight=getWeight(as.numeric(ShowMatch[,4])),nC=as.integer(ShowMatch[,3]))
  gg <- igraph::graph.data.frame(tab2,directed=FALSE)
  nods <- names(igraph::V(gg))
  
  #Colorcode circles wrt number of contributors
  cols <- rep("green",length( nods  )) #all reference profiles are green
  #cols[grepl("Unknown ",nods)] <- "cyan" #unknowns 
  cols[nods%in%tab2$to[tab2$nC==1]] <- "cyan" #assigned as single sources (becomes blue)
  cols[nods%in%tab2$to[tab2$nC==2]] <- "orange" 
  cols[nods%in%tab2$to[tab2$nC>2]] <- "red" 
  
  layout = igraph::layout_with_kk(gg,weights = igraph::E(gg)$weight) #selected layout: Alternative=layout_with_lgl
  ewidth = igraph::E(gg)$weight
  igraph::plot.igraph(gg,layout=layout,edge.width=ewidth,vertex.color=cols,vertex.size=10,vertex.label.color="black",vertex.label.cex=0.8,main=paste( L$matches , L$forr ,get("caseID",envir=nnTK)) )
  
  #Create interactive plots if h is different from NULL (this is only if "ShowMatchNetwork" button is clicked)
  if(createInteractive) {
    w0 = 1920 #width of plot
    h0 = 1200 #height of plot
    cols[cols=="green"] = "lime" 
    es <- as.data.frame(igraph::get.edgelist(gg),stringsAsFactors=FALSE)
    edge_shapes <- list()
    for(i in 1:nrow(es)) {
      v0 <- es[i,1]==nods #from ind
      v1 <- es[i,2]==nods #to ind
      edge_shapes[[i]] = list(type = "line",
                              line = list(color = "#030303", width =ewidth[i],dash="solid"),
                              opacity = 0.3,layer="below",x0 = layout[v0,1],y0 = layout[v0,2],x1 = layout[v1,1],y1 = layout[v1,2] )
    }
    conv = rep("",length(nods)) 
    for(vv in nods) { #for each vertex we indicate the links
      ind1 = es[,1]==vv #get from/to list  
      ind2 = es[,2]==vv #get from/to list  
      nam = unique(c(es[ind1,2],es[ind2,1]))
      conv[nods==vv] = paste0(vv,":\n",paste0(sort(nam),collapse="\n"))
    }
    ord = order(nods) #sort wrt this name
    df = data.frame(x=layout[ord,1],y=layout[ord,2],what=nods[ord],con=conv[ord],stringsAsFactors=FALSE)
    
    txtsz=12 #text size
    style = list(size = 40, color = cols[ord] )#, line = list(color = 'black', width = 1))
    txtfont = list(color = '#000000', size = txtsz)
    ax = list(showline=FALSE,showticklabels=FALSE,zeroline=FALSE,showgrid=FALSE,title="") #axis setup: none
    out <- plotly::plot_ly(df, x = ~x, y = ~y, mode = "markers+text", type="scatter", text = ~what,name= ~con , hoverinfo = "name",hoverlabel=list(font=list(size=txtsz),namelength=1000,bgcolor="white",bordercolor="black"), textfont=txtfont, marker=style )
    out <- plotly::layout(out,xaxis=ax,yaxis=ax,shapes = edge_shapes)
	out <- plotly::hide_legend(out)
	out <- plotly::config(out, scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0,height=h0)) 
    print( out ) 
  } #end if plotly
  return(TRUE)
}
