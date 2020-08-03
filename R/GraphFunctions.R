# 2016-11-04
# Graph functions

PlotDAGGraph <- function(sigDAG, enrichment.sum, jfontsize = 100, jtitle = "", amp.cutoff = 1){
  # BEGIN PLOTTING
  
  allTerm.df <- GetGONames(as.df = TRUE, filter.GO.ids=nodes(sigDAG))
  goid.to.name <- GetGONames(as.df = FALSE, filter.GO.ids=FALSE)
  
  graphAttrs <- getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL
  graphAttrs$node$shape <- 'ellipse'
  graphAttrs$node$fontsize <- '20'
  
  nodeAttrs <- list()
  edgeAttrs <- list()
  
  nodeAttrs$label <- sapply(nodes(sigDAG), function(goid) AssignHash(goid, goid.to.name, null.fill = NA))
  
  weightsList <- edgeWeights(sigDAG)
  to <- lapply(weightsList, names)
  from <- nodes(sigDAG)
  edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
  edge.weights <- unlist(weightsList)
  names(edge.weights) <- edge.names
  ##    0 for a is_a relation,  1 for a part_of relation
  edgeAttrs$color <- ifelse(edge.weights == 0, 'black', 'red')
  
  # nodeAttrs$label <- seq(length(nodeAttrs$label))
  nodeAttrs$fixedsize[nodes(sigDAG)] <- TRUE
  nodeAttrs$fontsize[nodes(sigDAG)] <- jfontsize
  # add node colors based on HSV
  # allTerm.df <- data.frame(goid = allDAGTerm, name = allDAGTerm.name, stringsAsFactors = FALSE)

  pval.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$pval)
  amp.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$amp)
  phase.hash <- hash(enrichment.sum$GO.ID, enrichment.sum$phase)
  allTerm.df$pval <- sapply(allTerm.df$goid, function(g) AssignHash(g, pval.hash, null.fill = NA))
  allTerm.df$amp <- sapply(allTerm.df$goid, function(g) AssignHash(g, amp.hash, null.fill = NA))
  allTerm.df$phase <- sapply(allTerm.df$goid, function(g) AssignHash(g, phase.hash, null.fill = NA))
  allTerm.df$Color <- PhaseAmpPvalToColor(phase = allTerm.df$phase, amp = allTerm.df$amp, pval = allTerm.df$pval, 
                                          rotate.hr = -8, amp.k = amp.cutoff, pval.k = Inf, method = "cutoff", black.to.white = TRUE)
  color.hash <- hash(allTerm.df$goid, allTerm.df$Color)
  nodeAttrs$fillcolor <- sapply(names(nodeAttrs$label), function(l) AssignHash(l, color.hash, null.fill = NA))
  plot(sigDAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, main = jtitle)
}

GetGONames <- function(as.df = TRUE, filter.GO.ids = FALSE){
  conn <- get("GO_dbconn")()
  .sql <- paste("select distinct go_id goid,term name from go_term where ontology='",
                toupper(ontology), "'", sep="")
  allTermName <-  dbGetQuery(conn,.sql)
  allTermName <- as.data.frame(allTermName)
  if (filter.GO.ids != FALSE){
    allTermName <- subset(allTermName, goid %in% filter.GO.ids)
  }
  if (!as.df){
    # as hash table
    allTermName <- hash(allTermName$goid, allTermName$name)
  }
  return(allTermName)
}

showSigNodes2 <-
  function(DAG, sigTerm, sigTerm_Local, sigTerm_Global, dagTermInfo) {
    
    require('Rgraphviz') || stop('package Rgraphviz is required')
    
    graphAttrs <- getDefaultAttrs(layoutType = 'dot')
    graphAttrs$cluster <- NULL
    graphAttrs$node$shape <- 'ellipse'
    graphAttrs$node$fontsize <- '20'
    
    nodeAttrs <- list()
    edgeAttrs <- list()
    
    allTerm <- as.character(dagTermInfo[,1])
    nodeAttrs$label[allTerm] <- allTerm
    
    rmLocalTerm <- setdiff(sigTerm, sigTerm_Local)
    nodeAttrs$color[rmLocalTerm] <- rep('red', length(rmLocalTerm))
    nodeAttrs$shape[rmLocalTerm] <- rep('circle', length(rmLocalTerm))
    
    rmGlobalTerm <- setdiff(sigTerm_Local, sigTerm_Global)
    nodeAttrs$color[rmGlobalTerm] <- rep('red', length(rmGlobalTerm))
    nodeAttrs$shape[rmGlobalTerm] <- rep('box', length(rmGlobalTerm))
    nodeAttrs$height[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
    nodeAttrs$width[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
    
    nodeAttrs$color[sigTerm_Global] <- rep('red', length(sigTerm_Global))
    nodeAttrs$shape[sigTerm_Global] <- rep('rectangle', length(sigTerm_Global))
    nodeAttrs$height[sigTerm_Global] <- rep('0.7', length(sigTerm_Global))
    nodeAttrs$width[sigTerm_Global] <- rep('1.1', length(sigTerm_Global))
    
    dagTermInfo[dagTermInfo[,5]<2.2E-16,5] <- 2.2E-16;
    dagTermInfo[dagTermInfo[,6]<2.2E-16,6] <- 2.2E-16;
    dagTermInfo$colorran <- round(log10(dagTermInfo[,6])-range(log10(dagTermInfo[,6]))[1] + 1)
    mm <- max(dagTermInfo$colorran)
    colorMap <- heat.colors(mm)
    nodeAttrs$fillcolor[allTerm] <- unlist(lapply(dagTermInfo$colorran, function(x) return(colorMap[x])))
    
    weightsList <- edgeWeights(DAG)
    to <- lapply(weightsList, names)
    from <- nodes(DAG)
    edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
    edge.weights <- unlist(weightsList)
    names(edge.weights) <- edge.names
    ##    0 for a is_a relation,  1 for a part_of relation
    edgeAttrs$color <- ifelse(edge.weights == 0, 'black', 'red')
    plot(DAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)
  }