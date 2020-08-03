# AnalyzeGeneEnrichment.R
# using topGO
# 25 Feb 2014
#
# if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# 
# NOTE: if you run functions here in parallel (mclapply) you need to source these scripts INSIDE the loop:
# https://support.bioconductor.org/p/38541/

#library(DBI)  # dbGetQuery() not found when loading topGO
#library(topGO)
#library(org.Mm.eg.db)

# source("scripts/functions/DataHandlingFunctions.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/DataHandlingFunctions.R")  # use across different projects

MergeGOTerms <- function(enrichment, go.terms, new.go.term){
  # After running GetGOEnrichment, merge some GO terms that may be "similar"
  # after merging, does not recalculate pvalue, simply takes the worst p-value
  enrichment.goterms <- subset(enrichment, GO.ID %in% go.terms)
  enrichment.toappend <- subset(enrichment, !GO.ID %in% go.terms)
  # merge go terms
  enrichment.merged.base <- tbl_df(data.frame(GO.ID = NA, Term = NA, Annotated = NA, Significant = NA, Expected = NA, classicFisher = NA, FDRadj = NA, N.genes = NA, genes = NA, minuslogpval = NA))
  enrichment.merged <- tbl_df(data.frame(GO.ID = paste(enrichment.goterms$GO.ID, collapse = ","),
                                         Term = new.go.term,
                                         Annotated = NA,
                                         Signficiant = NA,
                                         Expected = NA,
                                         classicFisher = NA,
                                         FDRadj = max(enrichment.goterms$FDRadj),
                                         N.genes = unique(enrichment.goterms$N.genes),
                                         genes = NA,  # add genes later
                                         # genes = unique(do.call(c, enrichment.goterms$genes)),
                                         minuslogpval = min(enrichment.goterms$minuslogpval)))
  enrichment.merged <- bind_rows(enrichment.merged, enrichment.merged.base)  # dplyr doesnt handle single row tables well
  # add genes
  genes <- unique(do.call(c, enrichment.goterms$genes))
  genes.lst <- list(genes, NA)
  enrichment.merged$genes <- sapply(genes.lst, function(g) return(g))
  enrichment.merged <- subset(bind_rows(enrichment.merged, enrichment.toappend), !is.na(GO.ID))
  return(enrichment.merged)
}

GetGOEnrichment <- function(genes.bg, genes.fg, fdr.cutoff, show.top.n = 8, ontology="BP", wd = "/home/yeung/projects/tissue-specificity", filter.GO.terms=FALSE){
  source(file.path(wd, "scripts/functions/AnalyzeGeneEnrichment.R"))
  enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg, FDR.cutoff = fdr.cutoff, which.ontology = ontology, return.GOdata = TRUE, filter.GO.terms = filter.GO.terms)
  enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
  
  enrichment <- tryCatch({
    enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
  }, error = function(e) {
    enrichment <- enrichment
  })
  
  # enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
  show.top.n.min <- min(nrow(enrichment), show.top.n)
  if (show.top.n.min == 0) return(NULL)
  enrichment <- enrichment[1:show.top.n.min, ]   # prevent taking more than you have enrichment
  # unload packages
  # sometimes topGO causes problems (unable to unload later), unload once you're done.
  detach(name = "package:topGO", unload = TRUE)
  detach(name = "package:org.Mm.eg.db", unload = TRUE)
  return(enrichment)
}

PlotEnrichmentGenes <- function(dat.freq, enrichment, max.genes, row.i = "max"){
  # take hit with most genes
  # if row.i is max, get row with most genes
  # otherwise use row.i as index for row
  if (row.i == "max"){
    row.i <- which(enrichment$Significant == max(enrichment$Significant, na.rm = TRUE))
  } else if (is.numeric(row.i)){
    row.i <- row.i
  } else {
    warning(paste("row.i is not max or numeric:", row.i))
  }
  go.genes <- enrichment$genes[[row.i]]
  go.term <- enrichment$Term[[row.i]]
  show.n.genes <- min(length(go.genes), max.genes)
  s <- SvdOnComplex(subset(dat.freq, gene %in% go.genes), value.var = "exprs.transformed")
  eigens <- GetEigens(s, period = 24, comp = comp, label.n = show.n.genes, eigenval = TRUE, adj.mag = TRUE, constant.amp = 4, peak.to.trough = TRUE, label.gene = FALSE)
  print(eigens$u.plot)
}

PlotGeneModuleWithGO <- function(dat.sub, enrichment, jtitle = "", text.size = 6, dot.size = 6, comp = 1, jsize = 22, legend.pos = "bottom", label.gene = c(), top.hits = 25, label.GO.terms.only=FALSE, unlabeled.alpha = 0.3){
  # show.top.n > 8 recycles colours
  # enrichment from GetGOEnrichment()
  # annotate dat.freq with GO terms for each gene (take most significant enrichment)
  txtgray <- "gray70"
  
  constant.amp <- text.size
  ylab <- ""
  xlab <- "Log2 Fold Change"
  
  # label clocks to based on GOterm
  go.hash <- hash()
  for (i in seq(nrow(enrichment))){
    genes.vec <- enrichment$genes[[i]]
    term <- as.character(enrichment$Term[[i]])
    for (g in genes.vec){
      if (is.null(go.hash[[g]])){
        go.hash[[g]] <- term
      } 
    }
  }
  
  # dat.sub <- subset(dat.freq, gene %in% genes.fg)
  s.sub <- SvdOnComplex(dat.sub, value.var = "exprs.transformed")
  
  # label dat.freq with goterm
  dat.sub$term <- sapply(as.character(dat.sub$gene), function(g){
    if (!is.null(go.hash[[g]])){
      return(go.hash[[g]])
    } else {
      return(NA)
    }
  })
  
  eig <- GetEigens(s.sub, period = 24, comp = comp, label.n = 25, eigenval = TRUE, adj.mag = TRUE, constant.amp = constant.amp, peak.to.trough = TRUE)
  
  omega <- 2 * pi / 24
  ampscale <- 2
  vec.complex <- eig$eigensamp 
  labels <- names(vec.complex)
  
  dat <- data.frame(amp = Mod(vec.complex) * ampscale,
                    phase = ConvertArgToPhase(Arg(vec.complex), omega = omega),
                    label = labels)
  
  # add term BEFORE filtering out labels
  dat$term <- as.factor(sapply(as.character(dat$label), function(g){
    if (g == ""){
      return("")
    }
    if (!is.null(go.hash[[g]])){
      return(go.hash[[g]])
    } else {
      return("")
    }
  }))
  # make "" the first term for levels (so you get colours of dots wihtout the labels 
  if (any(dat$term == "")){
    jterms <- as.character(enrichment$Term)
    dat$term <- factor(as.character(dat$term), levels = c("", jterms))
  }
  # make points transparent if not labeled to a GO term, otherwise use a color
  dat$jalpha <- sapply(as.character(dat$term), function(g){
    if (g == ""){
      return(unlabeled.alpha)
    } else {
      return(1)
    }
  })
  
  if (!label.GO.terms.only){
    top.amps <- as.character(head(dat[order(dat$amp, decreasing = TRUE), ], n = top.hits)$label)
  } else {
    dat.sub2 <- subset(dat, term != "")
    top.amps <- as.character(head(dat.sub2[order(dat.sub2$amp, decreasing = TRUE), ], n = top.hits)$label)
  }
  dat$label <- sapply(as.character(dat$label), function(l) ifelse(l %in% top.amps | l %in% label.gene, yes = l, no = ""))
  # label only top genes
  amp.max <- ceiling(max(dat$amp) * 2) / 2
  if (amp.max <= 1){
    amp.step <- 0.5
  } else {
    amp.step <- 1
  }
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#9F1C34", "#CC79A7")
  cbPalette <- cbPalette[1:length(levels(dat$term))]
  # change "" term to light gray
  cbPalette[levels(dat$term) == ""] <- txtgray
  
  m <- ggplot(data = dat, aes(x = amp, y = phase, label = label, colour = term, alpha = jalpha)) + 
    geom_point(size = dot.size, shape = 15) +
    coord_polar(theta = "y") + 
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(jtitle) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
    scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
    theme_bw(jsize) + 
    geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
    geom_hline(yintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position=legend.pos,
          panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank()) + 
    scale_alpha_identity()
  # add text
  df.txt <- subset(dat, label != "")
  if (constant.amp != FALSE){
    m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, label = label), size = constant.amp)
  } else {
    m <- m + geom_text_repel(data = df.txt, aes(x = amp, y = phase, size = amp, label = label))
  }
  return(m + scale_colour_manual(values=cbPalette) + ggtitle(jtitle))
}

CreateSym2Entrez <- function(){
  sym2entrez <- as.list(org.Mm.egALIAS2EG)
  sym2entrez <- sym2entrez[!is.na(sym2entrez)]
  return(sym2entrez)
}

CreateEntrez2Sym <- function(){
  entrez2sym <- as.list(org.Mm.egSYMBOL)
  entrez2sym <- entrez2sym[!is.na(entrez2sym)]
  return(entrez2sym)
}


CreateEntrez2GO <- function(){
  entrez2GO.obj <- org.Mm.egGO
  mapped.genes <- mappedkeys(entrez2GO.obj)
  entrez2GO.full <- as.list(entrez2GO.obj[mapped.genes])
  # collapse sublist of GO-IDs to vector
  entrez2GO <- lapply(entrez2GO.full, function(x){
    return(names(x))  # lapply is nice
  })
  return(entrez2GO)
}

CreateGO2AllEntrez <- function(){
  go2entrez <- as.list(org.Mm.egGO2ALLEGS)
  go2entrez <- go2entrez[!is.na(go2entrez)]
  return(go2entrez)
}

CreateGO2Entrez <- function(){
  go2entrez <- as.list(org.Mm.egGO2EG)
  go2entrez <- go2entrez[!is.na(go2entrez)]
  return(go2entrez)
}

ConvertGO2Entrez <- function(go.list, go2entrez){
  gene.list.entrez <- sapply(go.list, function(g){
    entrez <- go2entrez[[g]]
    if (!is.null(entrez)){
      return(entrez)
    }
  }, simplify = TRUE)
  names(gene.list.entrez) <- go.list
  # print(head(gene.list.entrez))
  gene.list.entrez <- gene.list.entrez[!sapply(gene.list.entrez, is.null)]
  # print(head(gene.list.entrez))
  return(unlist(gene.list.entrez, use.names=TRUE))

}

ConvertSym2Entrez <- function(gene.list, sym2entrez){
  
  # convert to entrez, auto filtering here
  gene.list.entrez <- sapply(gene.list, function(g){
    entrez <- sym2entrez[[g]]
    if (!is.null(entrez)){
      return(entrez)
    }
  }, simplify = TRUE)
  names(gene.list.entrez) <- gene.list
  # print(head(gene.list.entrez))
  gene.list.entrez <- gene.list.entrez[!sapply(gene.list.entrez, is.null)]
  # print(head(gene.list.entrez))
  return(unlist(gene.list.entrez, use.names=TRUE))
}

ConvertEntrez2Sym <- function(gene.list, entrez2sym){
  gene.list.sym <- sapply(gene.list, function(g){
    sym <- entrez2sym[[g]]
    if (!is.null(sym)){
      return(sym)
    }
  }, simplify = TRUE)
  names(gene.list.sym) <- gene.list
  # print(head(gene.list.entrez))
  gene.list.sym <- gene.list.sym[!sapply(gene.list.sym, is.null)]
  # print(head(gene.list.entrez))
  return(unlist(gene.list.sym, use.names=TRUE))
}

AnalyzeGeneEnrichment <- function(genes.bg, genes.hit, 
                                  sym2entrez,
                                  entrez2GO,
                                  convert.sym.to.entrez = TRUE,
                                  which.ontology = "BP", 
                                  write.path = FALSE,
                                  node.size = 5, 
                                  FDR.cutoff = 0.05,
                                  return.GOdata = FALSE,
                                  debug=FALSE,
                                  filter.GO.terms=FALSE){
  # Analyze gene enrichment given background and hit genes.
  #
  # INPUT:
  # genes.bg: vector of background genes
  # genes.hit: vector of hit genes
  # sym2entrez: Provide a mapping from gene symbol to entrez, see example via CreateSym2Entrez(). If missing, creates it via CreateSym2Entrez()
  # entrez2GO: mapping from entrez to GO terms, example CreateEntrez2GO(). If missing, creates it via CreateEntrez2Go()
  # convert.sym.to.entrez: convert gene symbols to entrez ID. Set to FALSE if genes already in entrez ID
  # write.path: path to write topGO results object. If FALSE, does not write to file
  # which.ontology: "BP", "MF", "CC"
  # node.size: prune GO hierarchy from terms which have less than node.size genes
  # FDR.cutoff: Fisher's exact test FDR-adjusted pvals cutoff to be significant
  # write.path: if false jsut return object, if true also write to file given by write.path
  #
  # OUTPUT:
  # all.res: topGO results object
  #
  if (missing(sym2entrez)){
    sym2entrez <- CreateSym2Entrez()
  }
  if (missing(entrez2GO)){
    entrez2GO <- CreateEntrez2GO()
  }
  # if (all(filter.GO.terms != FALSE)){
  #   entrez2GO <- lapply(entrez2GO, function(GOlist) GOlist[which(GOlist %in% filter.GO.terms)])
  # }
  
  if (convert.sym.to.entrez){
    genes.bg <- ConvertSym2Entrez(genes.bg, sym2entrez)
    genes.hit <- ConvertSym2Entrez(genes.hit, sym2entrez)
  }
  
  # Select genes in bg that are in hit, binary matrix.
  # used in topGO new() function
  sel.genes <- factor(as.integer(genes.bg %in% genes.hit))
  names(sel.genes) <- genes.bg
  
  # Get topGO object
  GOdata <- new("topGOdata",
                ontology = which.ontology,
                allGenes = sel.genes,
                nodeSize = node.size,
                annot = annFUN.gene2GO,
                gene2GO = entrez2GO)
  if (debug){
    return(GOdata)
  }
  
  # Run enrichment ----------------------------------------------------------
  
  result.fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  n.tests <- length(score(result.fisher))
  
  # Look at results ---------------------------------------------------------
  
  all.res <- GenTable(GOdata, 
                      classicFisher = result.fisher,
                      orderBy = "classicFisher",
                      ranksOf = "classicFisher",
                      topNodes = n.tests)
  # adjust pvalues
  all.res$FDRadj <- p.adjust(all.res$classicFisher, method = "BH")
  
  
  # Filter table by pvalue --------------------------------------------------
  
  print(paste(nrow(all.res), "GO.terms found enriched before filtering by pval"))
  all.res <- all.res[all.res$FDRadj <= FDR.cutoff, ]
  print(paste(nrow(all.res), "GO.terms found enriched."))
  
  # optionally filter by GO term
  if (all(filter.GO.terms != FALSE)){
    all.res <- subset(all.res, GO.ID %in% filter.GO.terms)
  }
  if (nrow(all.res) < 1){
    print("No enrichment found below cutoff")
    return(all.res)  # just skip it and return an empty DF
  } else if (nrow(all.res) == 1){
    # with just one row, we cannot add a vector of genes into a single
    # column
    all.res <- rbind(all.res, rep(NA, ncol(all.res)))
    all.res <- as.tbl(all.res)
  }
  
  # Add all genes that were considered (can recapitulate contTable) ---------
  N.genes <- length(genes(GOdata))
  
  
  all.res$N.genes <- N.genes
  
  # Write to file -----------------------------------------------------------
  if (write.path != FALSE){
    write.table(all.res, file = write.path, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE) 
  }
  if (!return.GOdata){
    return(all.res)
  } else {
    # https://support.bioconductor.org/p/65856/#66677
    # add gene names to all.res
    my.genes <- sapply(all.res$GO.ID, function(x){
      if (!is.na(x)){
        genes <- genesInTerm(GOdata, x) 
        genes[[1]][genes[[1]] %in% genes.hit]
      } else {
        return(NA)
      }
    })
    # convert to gene symbol
    entrez2sym <- as.list(org.Mm.egSYMBOL)
    entrez2sym <- entrez2sym[!is.na(entrez2sym)]
    all.res$genes <- sapply(my.genes, function(g){
      if (length(g) == 0){
        return(NA)
      }
      if (!is.na(g)){
        return(unlist(entrez2sym[g], use.names = FALSE))
      } else {
        return(NA)
      }
    })
    return(all.res)
    # return(list(res=all.res, godata=GOdata, genes=genes.hit))
  }
}




EnrichmentBinnedToFile <- function(sorted.hits, bin.vector, fname.base, sym2entrez, entrez2GO){
  # Take sorted hits and take top bin.vector[i] genes for enrichment
  # to see how enrichment of certain processes evolve over time
  # 
  # sorted.hits: list of genes sorted by "relevance" 
  # bin.vector: vector of integers to bin genes (take top 10 genes, then top 20 ... etc )
  # fname.base: full path of filename, we will add _i_Ontology.GOtop to end of file
  
  for (i in bin.vector){
    genes.hit <- sorted.hits[1:i]
    for (onto in c("BP", "MF", "CC")){
      fname.out <- paste0(fname.base, '_', i, '_', onto, '.GOtop')
      res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                   sym2entrez, entrez2GO, 
                                   convert.sym.to.entrez = TRUE, 
                                   which.ontology = onto, 
                                   write.path = fname.out, 
                                   node.size = 5, 
                                   FDR.cutoff = 0.05)
    }
  }
}

GetGoObject <- function(genes.hit, genes.bg, onto, node.size, sym2entrez, entrez2GO){
  # Run GO
  
}

GetEnrichment <- function(sorted.hits, bin, onto, sym2entrez, entrez2GO){
  # Same as getEnrichmentOverBins but for just one bin, that way we can parallelize it 
  genes.hit <- sorted.hits[1:bin]  
  res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                               sym2entrez, entrez2GO, 
                               convert.sym.to.entrez = TRUE, 
                               which.ontology = onto, 
                               write.path = FALSE, 
                               node.size = 5, 
                               FDR.cutoff = 1)
  return(res)
}

GetEnrichmentParallel <- function(sorted.hits, bin.vector, onto, sym2entrez, entrez2GO, n.cores = 4){
  # Run GetEnrichment in parallel over all bins in bin.vector
  library(parallel)
  print("Running GO (~4 minutes)")
  start <- Sys.time()
  res.split <- mclapply(bin.vector, function(bin){
    genes.hit <- sorted.hits[1:bin]
    res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                 sym2entrez, entrez2GO, 
                                 convert.sym.to.entrez = TRUE, 
                                 which.ontology = onto, 
                                 write.path = FALSE, 
                                 node.size = 5, 
                                 FDR.cutoff = 1)
    res$bin <- bin
    return(res)
  }, mc.cores = n.cores)
  print(Sys.time() - start)
  res.all <- do.call(rbind, res.split)
  return(res.all)
}

GetEnrichmentOverBins <- function(sorted.hits, bin.vector, onto, go.terms, sym2entrez, entrez2GO){
  # Iterate over bins and find enrichment for a specific go term
  # track this enrichment as we sweep across the bins
  
  # init output matrix: 8 columns ("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "N.genes", "bin")
  res.out <- matrix(NA, nrow = length(bin.vector) * length(go.terms), ncol = 9)
  colnames(res.out) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "N.genes", "bin")
  
  rowcount <- 1
  for (i in 1:length(bin.vector)){
    n.genes <- bin.vector[i]
    
    genes.hit <- sorted.hits[1:n.genes]
    
    res <- AnalyzeGeneEnrichment(sorted.hits, genes.hit, 
                                 sym2entrez, entrez2GO, 
                                 convert.sym.to.entrez = TRUE, 
                                 which.ontology = onto, 
                                 write.path = FALSE, 
                                 node.size = 5, 
                                 FDR.cutoff = 1)
    res.sub <- subset(res, Term %in% go.terms)
    res.sub$bin <- n.genes
    # make into matrix, makes life easier later
    res.sub <- as.matrix(res.sub)
    res.sub[, 3:ncol(res.sub)] <- as.numeric(res.sub[, 3:ncol(res.sub)])
    
    rowstart <- rowcount
    rowend <- rowstart + length(go.terms) - 1
    res.out[rowstart:rowend, ] <- res.sub
    rowcount <- rowcount + length(go.terms)
  }
#     # make to dataframe and unfactor things
#   res.out <- data.frame(noquote(res.out))
#   for (cname in c("Annotated", "Significant", "Expected", "classicFisher", "FDRadj")){
#     res.out[[cname]] <- as.numeric(res.out[[cname]])
#   }
#   res.out <- data.frame(res.out)
#   for (cname in c("Annotated", "Significant", "Expected", "classicFisher", "FDRadj", "bin")){
#     res.out[[cname]] <- as.numeric(as.character(res.out[[cname]]))
#   }
  return(res.out)
}
