# FitSaturationCurveFindProblemGenes.R

FitSaturationCurveFindProblemGenes <- function(slope0, fit.results.out, fit.select.results.out, array.adj.out, diagnostics.plot.out){
  # Fit saturation curve: 3 parameters --------------------------------------
  
  # x.factor: make fit not so close to saturation points
  # x.factor = 1: best fit
  # x.factor > 1: force saturation points farhter from data points
  x.factor <- 1.2
  # bg.factor: make fit not so close to background
  bg.factor <- 0.8
  
  # init vals: linear model
  # slope0 <- 0.3
  int0 <- 10
  
  # init out list
  fit.list <- vector(mode="list", length=length(common.genes))
  names(fit.list) <- common.genes
  
  for (gene in common.genes){
    R.exprs <- unlist(rna.seq.exprs.common.g[gene, ])
    M.exprs <- unlist(array.exprs.subset.common.g[gene, ])
    M.full <- unlist(array.exprs[gene, ])
    
    # lower and upper bounds: saturation model
    bmin <- 0
    kmin <- 0
    amin <- max(M.full) * x.factor
    bmax <- min(M.full) * bg.factor
    amax <- Inf
    kmax <- Inf
    # init vals: saturation model
    b0 <- 2^4
    k0 <- 2^6  # large to begin in linear regime
    a0 <- amin + 10  # expect it close to amin
    # lower and upper bounds: 
    slopemin <- 0
    intmin <- 0
    slopemax <- Inf
    intmax <- min(M.full) * bg.factor
    # get variance as weights
    M.var <- predict(fit.noise, M.exprs)
    # adjust M.var so all values less than 0 take on smallest non-negative number
    M.var.min <- min(M.var[which(M.var > 0)])
    # if M.var.min is infinity, probably RNA-Seq has no expression, default weights to 1
    if (is.infinite(M.var.min)){
      warning(paste("Gene:", gene, "...adjusting infinites to 1 for variance."))
      M.var.min <- 1
    }
    M.var[which(M.var < 0)] <- M.var.min
    weights <- 1 / M.var
    
    fits <- tryCatch({
      
      fit.saturation <- nls(M.exprs ~ b + (a * R.exprs) / (k + R.exprs),
                            algorithm = "port",
                            start=list(a=a0, 
                                       b=b0,
                                       k=k0),
                            lower=list(a=amin,
                                       b=bmin,
                                       k=kmin),
                            upper=list(a=amax,
                                       b=bmax,
                                       k=kmax),
                            weights=weights)
      fit.lm <- FitLmConstraint(M.exprs, R.exprs, weights, 
                                int0, slope0, 
                                intmin, slopemin, 
                                intmax, slopemax)
      
      fits <- list(saturation=fit.saturation,
                   lm=fit.lm)
      
    }, error = function(e) {
      fit.saturation <- NA
      # fit.lm <- lm(M.exprs ~ R.exprs)
      fit.lm <- FitLmConstraint(M.exprs, R.exprs, weights, 
                                int0, slope0, 
                                intmin, slopemin, 
                                intmax, slopemax)
      
      return(list(saturation=fit.saturation, 
                  lm=fit.lm))
      
    })
    fit.list[gene] <- list(fits)
  }
  
  # Save fit results to .RData ----------------------------------------------
  
  save(fit.list, file = fit.results.out)
  
  # F-test on saturation and linear fit ----------------------------------------
  
  fit.select.list <- vector(mode="list", length=length(common.genes))
  names(fit.select.list) <- common.genes
  for (gene in common.genes){
    # check if saturation fit was performed...
    fit.saturation <- fit.list[[gene]][["saturation"]]
    fit.lm <- fit.list[[gene]][["lm"]]
    # select either saturation or lm based on f.test or by default (i.e. saturation is NA)
    fit.select <- FTestSigmoidLinearModels(fit.lm, fit.saturation, pval=pval, complex.model="saturation")
    fit.select.list[[gene]] <- fit.select
  }
  
  save(fit.select.list, file = fit.select.results.out)
  
  # Plot clock genes: diagnostics -------------------------------------------
  
  pdf(diagnostics.plot.out)
  par(mfrow = c(2,1))
  for (gene in c(clockgenes, tissuegenes, problematicgenes)){
    fit.select <- fit.select.list[[gene]]
    fit.used <- fit.select$fit.used  # either saturation or lm
    myfit <- fit.select$myfit
    
    # Get vector of predicted values
    x <- GetFullR(gene, rna.seq.exprs, common.samples)
    y <- unlist(array.exprs[gene, ])
    # x.predict <- seq(min(x)*0.8, max(x)*1.2, length.out=10*length(x))
    y.predict <- seq(min(y)*0.8, max(y), length.out=10*length(y))
    if (fit.used == "saturation"){
      x.hat <- saturation.inv(coef(myfit), y.predict)
    } else if (fit.used == "lm"){
      x.hat <- linear.inv(coef(myfit), y.predict)
    } else {
      warning("Neither saturation nor lm")
    }
    
    # plot 
    params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
    symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
    sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
    plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
         xlab="RNA-Seq DESeq-normalized counts",
         ylab="Microarray normal scale")
    lines(x.hat, y.predict)
    plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
         xlab="RNA-Seq DESeq-normalized counts (log2)",
         ylab="Microarray log2")
    lines(log2(x.hat + 1), log2(y.predict))
    
  }
  dev.off()
  
  
  # Adjust microarray -------------------------------------------------------
  
  array.adj <- matrix(NA, nrow=nrow(array.exprs), ncol=ncol(array.exprs),
                      dimnames=list(rownames(array.exprs), 
                                    colnames(array.exprs)))
  
  for (gene in common.genes){
    fit.used <- fit.select.list[[gene]][["fit.used"]]
    myfit <- fit.select.list[[gene]][["myfit"]]  # either saturation or lm
    
    array.exprs.gene <- as.matrix(array.exprs[gene, ])
    
    if (fit.used == "saturation"){
      array.adj.gene <- saturation.inv(coef(myfit), array.exprs.gene)
      
    } else if (fit.used == "lm"){
      array.adj.gene <- linear.inv(coef(myfit), array.exprs.gene)
    }
    array.adj[gene, ] <- array.adj.gene
  }
  
  
  # Write to file -----------------------------------------------------------
  
  write.table(array.adj, file = array.adj.out, 
              quote=FALSE, sep='\t',
              row.names=TRUE, col.names=NA)
  
  
  # Check for problematic genes ---------------------------------------------
  
  # How many have negative values? ------------------------------------------
  
  negs <- apply(array.adj, 1, function(x){
    if (min(x) < 0){
      return(1)
    } else {
      return(0)
    }
  })
  
  problem.genes <- names(negs[which(negs == 1)])
  
  print(paste('slope:', slope0, 'has', length(problem.genes), 'problem genes.'))
  
  # Plot diagnostics for problem genes --------------------------------------
  
  pdf(diagnostics.plot.out)
  par(mfrow = c(2,1))
  for (gene in c(problem.genes)){
    fit.select <- fit.select.list[[gene]]
    fit.used <- fit.select$fit.used  # either saturation or lm
    myfit <- fit.select$myfit
    
    # Get vector of predicted values
    x <- GetFullR(gene, rna.seq.exprs, common.samples)
    y <- unlist(array.exprs[gene, ])
    # x.predict <- seq(min(x)*0.8, max(x)*1.2, length.out=10*length(x))
    y.predict <- seq(min(y)*0.8, max(y), length.out=10*length(y))
    if (fit.used == "saturation"){
      x.hat <- saturation.inv(coef(myfit), y.predict)
    } else if (fit.used == "lm"){
      x.hat <- linear.inv(coef(myfit), y.predict)
    } else {
      warning("Neither saturation nor lm")
    }
    
    # plot 
    params.str <- paste0(signif(as.vector(coef(myfit)), 2), collapse=",")
    symbols <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=8, obs=1)
    sizes <- GetUnobsObsSymbol(all.samples=colnames(array.exprs), common.samples, unobs=0.25, obs=1)
    plot(x, y, main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
         xlab="RNA-Seq DESeq-normalized counts",
         ylab="Microarray normal scale")
    lines(x.hat, y.predict)
    plot(log2(x + 1), log2(y), main=paste0("Gene=", gene, " Params=", params.str), pch=symbols, cex=sizes,
         xlab="RNA-Seq DESeq-normalized counts (log2)",
         ylab="Microarray log2")
    lines(log2(x.hat + 1), log2(y.predict))
    
  }
  dev.off()
  
  return(length(problem.genes))
}
