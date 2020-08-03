#library(Matrix)
#library(dplyr)
#library(hash)

UnzipMatrices <- function(chunk.path, outdir){
  # Unzip matrix and save with extension ".full.Robj"
  load(chunk.path)  # des.mats.list
  outname <- basename(chunk.path)
  outname <- paste(strsplit(outname, ".Robj")[[1]], collapse = ".")
  outname <- paste0(outname, ".unzipped.Robj")
  out.path <- file.path(outdir, outname)
  
  des.mats.list <- lapply(des.mats.list, function(des.mats){
    des.mats$mat <- as.matrix(des.mats$mat)
    return(des.mats)
  })
  return(des.mats.list)
  save(des.mats.list, file = out.path)
}

ConcatenateFits <- function(genedir){
  source("scripts/functions/AppendListFunctions.R")
  start <- Sys.time()
  fits.long.list <- expandingList()
  for (gene in list.files(genedir)){
    fitdir <- file.path(genedir, gene)
    fits <- LoadFitsFromModels(fitdir)
    fits.top <- GetTopNModelsFromFits(fits, top.n)
    fits.long.gene <- ListToLong(fits.top, genename = gene, top.n = top.n, period = 24)
    fits.long.list$add(fits.long.gene)
    #   break
  }
  fits.long.list <- fits.long.list$as.list()
  fits.long <- do.call(rbind, fits.long.list)
  return(fits.long)
}

LoadFitsFromModels <- function(fitdir){
  fits.list <- expandingList()
  fit.files <- list.files(fitdir)
  for (i in seq(length(fit.files))){
    load(file.path(fitdir, fit.files[[i]]))  # fits
    fits.list$add(fits)
  }
  return(fits.list$as.list())
}

AddDf <- function(obj){
  # add arbitrary object to last column of data frame
  return(obj)
}

FitToMat <- function(fit, genename, weight, weight.raw, period = 24){
  # Turn fit into a long data frame
  model.name <- CoefToModelName(fit)
  # vector to long
  params <- CoefToParams(fit)
  dat.out <- data.frame(gene = genename, model = model.name, weight = weight, weight.raw = weight.raw)  # init
  dat.params <- dat.out %>%
    group_by(gene, model, weight, weight.raw) %>%
    do(param.list = AddDf(params))
  return(dat.params)
}

ListToLong <- function(fits.list, genename, top.n, period = 24){
  # form list of fits (sorted by best to worst) into a long dataframe
  # should be easily used to append or edit a data frame
  mats <- list()
  for (i in seq(top.n)){
    mats[[i]] <- FitToMat(fits.list[[i]]$fit, genename, fits.list[[i]]$weight.norm, fits.list[[i]]$weight.raw, period)
  }
  return(do.call(rbind, mats))
}

GetTopNModelsFromFits <- function(fits.list, top.n){
  fits.list <- NormalizeWeightsFromFits(fits.list)  # handles pesky fit.weights.sum  
  top.n.weights <- vector("numeric", length = top.n)
  top.n.tuples <- matrix(data = 0, nrow = top.n, ncol = 2)  # 2 because i and j indices
  for (i in seq(length(fits.list))){
    # if (is.null(fits.list[[i]])) next
    for (j in seq(length(fits.list[[i]]))){  # should be list of n.top or more
      weight.j <- fits.list[[i]][[j]]$weight
      min.indx <- which.min(top.n.weights)
      # if (length(weight.j) == 0) next
      if (is.null(fits.list[[i]][[j]])) next
      # print(weight.j)
      # print(top.n.weights[min.indx])
      # print(fits.list[[i]][[j]])
      if (top.n.weights[min.indx] < weight.j){
        # update weights and indcies
        top.n.weights[min.indx] <- weight.j
        # print(top.n.tuples[min.indx, c(1, 2)])
        top.n.tuples[min.indx, c(1, 2)] <- c(i, j)
        # print(top.n.tuples[min.indx, c(1, 2)])
      }
    }
  }
  # order top.n.tuples by top.n.weights
  top.n.tuples <- top.n.tuples[order(top.n.weights, decreasing = TRUE), ]
  
  # get the top.n tuples
  fits.top <- expandingList()
  # print(top.n.tuples)
  apply(top.n.tuples, MARGIN = 1, function(tup){
    # print(tup)
    fits.top$add(fits.list[[tup[1]]][[tup[2]]])
  })
  return(fits.top$as.list())
}

NormalizeWeightsFromFits <- function(fits.list){
  fit.weight.sum <- SumWeights(fits.list)
  for (i in seq(length(fits.list))){
    fits.list[[i]][["fit.weights.sum"]] <- NULL
    for (j in seq(length(fits.list[[i]]))){
      if (is.null(fits.list[[i]][[j]])) next
      fits.list[[i]][[j]][["weight.norm"]] <- fits.list[[i]][[j]][["weight"]] / fit.weight.sum
    }
  }
  return(fits.list)
}

SumWeights <- function(fits.list){
  # Given list of fits, divide each weight by fit.weights.sum
  fit.weight.sum <- 0
  for (fit in fits.list){
    if (!is.null(fit$fit.weights.sum)){
      fit.weight.sum <- fit.weight.sum + fit$fit.weights.sum
    }
  }
  return(fit.weight.sum)
}

GetNrhythFromModel <- function(model){
  # from model, return number of tissues that are rhythmic 
  # e.g.: Mus;Adr,Kidney,Liver;Aorta,BFAT,Heart,Lung = 8 rhythmic tissues
  # strategy: replace ";" with ',', then get length of ',' split string
  n.rhyth <- gsub(pattern = ";", replacement = ",", x = model)
  return(length(strsplit(n.rhyth, ",")[[1]]))
}

GetAvgAmpFromParams <- function(params, by.model=FALSE){
  # look at param names to get average amplitudes
  # option: by.model: go strictly by model (used if model was filtered by some amp criteria beforehand)
  # by.model is a model name e.g., "Kidney;Liver;Adr,Lung" or FALSE
  amps <- params[grepl("amp", names(params))]
  if (by.model){
    if (by.model == "") return(NA)  # non rhyth model
    grepstr = gsub(pattern = ";", replacement = "|", as.character(by.model))
    amps <- amps[grepl(grepstr, names(amps))]
  }
  if(length(amps) == 0) return(NA)  # non rhyth model
  
  # get number of tissues rhythmic for each rhythmic paramter (because tissues can
  # share paramters)
  n.tiss <- sapply(names(amps), function(tiss.id) length(strsplit(tiss.id, ",")[[1]]))
  amps.weighted <- sum(amps * n.tiss)
  # return weighted mean
  return(amps.weighted / sum(n.tiss))
}

GetSdPhaseFromParams <- function(params, by.model=FALSE){
  # look at param names to get average amplitudes
  # option: by.model: go strictly by model (used if model was filtered by some amp criteria beforehand)
  # by.model is a model name e.g., "Kidney;Liver;Adr,Lung" or FALSE
  phases <- params[grepl("phase", names(params))]
  if (by.model){
    if (by.model == "") return(NA)  # non rhyth model
    grepstr = gsub(pattern = ";", replacement = "|", as.character(by.model))
    phases <- phases[grepl(grepstr, names(phases))]
  }
  if(length(phases) == 0) return(NA)  # non rhyth model
  
  # get number of tissues rhythmic for each rhythmic paramter (because tissues can
  # share paramters)
  n.tiss <- sapply(names(phases), function(tiss.id) length(strsplit(tiss.id, ",")[[1]]))
  weights = n.tiss / sum(n.tiss)
  phases.avg = weighted.mean(phases, weights)
  # phases.weighted = sum(phases * n.tiss); phases.avg = phases.weighted / sum(n.tiss)  # equivalent but less readable
  
  phases.var =  sum(weights * (DiffPhase(phases, phases.avg)) ^ 2)  # denominator is not N-1 it is N
  return(sqrt(phases.var))
}

GetPhaseFromParams <- function(params, by.model=FALSE){
  # works only for n.params == 1
  phases <- params[grepl("phase", names(params))]
  if (by.model){
    if (by.model == "") return(NA)  # non rhyth model
    grepstr = gsub(pattern = ";", replacement = "|", as.character(by.model))
    phases <- phases[grepl(grepstr, names(phases))]
  }
  if(length(phases) == 0) return(NA)  # non rhyth model
  return(mean(phases))
}

DiffPhase <- function(phase1, phase2, period=24){
  # Get phase1 and phase2 difference, modulo period
  jdiff <- abs(diff(c(phase1, phase2)))
  return(min(period - jdiff, jdiff))
}

MaxDiffPhaseVector <- function(phases, period=24){
  phasediff.max <- 0
  for (phase.i in phases){
    for (phase.j in phases){
      jdiff <- DiffPhase(phase.j, phase.i, period)
      if (jdiff > phasediff.max){
        phasediff.max <- jdiff
      }
    }
  }
  return(phasediff.max)
}

GetMaxPhaseDiffFromParams <- function(params, by.model=FALSE, period=24){
  # Get maximum phase difference between tissues
  phases <- params[grepl("phase", names(params))]
  if (by.model != FALSE){
    if (by.model == "") return(NA)  # non rhyth model
    grepstr = gsub(pattern = ";", replacement = "|", as.character(by.model))
    phases <- phases[grepl(grepstr, names(phases))]
  }
  if(length(phases) == 0) return(NA)  # non rhyth model
  phasediff.max <- MaxDiffPhaseVector(phases, period)
  return(phasediff.max)
}

FilterModelByAmp <- function(model.name, params, amp.cutoff = 0.5){
  # http://stackoverflow.com/questions/10049402/calculating-weighted-mean-and-standard-deviation
  # Get amplitude of each tissue in model name.
  # model name: Kidney;Liver;Adr,Lung 
  # params: fit from nconds2: amp and phase are called tissue1.amp tissue1.phase ... 
  # return:
  #   model.filt: Kidney:Liver if Adr,Lung amp parameter is less than amp.cutoff.
  model.vec <- strsplit(as.character(model.name), split = ";")[[1]]
  amps <- params[grepl("amp", names(params))]
  # remove tissues from model with low amplitude
  amps.filt <- amps[which(amps >= amp.cutoff)]
  
  # create new model name based on filtering
  # from c("Mus.amp", "Adr,Kidney,Liver.amp" ... ) -> c("Mus", "Adr,Kidney,Liver" ... )
  models.filt <- sapply(names(amps.filt), function(name) strsplit(name, ".amp")[[1]])
  models.filt <- paste(models.filt, collapse = ";")
  return(as.factor(models.filt))
}


MakeDesMatFromModelName <- function(dat.gene, model.name, tissues, w = 2 * pi / 24){
  # Given model.name, return design matrix
  # model.name: Aorta,BFAT,Liver;Adr,BS,Cere,Heart,Kidney
  des.mat <- GetFlatModel(dat.gene)  # cbind to this to create matrix
  tiss.combos <- GetAllCombos(tissues)
  des.mat.sinhash <- GetSinCombos(dat.gene, w, tissues, tiss.combos)
  des.mat.coshash <- GetCosCombos(dat.gene, w, tissues, tiss.combos)
  
  model.name.vec <- strsplit(model.name, ";")[[1]]
  for (rhyth.tiss in model.name.vec){
    des.mat <- cbind(des.mat, AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, rhyth.tiss))
  }
  return(des.mat)
}

LoadDesMatDatGeneRunFits <- function(dat.gene, mat.chunk, criterion = "BIC", normalize.weights = TRUE, top.n = 10, sparse = TRUE){
  # Run fits given dat.gene and chunk of mat, can be loaded from file 
  # track top.n
  # each element in list of mat.chunk has $mat, $rhyth.tiss etc. So use $mat for fitting
  
  # store my fits
  # to save memory, only store top n fits
  # fits <- expandingList()
  fits <- list()
  fit.count <- 0
  
  # track weights to normalize afterwards
  fit.weights.sum <- 0
  fit.weights <- vector(mode="numeric", length=top.n)
  
  # fit my model
  for (i in seq(length(mat.chunk))){
    des.mat <- mat.chunk[[i]]
    if (sparse){
      fit <- FitModel(dat.gene, as.matrix(des.mat$mat), get.criterion = criterion, condensed = TRUE)  # condensed will save some memory
    } else {
      fit <- FitModel(dat.gene, des.mat$mat, get.criterion = criterion, condensed = TRUE)  # condensed will save some memory
    }
    fit.weights.sum <- fit.weights.sum + fit$weight  # track even if it is a bad model, helps normalization
    fit.weights.lst <- UpdateFitWeights(fit$weight, fit.weights)
    if (!is.na(fit.weights.lst$i)){
      fits[[fit.weights.lst$i]] <- fit
      fit.weights <- fit.weights.lst$weights
    }
    fit.count <- fit.count + 1
  }
  
  if (normalize.weights){
    # print("Normalizing weights...")
    # print(dat.gene$gene[[1]])
    # fits.list <- NormalizeWeights(fits.list, cutoff)
    for (i in seq(length(fits))){
      fits[[i]]$weight.norm <- fits[[i]]$weight / fit.weights.sum
    }
  }
  fits$fit.weights.sum <- fit.weights.sum
  return(fits)
}

MakeDesMatChunks <- function(dat.gene, out.dir, tissues, n.rhyth.max, w = 2 * pi / 24, sparse = TRUE, chunks=10000, only.n.params = FALSE, count.mode=FALSE){
  # dat.gene: long format of gene expression, time and
  # conditions. Can include additional factors such as 
  # experiment.
  # to reduce computation, savee matrices into chunks (user defined, let's say 10000 models a chunk) so
  # if we run genome-wide, then we load the chunk and fit.
  #
  # only.n: store models only with exactly n parameters
  # count.mode: only care about counts, dont save anything
  
  if (missing(tissues)){
    tissues <- unique(as.character(dat.gene$tissue))
  }
  
  if (missing(n.rhyth.max)){
    n.rhyth.max <- length(tissues)
  } else if (n.rhyth.max < 1){
    print(n.rhyth.max)
    warning("N rhyth max cannot be less than 2")
  }
  
  tiss.combos <- GetAllCombos(tissues, ignore.full = FALSE)
  my_mat.queue <- new.queue()
  
  # BEGIN: init with flat model
  des.mat.flat <- GetFlatModel(dat.gene)
  
  if (sparse){
    # make sparse using Matrix package
    des.mat.flat <- Matrix(des.mat.flat)
  }
  
  # get rhythmic parameters which will be used for adding later: hash structure has fast lookup
  des.mat.sinhash <- GetSinCombos(dat.gene, w, tissues, tiss.combos)
  des.mat.coshash <- GetCosCombos(dat.gene, w, tissues, tiss.combos)
  
  rhyth.tiss <- list(character(0))  # needs to track shared and independent parameters, e.g.: c("Liver,Kidney", "Adr") no duplicates allowed
  # n.rhyth <- NRhythmicFromString(rhyth.tiss)  # number of independent rhythmic parameters perhaps? Do this later for speed?
  complement <- FilterCombos(tiss.combos, rhyth.tiss)  # given current matrix, you will know which tissues to iterate
  des.mat.list <- list(mat=des.mat.flat, rhyth.tiss=rhyth.tiss, complement = complement)
  # END: init with flat model
  
  # load up my queue
  enqueue(my_mat.queue, des.mat.list)
  
  # uncomment if you want to store all the matrices were used
  des.mats <- expandingList() 
  des.mats$add(des.mat.list)
  total.count <- 1  # mats total, to track how many models
  chunk.id <- 1  # increment every time we save a model
  
  # need to track models that we have done, so we eliminate "permutations" like c("Liver", "Kidney") and c("Kidney", "Liver) models
  # use hash for speed
  models.done <- hash()
  
  # generate matrix by adding combinations of columns and adding
  # those matrices into the queue
  while (! is.empty(my_mat.queue)) {
    des.mat.list <- dequeue(my_mat.queue)
    
    # check that this matrix is not already the maximum complexity (n.rhyth.max),
    # if it is already as complex as we want, then ignore it because 
    # we dont want to add another rhythmic column to this.
    #
    # strictly greater than should handle the case of "no rhytmic tissues"
    if (length(des.mat.list$rhyth.tiss) > n.rhyth.max){
      # print(paste("Skipping", des.mat.list$rhyth.tiss))
      next  # should work even in n.rhyth.max == length(tissues)
    }
    
    # determine tissue combinations that need to be added based on rhyth.tiss
    # e.g., no need to add Liver twice, they can't have two rhythmic paramters  
    for (tiss.comb in des.mat.list$complement){
      # add column for each tissue combination
      tiss.key <- paste(tiss.comb, collapse = ",")
      
      # append tiss.key to rhyth.tiss
      rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)  # form list("Adr,Kidney", "Mus")
      
      # check if this tissue combination has been already submitted into queue (but in different permutation)
      # track models we have done globally
      modelname <- MakeModelName(rhyth.tiss)
      if (! is.null(models.done[[modelname]])){
        # this is a permutation of an already done combo, skip
        #       print(rhyth.tiss)
        #       print(paste('Skipping', modelname))
        next
      }
      
      col.new <- AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, tiss.key)
      
      #     rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)
      
      # further remove complement after having
      tiss.complement.new <- FilterCombos(des.mat.list$complement, tiss.comb)
      
      # make new matrix, put it into queue
      mat.new <- cbind(des.mat.list$mat, col.new)
      des.mat.list.new <- list(mat=mat.new, rhyth.tiss = rhyth.tiss, complement = tiss.complement.new)
      enqueue(my_mat.queue, des.mat.list.new)
      models.done[[modelname]] <- TRUE  # we dont want to redo permutations of same models
      total.count <- total.count + 1
      # add to list unless only.n is not FALSE
      if (only.n.params == FALSE){
        des.mats$add(des.mat.list.new)
      } else {
        # only add if only.n.params matches n.rhyth.params
        n.rhyth.params <- length(rhyth.tiss) - 1  # rhyth.tiss contains a flat model, remove 1 to consider only n.rhyth.params
        if (n.rhyth.params == only.n.params){
          des.mats$add(des.mat.list.new)
        }
      }
      if (total.count %% chunks == 0){
        if (count.mode == FALSE){
          SaveChunk(chunk.id, out.dir, des.mats)
        }
        # start a new des.mats
        des.mats <- expandingList()
        chunk.id <- chunk.id + 1
      }
    }  
  }
  # save last batch
  if (count.mode == FALSE){
    SaveChunk(chunk.id, out.dir, des.mats)
  }
  print(paste("Total models:", total.count))
  return(NULL)
}

SaveChunk <- function(chunk.id, out.dir, des.mats){
  fname <- paste0("chunk.", chunk.id, ".Robj")
  des.mats.list <- des.mats$as.list()
  save(des.mats.list, file = file.path(out.dir, fname))
  return(NULL)
}

LoadChunkRunNconds <- function(chunk.path, genename, tissues, n.rhyth.max, w, criterion, normalize.weights, cutoff, top.n, sparse){
  # given path to a dat.gene chunk, run nconds2
  if (missing(genename)){
    genename <- strsplit(basename(chunk.path), ".Robj")[[1]][[1]]
  }
  load(chunk.path)
  return(MakeDesMatRunFit(dat.gene, genename, tissues, n.rhyth.max, w, criterion, normalize.weights, cutoff, top.n, sparse))
}

ChunkDatGenesToFile <- function(dat.long, write.dir){
  # Split dat.long by gene, save each as an genename.Robj
  dat.long %>%
    group_by(gene) %>%
    do(SaveDatGenes(., write.dir))
  return(NULL)
}

SaveDatGenes <- function(dat.gene, write.dir){
  gene <- dat.gene$gene[[1]]
  fname <- paste0(gene, ".Robj")
  save(dat.gene, file = file.path(write.dir, fname))
  return(data.frame(NULL))
}

GetBestModel <- function(fits){
  fit.bestweight <- 0
  for (i in seq(length(fits))){
    if (length(fits[[i]]) == 1){
      # skip because one of the parameters is just gene name
      next
    }
    if (fits[[i]]$weight.norm > fit.bestweight){
      fit.best <- fits[[i]]
      fit.bestweight <- fits[[i]]$weight.norm
    }
  }
  return(fit.best)
}

MakeDesMatRunFitEnv <- function(env, genename, tissues, n.rhyth.max, w = 2 * pi / 24, criterion = "BIC", normalize.weights = FALSE, cutoff = 1e-3, top.n = 3, sparse = TRUE){
  # wrapper function to grab dat.gene from environment
  # genename is gene name found from ls(env)
  dat.gene <- get(genename, envir = env)
  return(MakeDesMatRunFit(dat.gene, genename, tissues, n.rhyth.max, w, criterion, normalize.weights, cutoff, top.n))
}

MakeDesMatRunFit <- function(dat.gene, gene, tissues, n.rhyth.max, w = 2 * pi / 24, criterion = "BIC", normalize.weights = FALSE, cutoff = 1e-3, top.n = 3, sparse = TRUE){
  # dat.gene: long format of gene expression, time and
  # conditions. Can include additional factors such as 
  # experiment.
  # To save memory, generate model and then run fit immediately afterwards.
  # optionally allow stopping after models reach a certain complexity
  
  if (missing(tissues)){
    tissues <- unique(as.character(dat.gene$tissue))
  }
  if (missing(gene)){
    gene <- as.character(dat.gene$gene[[1]])
  }
  
  if (missing(n.rhyth.max)){
    n.rhyth.max <- length(tissues)
  } else if (n.rhyth.max < 1){
    print(n.rhyth.max)
    warning("N rhyth max cannot be less than 2")
  }
  tiss.combos <- GetAllCombos(tissues, ignore.full = FALSE)
  my_mat.queue <- new.queue()
  
  # BEGIN: init with flat model
  des.mat.flat <- GetFlatModel(dat.gene)
  
  if (sparse){
    # make sparse using Matrix package
    des.mat.flat <- Matrix(des.mat.flat)
  }
  
  # get rhythmic parameters which will be used for adding later: hash structure has fast lookup
  des.mat.sinhash <- GetSinCombos(dat.gene, w, tissues, tiss.combos)
  des.mat.coshash <- GetCosCombos(dat.gene, w, tissues, tiss.combos)
  
  rhyth.tiss <- list(character(0))  # needs to track shared and independent parameters, e.g.: c("Liver,Kidney", "Adr") no duplicates allowed
  # n.rhyth <- NRhythmicFromString(rhyth.tiss)  # number of independent rhythmic parameters perhaps? Do this later for speed?
  complement <- FilterCombos(tiss.combos, rhyth.tiss)  # given current matrix, you will know which tissues to iterate
  des.mat.list <- list(mat=des.mat.flat, rhyth.tiss=rhyth.tiss, complement = complement)
  # END: init with flat model
  
  # load up my queue
  enqueue(my_mat.queue, des.mat.list)
  
  # uncomment if you want to store all the matrices were used
  # des.mats <- expandingList() 
  # des.mats$add(des.mat.list)
  
  # store my fits
  # to save memory, only store top n fits
  # fits <- expandingList()
  fits <- list()
  fit.count <- 0
  
  # track weights to normalize afterwards
  fit.weights.sum <- 0
  if (!is.null(top.n)){
    fit.weights <- vector(mode="numeric", length=top.n)
  }
  
  # need to track models that we have done, so we eliminate "permutations" like c("Liver", "Kidney") and c("Kidney", "Liver) models
  # use hash for speed
  models.done <- hash()
  
  # generate matrix by adding combinations of columns and adding
  # those matrices into the queue
  while (! is.empty(my_mat.queue)) {
    des.mat.list <- dequeue(my_mat.queue)
    
    # fit my model
    if (sparse){
      fit <- FitModel(dat.gene, as.matrix(des.mat.list$mat), get.criterion = criterion, condensed = TRUE)  # condensed will save some memory
    } else {
      fit <- FitModel(dat.gene, des.mat.list$mat, get.criterion = criterion, condensed = TRUE)  # condensed will save some memory
    }
    fit.count <- fit.count + 1
    fit.weights.sum <- fit.weights.sum + fit$weight  # track even if it is a bad model, helps normalization
    if (!is.null(top.n)){
      fit.weights.lst <- UpdateFitWeights(fit$weight, fit.weights)
      if (!is.na(fit.weights.lst$i)){
        fits[[fit.weights.lst$i]] <- fit
        fit.weights <- fit.weights.lst$weights
      }
    } else {
      fits[[fit.count]] <- fit
    }
    
    
    # check that this matrix is not already the maximum complexity (n.rhyth.max),
    # if it is already as complex as we want, then ignore it because 
    # we dont want to add another rhythmic column to this.
    #
    # strictly greater than should handle the case of "no rhytmic tissues"
    if (length(des.mat.list$rhyth.tiss) > n.rhyth.max){
      # print(paste("Skipping", des.mat.list$rhyth.tiss))
      next  # should work even in n.rhyth.max == length(tissues)
    }
    
    # determine tissue combinations that need to be added based on rhyth.tiss
    # e.g., no need to add Liver twice, they can't have two rhythmic paramters  
    for (tiss.comb in des.mat.list$complement){
      # add column for each tissue combination
      tiss.key <- paste(tiss.comb, collapse = ",")
      
      # append tiss.key to rhyth.tiss
      rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)  # form list("Adr,Kidney", "Mus")
      
      # check if this tissue combination has been already submitted into queue (but in different permutation)
      # track models we have done globally
      modelname <- MakeModelName(rhyth.tiss)
      if (! is.null(models.done[[modelname]])){
        # this is a permutation of an already done combo, skip
        #       print(rhyth.tiss)
        #       print(paste('Skipping', modelname))
        next
      }
      
      col.new <- AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, tiss.key)
      
      #     rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)

      # further remove complement after having
      tiss.complement.new <- FilterCombos(des.mat.list$complement, tiss.comb)
      
      # make new matrix, put it into queue
      mat.new <- cbind(des.mat.list$mat, col.new)
      des.mat.list.new <- list(mat=mat.new, rhyth.tiss = rhyth.tiss, complement = tiss.complement.new)
      enqueue(my_mat.queue, des.mat.list.new) 
      models.done[[modelname]] <- TRUE  # we dont want to redo permutations of same models
      # des.mats$add(des.mat.list.new)
    }  
  }   
  # unpack my fits 
  # fits.list <- fits$as.list()
  # print(paste("Models fitted:", fit.count))
  fits.list <- fits
  
  if (normalize.weights){
    # print("Normalizing weights...")
    # print(dat.gene$gene[[1]])
    # fits.list <- NormalizeWeights(fits.list, cutoff)
    for (i in seq(length(fits.list))){
      fits.list[[i]]$weight.norm <- fits.list[[i]]$weight / fit.weights.sum
    }
  }
  fits.list$gene <- gene
  # print(fit.weights.sum)
  return(fits.list)
}

UpdateFitWeights <- function(weight, weights){
  # check min weights, if weight is larger than min(weights), replace min(weights) with weight
  which.min.i <- which.min(weights)
  min.weight <- weights[which.min.i]
  if (weight <= min.weight){
    return(list(weights = NA, i = NA))  # do nothing
  }
  else if (weight > min.weight){
    # replace it
    weights[which.min(weights)] <- weight
    return(list(weights = weights, i = which.min.i))
  }
}

NormalizeWeights <- function(fit.list, cutoff = 1e-3){
  # take list with $weight and 
  # add $weight.norm
  weight.sum <- sum(unlist(lapply(fit.list, function(fit){
    fit$weight
  })))
#   fit.list.norm <- Filter(function(fit){
#     fit$weight / weight.sum >= cutoff
#   }, fit.list)

  fit.list.norm <- lapply(fit.list, function(fit){
    fit$weight.norm <- fit$weight / weight.sum
    if (fit$weight.norm < cutoff){
      return(NULL)
    } else {
      return(list(fit = fit$fit, weight.norm = fit$weight.norm))
    }
  })
  return(fit.list.norm)
}

CoefToModelName <- function(coef){
  # Given set of coef, with sin and cos to designate rhythmic parameters, extract
  # a rhythmic model name e.g., Adr;Liver,Kidney;Mus means 3 rhythmic parameters, liver and kidney share 
  # rhythmic params
  rhyth.names <- names(coef[grepl(":sin", names(coef))])
  rhyth.names <- sapply(rhyth.names, function(jname) strsplit(jname, ":")[[1]][[1]])
  return(paste(rhyth.names, collapse=";"))
}

CoefToParams <- function(coef, period = 24){
  w = 2 * pi / 24
  flat.params <- coef[grepl("tissue", names(coef))]
  
  # Convert linear coefficients sin cos to amplitude and phase.
  rhyth.params.sin <- coef[grepl(":sin", names(coef))]
  rhyth.params.cos <- coef[grepl(":cos", names(coef))]
  rhyth.params.names <- names(coef)[grepl(":sin", names(coef))]
  
  if (length(rhyth.params.names) > 0){
    has.rhyth <- TRUE
  } else {
    has.rhyth <- FALSE
  }
  # TODO: check for case with no rhythmic param
  
  if (has.rhyth){
    # Liver:sin(w * time) -> Liver (can be Adr,Liver)
    rhyth.tiss <- sapply(rhyth.params.names, function(rhyth.param) strsplit(rhyth.param, ":")[[1]][[1]])
    
    # create amplitude and amp for each rhyth.tiss
    amps <- mapply(function(a, b) return(sqrt(a^2 + b^2)), rhyth.params.sin, rhyth.params.cos)
    phases <- mapply(function(a, b){
      phase.rad <- atan2(a, b)
      return((phase.rad / w) %% period)
    }, rhyth.params.sin, rhyth.params.cos)
    
    # name my new params
    if (length(amps) == 0 & length(phases) == 0) return(NA)
    names(amps) <- paste0(rhyth.tiss, ".amp")
    names(phases) <- paste0(rhyth.tiss, ".phase")
  } else {
    return(flat.params)
  }
  # return as numeric vector
  rhyth.params.out <- c(amps, phases)
  out.params <- c(flat.params, rhyth.params.out)
  return(out.params)
}

BICFromLmFit <- function(coefficients, residuals){
  # from vector of coefs and residuals, calculate BIC
  n <- length(residuals)
  RSS <- sum(residuals ^ 2)
  k <- length(coefficients)
  criterion <- n * log(RSS / n) + k * log(n)
  return(criterion)
}

AICFromLmFit <- function(coefficients, residuals){
  # from vector of coefs and residuals, calculate AIC
  n <- length(residuals)
  RSS <- sum(residuals ^ 2)
  k <- length(coefficients)
  criterion <- n * log(RSS / n) + 2 * k
  return(criterion)
}

FitModels <- function(dat.gene, my_mat, get.criterion = "BIC", normalize.weights = TRUE){
  # Fit many models with lm.fit() which is faster than lm()
  fits <- lapply(my_mat, function(mat) FitModel(dat.gene, mat, get.criterion))
  # calculate weights.sum
  weights.sum <- sum(unlist(lapply(fits, function(fit){
    return(fit$weight)
  })))
  
  # normalize weights so sum = 1
  if (normalize.weights){
    fits <- lapply(fits, function(fit){
      fit$weight.norm <- fit$weight / weight.sum
      return(fit)
    })
  }
  return(fits)
}

FitModel <- function(dat.gene, mat, weight.sum, get.criterion="BIC", condensed=FALSE){
  # Subroutine for fitting many models because
  # only one matrix, cannot normalize weights
  # 
  # weight.sum: track weight.sum 
  #   print(head(mat))
  #   print(dat.gene$exprs)
  fit <- lm.fit(y = dat.gene$exprs, x = mat)
  if (get.criterion == "BIC"){
    criterion <- BICFromLmFit(fit$coefficients, fit$residuals)
    weight <- exp(-0.5 * criterion)
    weight.raw <- criterion
  } else if (get.criterion == "AIC"){
    criterion <- AICFromLmFit(fit$coefficients, fit$residuals)
    weight <- exp(-0.5 * criterion)
    weight.raw <- criterion
  } else {
    rsquared <- GetRSquaredFromFits(dat.gene$exprs, fit$residuals)
    #     # sanity check
    #         rsquared.check <- summary(lm(dat.gene$exprs ~ mat))$r.squared
    #         print(paste("RSquared calculated: ", rsquared))
    #         print(paste("RSquared real: ", rsquared))
    N <- length(dat.gene$exprs)
    p <- length(fit$coefficients)
    bf <- tryCatch({
      bf <- GetBayesFactor(N, p, rsquared, method = get.criterion, plot.integrand = FALSE)  # output is log
    }, error = function(e) {
      print(as.character(dat.gene$gene[[1]]))
      print(e)
      bf <- NA
    })
    weight <- exp(bf)
    weight.raw <- -bf
  }
  if (condensed){
    return(list(fit = fit$coefficients, weight = weight, weight.raw = weight.raw, method = get.criterion))
  } else {
    return(list(fit = fit$coefficients, residuals = fit$residuals, weight = weight, method = get.criterion))
  }
}

GetBayesFactor <- function(N, p, rsquared, method = "zf", plot.integrand=FALSE){
  # Bayes Factor a la Liang et al 2008 mixture of g priors
  # http://www.tandfonline.com/doi/abs/10.1080/00273171.2012.734737
  # http://amstat.tandfonline.com/doi/abs/10.1198/016214507000001337#.V017rTdWeMM  
  # 
  # Taken from Richard Morey Bayes Factor package
  # github code: https://github.com/richarddmorey/BayesFactor/blob/a90b5929a8f44de1ef853b49de1d49b3c8095285/pkg/BayesFactor/R/regressionBF-utility.R
  #
  # Methods:
  # zf = "Zellner-Siow priors"
  # hyperg = "hyper-g priors"
  # zf_laplace = "Zellner-Siow priors approximated with Laplace"
  # hyperg_laplace = "hyper g approximated with Laplace
  # eb = "Empirical Bayes" use mode of g 
  
  if (method == "zf"){
    # get Marginal Likelihood of the model
    debug <- FALSE
    if (debug){
      g.mode <- ModeG(N, p, rsquared)
      gvec <- seq(0, 500000, length.out = 1000)
      bf.real <- exp(linearReg.R2stat(N=N, p=p, R2=rsquared, rscale = 1)[['bf']])  # needs BayesFactor package
      gvec.log <- seq(-50, 50, length.out = 1000)
      plot(gvec, ModelLikelihood(gvec, N=N, p=p, R2=rsquared, .log=TRUE, log.const=0, return.log=FALSE), type = "l", main=paste("log(g) vs transformed integrand. Mode:", signif(log(g.mode), digits = 2)))
      plot(gvec.log, ModelLikelihood2(gvec.log, N=N, p=p, R2=rsquared, shift = log(g.mode)), type = "l", main=paste("log(g) vs transformed integrand. Mode:", signif(log(g.mode), digits = 2)))
      # test my modellikelihood
      bf.old <- ModelLikelihood(g.mode, N=N, p=p, R2=rsquared, .log = TRUE, log.const = 0, return.log = FALSE)
      bf.bylog <- integrate(ModelLikelihood2, lower=-Inf, upper=Inf, N=N, p=p, R2=rsquared, log.const = 0, shift = log(g.mode))
      print(paste("BF old:", bf.old[[1]]))
      print(paste("BF real:", bf.real))
      print(paste("BF by log:", bf.bylog[[1]]))
    }
    # bf <- integrate(ModelLikelihood, lower=-Inf, upper=Inf, N=N, p=p, R2=rsquared, .log = TRUE, log.const = 0, return.log = FALSE)  # ModelLikelihood ratio with NULL model
    g.mode <- ModeG(N, p, R2 = rsquared)
    log.const <- ModelLikelihood2(tau = 0, N = N, p = p, R2 = rsquared, shift = log(g.mode), return.log=TRUE)
    # bf <- integrate(ModelLikelihood, lower=0, upper=Inf, N=N, p=p, R2=rsquared, .log = TRUE, log.const = 0, return.log = FALSE)  # ModelLikelihood ratio with NULL model
    bf <- integrate(ModelLikelihood2, lower=-Inf, upper=Inf, N=N, p=p, R2=rsquared, log.const = log.const, shift = log(g.mode))  # shift to centre
    if (bf[[4]] == "OK"){
      bf <- log(bf[[1]]) + log.const
    } else {
      warning(paste("Integration returned non-OK message", bf[[4]]))
      bf <- log(bf[[1]])
    }
  } else if (method == "zf_laplace"){
    warning("zf_laplace still buggy")
    # Approximate Marginal Likelihood of the model
    approx <- LaplaceApprox(N, p, rsquared)
    print(paste("mode:", approx$g.mode))
    bf <- approx$bf
  } else if (method == "hyperg"){
    a <- 3  # 2 to 4 is OK. 3 is recommended by Liang et al 2008
    bf <- log(((a - 2) / (p + a - 2)) * Re(hypergeo::hypergeo(A = (N - 1 / 2), B = 1, C = (p + a) / 2, z = rsquared)))
  } else if (method == "eb"){
    # proper way is to do EM algorithm, which I will not do.
    g.mode <- 1000
    bf <- log(ModelLikelihood(g = g.mode, N = N, p = p, R2 = rsquared, .log = TRUE, log.const = 0, return.log = FALSE))
  } else {
    # guess that user set method as "g"
    if (!is.numeric(method)){
    	# if not numeric, expect form g=1000
    	g <- as.numeric(strsplit(method, "=")[[1]][[2]])
    } else {
    	g <- method
    }
    # no prior needed (Equation 8 of Liang 2008)
    bf <- log(ModelLikelihood(g = g, N = N, p = p, R2 = rsquared, .log = TRUE, log.const = 0, return.log = FALSE, add.prior = FALSE))
  }
  if (plot.integrand){
    par(mar = c(5,5,2,5))
    g <- seq(0, 1000)
    cat(paste("N:", N, "\np:", p, "\nR2:", rsquared, "\n"))
    plot(g, ModelLikelihood(g, N = N, p = p, R2 = rsquared), type = "l", main = paste("N:", N, "\np:", p, "\nR2:", rsquared, "\n"))
    if (method == "zf_laplace"){
      abline(v = approx$g.mode, lty = 3)
      par(new = T)
      plot(g, SecondDerivG(g, N, p, rsquared), type = "l", lty = 2, axes=F, xlab=NA, ylab=NA)
      axis(side = 4)
      mtext(side = 4, line = 3, 'd2h')
    }
  }
  return(bf)
}

LaplaceApprox <- function(N, p, R2){
  # Approximate the marginal likelihood using LaPlace
  # Appendix A of Liang et al 2008
  # 
  ## Compute second derivative
  g.mode <- ModeG(N, p, R2)
  d2h <- SecondDerivG(g.mode, N, p, R2)  # Equation A.2
  print(paste("d2h:", d2h))
  sig <- (-1 * d2h) ^ -0.5
  print(paste("sig:", sig))
  bf <- sqrt(2 * pi)* sig * ModelLikelihood(g.mode, N, p, R2, .log = TRUE, log.const = 0, return.log = TRUE)  # Equation A.1
  return(list(bf = bf, g.mode = g.mode))
}

ModeG <- function(N, p, R2){
  # approximate with Laplace integral
  ### Compute approximation to posterior mode of g
  ### Liang et al Eq. A.3, assuming a=b=0
  g3 = -(1 - R2) * (p + 3) #* g^3 
  g2 = (N - p - 4 - 2 * (1 - R2)) #* g^2
  g1 = (N * (2 - R2) - 3) #*g
  g0 = N
  sol = polyroot(c(g0, g1, g2, g3))  # solutions to cubic equation
  ## Pick the real solution
  g.mode = Re(sol[which.min(Im(sol)^2)])
  if (g.mode <= 0){
    g.mode <- N / 20  # somewhere close to 0, can't do log of negs
  }
  return(g.mode)
}

SecondDerivG <- function(g, N, p, R2, a = 0, b = 0, verbose=FALSE){
  # Liang et al Eq. A. 3
  t1 <- (N - 1) * (1 - R2) ^ 2 / (1 + (g * (1 - R2))) ^ 2
  t2 <- -1 * (N - p - 1) / ((1 + g) ^ 2)
  t3 <- (3 - 2 * a) / (g ^ 2)
  t4 <- -1 * (2 * N / g ^ 3)
  d2h <- 0.5 * (t1 + t2 + t3 + t4)
  if (verbose){
    cat(paste("t1: ", t1, "\nt2", t2, "\nt3", t3, "\nt4", t4, "\n"))
  }
#   d2h <- 0.5 * ( 
#     (((N - 1)*(1 - R2)) / ((1 + g * (1 - R2)) ^ 2)) - 
#     ((N - p - 1) / ((1 + g) ^ 2)) + 
#     ((3 - 2 * a) / (g ^ 2)) - 
#     ((2 * N) / (g ^ 3))
#     )
  return(d2h)
}

ModelLikelihood <- function(g, N, p, R2, .log=FALSE, log.const=0, return.log = FALSE, method = "zf", add.prior=TRUE){
  if (!.log){
    if (add.prior){
      L <- (((1 + g) ^ ((N - p - 1) / 2)) * (1 + (1 - R2) * g) ^ -((N - 1) / 2)) * PriorG(N, g, .log = FALSE, method = "zf")
    } else {
      L <- (((1 + g) ^ ((N - p - 1) / 2)) * (1 + (1 - R2) * g) ^ -((N - 1) / 2))
    }
  } else {
    debug <- FALSE
    if (debug){
      # print((N - p -1) * log1pExp(g))
      # print(g + log(1 - R2))
      # print(.5 * ((N - p - 1 ) * log1pExp(g) - (N - 1) * log1pExp(g + log(1 - R2))))
      # a2 <- .5 * ((N - p - 1 ) * log1pExp(g) - (N - 1) * log1pExp(g + log(1 - R2)))  
      # a1 <- (((N - p - 1) / 2) * log(1 + g)) + (-(N - 1) / 2) * log(1 + (1 - R2) * g)
      a2 <- (- (N - 1)/2 * log1pExp(g + log(1 - R2)))
      a1 <-  (-(N - 1) / 2) * log(1 + (1 - R2) * g)
      print(paste("a from jake:", a1))
      print(paste("a from morey", a2))
      print(paste("a from jake + correction", a1 + log(g)))
    }
    # two lines should be equal
    # L <- (((N - p - 1) / 2) * log1pExp(log(g)) + (-(N - 1) / 2) * log1pExp(log((1 - R2) * g)) + PriorG(N, g, .log = TRUE, method = "zf"))
    if (add.prior){
      L <- (((N - p - 1) / 2) * log(1 + g)) + (-(N - 1) / 2) * log(1 + (1 - R2) * g) + PriorG(N, g, .log = TRUE, method = "zf")
    } else {
      L <- (((N - p - 1) / 2) * log(1 + g)) + (-(N - 1) / 2) * log(1 + (1 - R2) * g)
    }
    # exponentiate at the end
    if (!return.log){
      L <- exp(L)
    }
  }
  return(L)
}

ModelLikelihood2 <- function(tau, N, p, R2, log.const=0, shift=0, return.log=FALSE){
  # integrate from -Inf to Inf, function of log(g)
  # allow shift to move tau where peak is centre
  tau <- tau + shift
  a <- ((N - 1 - p) / 2) * log1pExp(tau) +  # log1pExp works when tau is large
   (-(N - 1) / 2) * log1pExp(tau + log(1 - R2))  # log1pExp(tau + log(1 - R2)) == log(1 + (1 - R2) * exp(tau))
  # a = .5 * ((N - p - 1 ) * log1pExp(tau) - 
  #             (N - 1) * log1pExp(tau + log(1 - R2)))
  L <- a + PriorG(N, tau, .log=TRUE, .logx = TRUE, method="zf") - log.const + tau
  if (!return.log){
    return(exp(L))
  } else {
    return(L)
  }
}

log1pExp <- Vectorize(function(x){
  # return log(1 + exp(x)), preventing Infs
  if (x > -log(.Machine$double.eps)){
    # log(1 + exp(x)) == x, x too large for computer
    return(x)
  } else{
    return(log(1 + exp(x)))
  }
}, "x")



PriorG <- function(N, g, .log=TRUE, .logx=FALSE, method = "zf"){
  # Evaluate prior distribution of g
  # uses Inv-Gamma(1/2, N/2)
  if (method == "zf"){
    shape <- 1/2
    scale <- N/2  # differs from BayesFactor implementation
    if (!.logx){
      prior.g <- InvGamma(g, shape, scale, .log = TRUE, .logx=.logx)
    } else {
      prior.g <- InvGamma(g, shape, scale, .log = TRUE, .logx=.logx)
    }
    # prior.g <- MCMCpack::dinvgamma(g, shape, scale)
  } else if (method == "hg"){
    a <- 3  # range from 2 to 4 is reasonable. Equation 16 of Liang et al 2008
    prior.g <- (a - 2)
  } else {
    warning("Unknown method")
  }
  return(prior.g)
}

InvGamma <- function(g, shape, scale, .log=FALSE, .logx=FALSE){
  const <- ((scale ^ shape) / gamma(shape))
  if (!.log){
    if (.logx) warning("Log x not coded for when .log=FALSE")
    d <- const * g ^ (-shape - 1) * exp(-scale / g)
  } else {
    if (!.logx){
      d <- log(const) + (-shape - 1) * log(g) + (-scale / g)
    } else {
      # g is tau
      d <- log(const) + (-shape - 1) * g + (-scale * exp(-g))
    }
  }
  return(d)
}

GetRSquaredFromFits <- function(y, residuals){
  TSS <- sum((y - mean(y)) ^ 2)
  RSS <- sum(residuals^2)
  return(1 - (RSS / TSS))
}

GetSelectionCriterion <- function(fits, model.selection = "BIC"){
  # given list of fits, calculate BIC weights
  # get coefficients and residuals, get model selection by BIC or AIC
  if (model.selection == "BIC"){
    # https://en.wikipedia.org/wiki/Bayesian_information_criterion#Definition
    # BIC = -2 log(RSS/n) + k * log(n)
    bics <- unlist(lapply(fits, function(fit){
      n <- length(fit$residuals)
      RSS <- sum(fit$residuals ^ 2)
      k <- length(fit$coefficients)
      criterion <- -2 * log(RSS / n) + k * log(n)
      bicweight <- exp(-0.5 * criterion)
      }))
    # normalize
    bics.weight <- bics / sum(bics)
  } else {
    warning("Only BIC is currently implemented.")
  }
}

MakeRhythmicDesignMatrices <- function(dat.gene, w = 2 * pi / 24, simplify=FALSE){
  # dat.gene: long format of gene expression, time and
  # conditions. Can include additional factors such as 
  # experiment.
  # simplify: return only design matrix rather than full matrix containing meta data (FALSE is debug mode)
  
  tissues <- unique(as.character(dat.gene$tissue))
  
  tiss.combos <- GetAllCombos(tissues, ignore.full = FALSE)
  my_mat.queue <- new.queue()
  
  # BEGIN: init with flat model
  des.mat.flat <- GetFlatModel(dat.gene)
  # get rhythmic parameters which will be used for adding later: hash structure has fast lookup
  des.mat.sinhash <- GetSinCombos(dat.gene, w, tissues, tiss.combos)
  des.mat.coshash <- GetCosCombos(dat.gene, w, tissues, tiss.combos)
  
  rhyth.tiss <- list(character(0))  # needs to track shared and independent parameters, e.g.: c("Liver,Kidney", "Adr") no duplicates allowed
  # n.rhyth <- NRhythmicFromString(rhyth.tiss)  # number of independent rhythmic parameters perhaps? 
  n.rhyth <- NRhythmicFromVector(rhyth.tiss)  # do that later? naw faster if we do it now
  complement <- FilterCombos(tiss.combos, rhyth.tiss)
  des.mat.list <- list(mat=des.mat.flat, rhyth.tiss=rhyth.tiss, n.rhyth=n.rhyth, complement = complement)
  # END: init with flat model
  
  # load up my queue
  enqueue(my_mat.queue, des.mat.list)
  
  n.mat.submitted <- 1
  des.mats <- expandingList() 
  des.mats$add(des.mat.list)
  
  # need to track models that we have done, so we eliminate "permutations" like c("Liver", "Kidney") and c("Kidney", "Liver) models
  # use hash for speed
  models.done <- hash()
  
  # generate matrix by adding combinations of columns and adding
  # those matrices into the queue
  while (! is.empty(my_mat.queue)) {
    des.mat.list <- dequeue(my_mat.queue)
    # determine tissue combinations that need to be added based on rhyth.tiss
    # e.g., no need to add Liver twice, they can't have two rhythmic paramters
    
    for (tiss.comb in des.mat.list$complement){
      # add column for each tissue combination
      tiss.key <- paste(tiss.comb, collapse = ",")
      
      # append tiss.key to rhyth.tiss
      rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)  # form list("Adr,Kidney", "Mus")
      
      # check if this tissue combination has been already submitted into queue (but in different permutation)
      # track models we have done globally
      modelname <- MakeModelName(rhyth.tiss)
      if (! is.null(models.done[[modelname]])){
        # this is a permutation of an already done combo, skip
        #       print(rhyth.tiss)
        #       print(paste('Skipping', modelname))
        next
      }
      
      col.new <- AddRhythmicColumns(des.mat.sinhash, des.mat.coshash, tiss.key)
      
      #     rhyth.tiss <- c(des.mat.list$rhyth.tiss, tiss.key)
      
      
      # further remove complement after having
      tiss.complement.new <- FilterCombos(des.mat.list$complement, tiss.comb)
      
      # add meta data: makes finding models easier
      n.rhyth <- des.mat.list$n.rhyth + length(tiss.comb)
      
      # make new matrix, put it into queue
      mat.new <- cbind(des.mat.list$mat, col.new)
      des.mat.list.new <- list(mat=mat.new, rhyth.tiss = rhyth.tiss, n.rhyth=n.rhyth, complement = tiss.complement.new)
      enqueue(my_mat.queue, des.mat.list.new) 
      models.done[[modelname]] <- TRUE  # we dont want to redo permutations of same models
      n.mat.submitted <- n.mat.submitted + 1
      des.mats$add(des.mat.list.new)
    }  
  } 
  print(paste("Number of matrices generated:", n.mat.submitted))
  des.mats.list <- des.mats$as.list()
  if (simplify){
    des.mats.list <- lapply(des.mats.list, function(x) return(x$mat))
  }
  return(des.mats.list)
}

MakeModelName <- function(rhyth.tiss.lst, delim = ";"){
  # From a list of rhythmic tissues, generate a hash key to track models
  # therefore: Adr,Kidney;Liver means Adr,Kidney same param, Liver independent param
  lst <- sort(unlist(rhyth.tiss.lst))
  return(paste(lst, collapse = delim))
}

AddRhythmicColumns <- function(des.mat.sinhash, des.mat.coshash, tiss.key){
  # cbind rhythnic columns, rename to track which tissues share this parameter
  col.new <- cbind(des.mat.sinhash[[tiss.key]], des.mat.coshash[[tiss.key]])
  col.new.namebase <- paste(tiss.key, sep = ",")
  colnames(col.new) <- c(paste0(col.new.namebase, ":sin(w * time)"), paste0(col.new.namebase, ":cos(w * time)"))
  return(col.new)
}

GetSinCombos <- function(dat.gene, w, tissues, combos){
  des.mat.sin <- GetRhythModel.sin(dat.gene, w)
  des.mat.sin.hash <- MakeHash(des.mat.sin, dat.gene, tissues, combos)
  return(des.mat.sin.hash)
}

GetCosCombos <- function(dat.gene, w, tissues, combos){
  des.mat.cos <- GetRhythModel.cos(dat.gene, w)
  des.mat.cos.hash <- MakeHash(des.mat.cos, dat.gene, tissues, combos)
  return(des.mat.cos.hash)
}

NRhythmicFromString <- function(rhyth.tiss, jsep = ","){
  # How many tissues are rhythmic given a comma separated stringS
  # "" suggests 0
  if (length(rhyth.tiss) == 0) return(0)
  if (rhyth.tiss == ""){
    return(0) 
  } else {
    return(length(strsplit(rhyth.tiss, jsep)[[1]]))
  }
}

NRhythmicFromVector <- function(rhyth.tiss.vec, jsep = ","){
  # if c("Adr,Liver", "Kidney") output 3 rhythmic tisues
  n.rhyth <- 0
  for (tiss in rhyth.tiss.vec){
    n.rhyth <- n.rhyth + NRhythmicFromString(tiss)
  }
  return(n.rhyth)
}

FilterCombos <- function(combos, combos.sub){
  # Filter out combos given a combos.sub, filters
  # any combo that contains even one of the elements in combo.sub
  # always removes an empty set no matter what
  combos.subl <- sapply(combos, function(comb){
    if (length(comb) == 0) return(FALSE)  # empty set
    inrhyth <- TRUE  # say it is true, loop through conditions to set to false if it is a duplicate
    for (ci in comb){
      if (ci %in% combos.sub){
        inrhyth <- inrhyth * FALSE
      } 
    }
    return(as.logical(inrhyth))
  })
  return(combos[combos.subl])
}

GetRhythmicFormula <- function(exprs, time, intercepts, w = 2 * pi / 24, with.intercept = FALSE){
  # Get formula for full rhythmic model
  # tissue: genotype, or tissues
  # x: time
  # y: expression
  # experiment: column name of experiment (can be blank I think "")
  # with.intercept: if you want tissue and experiments to be your intercepts, keep it False
  
  sinterm <- paste0(tissue, " : ", "sin(", w, " * ", time, ")")
  costerm <- paste0(tissue, " : ", "cos(", w, " * ", time, ")")
  rhyth.terms <- paste(sinterm, costerm, sep = " + ")
  
  if (with.intercept){
    intercept <- "1"
  } else {
    intercept <- "0"
  }
  
  intercepts <- paste(intercepts, collapse = "+")
  
  form <- as.formula(paste0(exprs, "~", paste(intercept, tissue, experiment, rhyth.terms, sep = "+")))
  return(form)
}

GetFlatModel <- function(dat.gene){
  des.mat.flat <- tryCatch({
    model.matrix(exprs ~ 0 + tissue + tissue:experiment, dat.gene)
  }, error = function(e) {
    # error, try without experiment
    model.matrix(exprs ~ 0 + tissue, dat.gene)
  })
}

GetRhythModel <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:sin(w * time) + tissue:cos(w * time), dat.gene)
}

MakeHash <- function(des.mat.rhyth, dat.gene, tissues, combos){
  # return as hash object requires hash library
  if (missing(tissues)) tissues <- unique(as.character(dat.gene$tissue))
  # init with single tissues
  des.mat.hash <- hash()
  for (tiss.i in seq(tissues)){
    tiss <- tissues[tiss.i]
    des.mat.hash[[tiss]] <- des.mat.rhyth[, colnames(des.mat.rhyth)[grepl(tiss, colnames(des.mat.rhyth))]]
  }
  # add combos only if they are composed of >1 tissues
  for (comb in combos){
    if (length(comb) <= 1) next
    # init
    vec <- numeric(length = nrow(des.mat.rhyth))
    for (tiss in comb){
      vec <- vec + des.mat.hash[[tiss]]
    }
    # put combo into hash as comma separated string
    comb.key <- paste(comb, collapse = ",")
    des.mat.hash[[comb.key]] <- vec
  }
  return(des.mat.hash)
}


GetRhythModel.sin <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:sin(w * time), dat.gene)
  return(des.mat.rhyth)
}

GetRhythModel.cos <- function(dat.gene, w = 2 * pi / 24){
  des.mat.rhyth <- model.matrix(exprs ~ 0 + tissue:cos(w * time), dat.gene)
}

GetRhythmicFormula.Shared <- function(exprs, time, rhythmic.factor, intercepts, w = 2 * pi / 24, with.intercept = FALSE){
  # Get formula for full rhythmic model
  # exprs: expression col
  # time: time
  # rhythmic.factor: factor describing if rhythmic or not. 
  # with.intercept: if you want tissue and experiments to be your intercepts, keep it False
  
  # get rhythmic terms
  sinterm <- paste0(rhythmic.factor, " : ", "sin(", w, " * ", time, ")")
  costerm <- paste0(rhythmic.factor, " : ", "cos(", w, " * ", time, ")")
  rhyth.terms <- paste(sinterm, costerm, sep = " + ")
  
  if (with.intercept){
    intercept <- "1"
  } else {
    intercept <- "0"
  }
  
  intercepts <- paste(intercepts, collapse = "+")
  
  form <- as.formula(paste0(exprs, "~", paste(intercept, tissue, experiment, rhyth.terms, sep = "+")))
  return(form)
}

SubsetFullDesignMatrix <- function(des.mat, rhythmic.tissues){
  # Remake design matrix to include only "rhythmic.tissues" given the full des.mat from GetRhythmicFormula
  
  # keep intercept terms
  cols.int <- !grepl(":sin|:cos", colnames(des.mat))  # all terms without :sin or :cos are intercept terms
  
  if (length(rhythmic.tissues) == 0){
    # no rhytmic tissues
    return(des.mat[, cols.int])
  }
  
  # keep subset of rhythmic terms, depending on rhythmic tissues
  grep.str <- paste0(rhythmic.tissues, ":")  # greps sin and cos terms
  grep.str <- paste0(grep.str, collapse = "|")  # greps OR
  
  cols.rhyth <- grepl(grep.str, colnames(des.mat))
  # take intercept and subset of rhythmics as subset
  des.mat <- des.mat[, cols.int | cols.rhyth]
  return(des.mat)
}

GetDesignMatrices <- function(dat, formula){
  # Return all possible design matrices from dat
}

GetAllCombos <- function(tissues, ignore.full = TRUE){
  # get subsets of tissues possible
  # ignore.full: does not consider the full model as a combination
  tissue.combos.lst <- list()
  model.i <- 1
  
  if (ignore.full){
    max.rhyth <- length(tissues) - 1  # ignore full model
  } else {
    max.rhyth <- length(tissues)
  }
  
  for(i in 0:max.rhyth){
    tissues.combo <- combn(tissues, i)
    if (is.null(ncol(tissues.combo))){
      # only one choice, put that choice into list
      tissue.combos.lst[[model.i]] <- tissues.combo
      model.i <- model.i + 1
    }
    for(j in seq(ncol(tissues.combo))){
      tissue.combos.lst[[model.i]] <- tissues.combo[, j]
      model.i <- model.i + 1
    }
  }
  return(tissue.combos.lst)
}

BICWeight <- function(bic.vec){
  # Vector of BIC values, convert to weights "probability"
  # BIC = -2 log (BayesFactor)
  # therefore, BICw = exp(-0.5 * BIC) normalized across all BICs
  bic.vec.weight <- exp(-0.5 * bic.vec)
  bic.vec.weight.frac <- bic.vec.weight / sum(bic.vec.weight)
}

FitCombinations <- function(dat.gene, tiss.combos, N=3, n.cores=30){
  # Return dataframe of combinations, fit, BIC for each possible model (2^N)
  # tiss.combos: list of tissue combinations to run lienar model
  # N: return only a subset of BIC models (top 3 by default)
  tissues <- unique(dat.gene$tissue)
  gene <- dat.gene$gene[[1]]
  n.models <- 2 ^ length(tissues)
  
  fits.lst <- list()
  bics.lst <- vector(length = n.models)
  
  #   print(paste("Number of models to fit:", n.models))
  
  # BEGIN: FULL DESIGN MATRIX
  form <- GetRhythmicFormula("exprs", "time", c("tissue", "experiment"), with.intercept = FALSE)
  des.mat.full <- model.matrix(form, dat.gene)
  # END: FULL DESIGN MATRIX
  
  # BEGIN: INIT FITTING WITH FULL MODEL
  # init with fitting full model, then just update with new formula
  des.mat.sub <- des.mat.full  # rename to keep fit output names consistent
  # fit.full <- lm(exprs ~ 0 + des.mat.sub, data = dat.gene)
  # END: INIT FITTING
  
  # BEGIN: FIT DIFFERENT COMBOS
  
  # PARALLEL
  #   out.lst <- mclapply(tiss.combos, function(tiss.combo){
  #     des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
  #     fit.long <- dat.gene %>% do(mod = lm(exprs ~ 0 + des.mat.sub, data = .)) %>%
  #       mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","), 
  #              bic = BIC(mod[[1]]))
  #     return(fit.long)
  #   }, mc.cores = n.cores)
  
  # SERIAL
  #   out.lst <- lapply(tiss.combos, function(tiss.combo){
  #     des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
  #     fit.long <- dat.gene %>% do(mod = lm(exprs ~ 0 + des.mat.sub, data = .)) %>%
  #       mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","), 
  #              bic = BIC(mod[[1]]))
  #     return(fit.long)
  #   })
  
  out.lst <- lapply(tiss.combos, function(tiss.combo){
    des.mat.sub <- SubsetFullDesignMatrix(des.mat.full, tiss.combo)
    fit.long <- dat.gene %>% 
      do(mod = FitRhythmicDesMat(., des.mat.sub)) %>%
      mutate(rhyth.tiss = paste0(tiss.combo, collapse = ","))
    return(fit.long)
  })
  
  # END: FIT DIFFERENT COMBOS
  out.df <- do.call(rbind, out.lst)
  # out.df$bicweight <- BICWeight(out.df$bic)
  out.df$bicweight <- BICWeight(sapply(out.df$mod, function(x) x$bic))
  
  # get top 3
  out.df <- out.df[order(out.df$bicweight, decreasing = TRUE), ][1:N, ]
  out.df$gene <- gene
  # gc()  # too slow
  return(out.df)
}

FitRhythmicDesMat <- function(dat, des.mat){
  # Fit rhythmic with design matix
  mod <- lm(exprs ~ 0 + des.mat, data = dat)
  bic <- BIC(mod)
  return(list(fit = coef(mod), bic = bic))
}
