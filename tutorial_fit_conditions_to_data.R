# Jake Yeung
# Date of Creation: 2020-08-03
# File: ~/projects/CircadianRNASeq/test_nconds_fits.R
# Test nconds fits


library(CircadianRNASeq)
# library(Matrix)
# library(hash)
# library(numbers)  # easy bell number calculation
# library(ggplot2)


#' ## Load data
data(dat_long_WTKO, verbose=TRUE)




#' # Run fits

#' Subset data to just two genes for the tutorial

dat.long <- subset(dat_long_WTKO, gene %in% c("Nr1d1", "Insig2"))

dat.long.forplot <- dat.long %>%
  rowwise() %>%
  mutate(geno = strsplit(as.character(tissue), "_")[[1]][[2]],
         tissue = strsplit(as.character(tissue), "_")[[1]][[1]])

#' Expression of Insig2, a gene that is rhythmic in only one of the four conditions (Liver WT):


PlotGeneTissuesWTKO(subset(dat.long.forplot, gene == "Insig2"))

tissues.uniq <- unique(as.character(dat.long$tissue))
dat.env <- DatLongToEnvironment(dat.long)

method <- "g=1000"  # Method for penalizing model complexity: other possibilities: BIC, or g={integer}
jstart <- Sys.time()
print("Running fits...")
fits.all <- lapply(ls(dat.env), function(gene){
  MakeDesMatRunFitEnv(dat.env, gene, tissues.uniq,
                      n.rhyth.max = length(tissues.uniq), w = 2 * pi / 24,
                      criterion = method, normalize.weights = TRUE,
                      cutoff = 1e-5, top.n = NULL, sparse = FALSE)
})
print("Running fits... done")
print(Sys.time() - jstart)

n.combos <- numbers::bell(length(tissues.uniq) + 1)
fits.all.long <- lapply(fits.all, function(x){
  gene <- x$gene
  x$gene <- NULL
  fits.temp.long <- ListToLong(x, gene, top.n = n.combos)
})
fits.all.long <- do.call(rbind, fits.all.long)


fits.long.filt.subset <- fits.all.long %>%
  group_by(gene) %>%
  filter(weight.raw == min(weight.raw))

print(fits.long.filt.subset$param.list[[1]])

#' The model allows automatic separation between oscillatory conditions (here Insig2 oscillates in only Liver WT) versus "flat" conditions (all other conditions):

jgene <- "Insig2"


#' Model separates between oscillating conditions versus flat conditions



PlotGeneByRhythmicParameters(fits.long.filt.subset, dat_long_WTKO, jgene) + xlab("Time (ZT)")




