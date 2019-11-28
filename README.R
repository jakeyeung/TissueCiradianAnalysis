# Jake Yeung
# Date of Creation: 2019-11-27
# File: ~/projects/CircadianRNASeq/README.R
# Test things out before moving to Rmd

library(reshape2)
library(dplyr)
library(ggplot2)
library(hash)
library(CircadianRNASeq)
library(here())

setwd(here())


#' ## Load data

inf <- "data/GR_2018_Primetime_Objects.Rdata"
load(inf, v=T)



#' ## Plot some genes

jgene <- "Dbp"
dat.sub <- subset(dat.long, gene == jgene)

#' **Expression of Dbp across tissues and time from Hogenesch data: ** `#+ fig.width=4, fig.height=4, dpi=50`

PlotGeneAcrossTissues(dat.sub) + theme_bw()  + theme(aspect.ratio = 1, legend.position = "none")

#' we do model selection to identify rhythmic parameters shared or differenet across tissues

#' **Model selection output automatically groups tissues by rhythmic parameters: ** `#+ fig.width=4, fig.height=4, dpi=50`
#'
jgene.byparam <- "Ube2u"
PlotGeneByRhythmicParameters(fits.long, subset(dat.long, experiment == "array"),
                                  jgene.byparam, amp.filt = 0.2, jtitle=jgene.byparam, facet.rows = 1, jcex = 8,
                                  pointsize = 0)

#' Here we show a liver-specific rhythmic gene. We plot just the microarray data for visualization purposes, although the model is fit on both microarray and RNA-seq together.

#' ### Plot genes for WT and KO to see whether a gene is clock-controlled or clock-independent (e.g., driven by feeding rhythms)
#'
#' Dbp is clock controlled, because it is flat in Bmal1 KO

#' **WT vs KO DBP: ** `#+ fig.width=4, fig.height=4, dpi=50`

PlotGeneTissuesWTKO(subset(dat.wtko, gene == jgene), jtitle = jgene, jsize = 10, single.day = TRUE)

#' The model selection for dat.wtko object is called `fits.long.filt``, here I show that Dbp is rhythmic in liver and kidney wild type, but not KO

print(subset(fits.long.filt, gene == jgene))


#' ## Complex-valued SVD


#' summarize tissue-wide genes using complex-valued SVD, a compact way to visualize
#' how modules of genes oscillate across tissues


svdcomponent <- 1  # look at first singular values, which captures the most variance
genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
genes.tw.wtko <- as.character(subset(fits.long.filt, model %in% c("Liver_SV129,Kidney_SV129"))$gene)

# on Hogenesch data
s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")
# on Naef data
s.tw.wtko <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw.wtko), value.var = "exprs.transformed")

eigens.tw <- GetEigens(s.tw, period = 24, comp = svdcomponent, add.arrow = TRUE, jsize = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

#' **Tissue-wide SVD module: ** `#+ fig.width=4, fig.height=4, dpi=50`
multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

#' The gene loadings show how the phase and amplitude relates to each other. For example, we find /Arntl/ to oscillate in phase with /Npas2/ across tissues, but antiphasic with /Dbp/.
#'
#' The tissue loadings show how the oscillations of each gene relates across tissues. Here we see in this tissue-wide module that these genes oscillate in nearly all tissues, with coherent phases. However, we find the amplitudes vary (we often see brain tissues have lower amplitudes than other tissues like liver)

eigens.tw.wtko <- GetEigens(s.tw.wtko, period = 24, add.arrow = TRUE, comp = svdcomponent, jsize = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

#' **Tissue-wide SVD module on WT and KO data ** `#+ fig.width=4, fig.height=4, dpi=50`
multiplot(eigens.tw.wtko$u.plot, eigens.tw.wtko$v.plot, layout = jlayout)


#' Here we find the gene loadings of the clock-controlled genes to be very comparable to the Hogenesch data (Hogenesch data will contain both clock-controlled and clock-independent genes).
#' The addition of WT and /Bmal1/ KO data now shows that indeed these genes are clock-controlled: we see the amplitudes of these genes in kidney and liver /Bmal1/ KO data.









