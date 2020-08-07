#' ## Installation instructions to install data and functions to use models developed and in Yeung 2018
#' 
#'  on github (need to download GR_2018_Primetime_Objects.RData from our remote server)

# devtools::install_github("https://github.com/jakeyeung/TissueCiradianAnalysis.git")

#'  on bitbucket (GR_2018_Primetime_Objects.RData included)

# devtools::install_bitbucket("https://jakeyeung@bitbucket.org/jakeyeung/circadianrnaseq")

#' This is a quick tutorial on how to load the output of the statistical model and visualize the results. 
#' For a tutorial on applying the model selection method to the RNA-seq data, see `tutorial_fit_conditions_to_data.md`



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

# Load this .RData object: http://upnaesrv1.epfl.ch/JakeYeung_TissuePaper_RData_Objects/GR_2018_Primetime_Objects.RData
inf <- "/home/yeung/projects/CircadianRNASeq/data/GR_2018_Primetime_Objects.RData"
suppressMessages(load(inf, v=F))
# alternatively, if downloaded from bitbucket: 
# data(GR_2018_Primetime_Objects, verbose=TRUE)  # data is large >500MB, so not installed by default. If package downloaded from https://bitbucket.org/jakeyeung/circadianrnaseq, then you can load this data directly. Otherwise downloda the .RData and load it


#' ## Plot some genes

jgene <- "Dbp"
dat.sub <- subset(dat.long, gene == jgene)

#' **Expression of Dbp across tissues and time from Hogenesch data: **


PlotGeneAcrossTissues(dat.sub) + theme_bw()  + theme(aspect.ratio = 1, legend.position = "none")


#' we do model selection to identify rhythmic parameters shared or differenet across tissues


#' Model selection output automatically groups tissues by rhythmic parameters


jgene.byparam <- "Ube2u"
PlotGeneByRhythmicParameters(fits.long, subset(dat.long, experiment == "array"),
                                  jgene.byparam, amp.filt = 0.2, jtitle=jgene.byparam, facet.rows = 1, jcex = 8,
                                  pointsize = 0)


#' Here we show a liver-specific rhythmic gene. We plot just the microarray data for visualization purposes, although the model is fit on both microarray and RNA-seq together.


#' Plot genes for WT and KO to see whether a gene is clock-controlled or clock-independent (e.g., driven by feeding rhythms)


#' *Dbp* is clock controlled, because it is flat in Bmal1 KO


#' WT vs KO DBP:


PlotGeneTissuesWTKO(subset(dat.wtko, gene == jgene), jtitle = jgene, jsize = 10, single.day = TRUE)


#' The model selection for dat.wtko object is called `fits.long.filt``, here I show that Dbp is rhythmic in liver and kidney wild type, but not KO

print(subset(fits.long.filt, gene == jgene))


#' Check out `tutorial_fit_conditions_to_data.md` to a small example of how this model is applied to RNAseq data.
#' 

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

eigens.tw <- GetEigens(s.tw, period = 24, comp = svdcomponent, add.arrow = TRUE, jsize = 12, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)


#' ## Tissue-wide SVD module


multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

#' The gene loadings show how the phase and amplitude relates to each other. For example, we find *Arntl* to oscillate in phase with *Npas2* across tissues, but antiphasic with *Dbp*.
#'
#' The tissue loadings show how the oscillations of each gene relates across tissues. Here we see in this tissue-wide module that these genes oscillate in nearly all tissues, with coherent phases. However, we find the amplitudes vary (we often see brain tissues have lower amplitudes than other tissues like liver)

eigens.tw.wtko <- GetEigens(s.tw.wtko, period = 24, add.arrow = TRUE, comp = svdcomponent, jsize = 12, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))
jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)


#' ## Tissue-wide SVD module on WT and KO data


multiplot(eigens.tw.wtko$u.plot, eigens.tw.wtko$v.plot, layout = jlayout)


#' Here we find the gene loadings of the clock-controlled genes to be very comparable to the Hogenesch data (Hogenesch data will contain both clock-controlled and clock-independent genes).
#' The addition of WT and /Bmal1/ KO data now shows that indeed these genes are clock-controlled: we see the amplitudes of these genes in kidney and liver /Bmal1/ KO data.








