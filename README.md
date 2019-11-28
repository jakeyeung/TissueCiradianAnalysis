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

Load data
---------

    inf <- "data/GR_2018_Primetime_Objects.Rdata"
    load(inf, v=T)

    ## Loading objects:
    ##   ModelStrToModel
    ##   plot.dir
    ##   NormalizeWeightsFromFits
    ##   UnzipMatrices
    ##   .sql
    ##   remove.kidney.outliers
    ##   CollapseTissueGeno
    ##   PlotComplex
    ##   GetPhi
    ##   MakeRhythmicDesignMatrices
    ##   PlotRnaMicroarrayFit
    ##   UpdateFitWeights
    ##   dat.freq
    ##   SubsetGenoSignalTimeType
    ##   LoadActivitiesLongDhs
    ##   ListToLong
    ##   dat.long
    ##   OscillateRelamp
    ##   StaggeredTimepointsLivKid
    ##   PlotHeatmapAmpPhasePval
    ##   PlotEigensamp
    ##   GetAvgAmpFromParams
    ##   .Random.seed
    ##   PlotAgainstRnaSeq
    ##   PlotEncodeRnaseq
    ##   PlotComplex2
    ##   .First
    ##   DoFourier.se
    ##   SameTimepointsLivKid
    ##   ModelLikelihood
    ##   GetSdPhaseFromParams
    ##   PlotLoadings
    ##   GetEigens
    ##   PlotComplexCircle
    ##   CountModels
    ##   ScaleSignal
    ##   BICFromLmFit
    ##   GetUnobsObsSymbol
    ##   RunPldaSystemsClock
    ##   CalculatePeriodogramLong
    ##   zscore.min
    ##   PlotPolarHistogram
    ##   ProjectWithZscore
    ##   GetFullR
    ##   ProjectToPeriodicTime
    ##   GetTissuesFromCoefFit
    ##   is_outlier
    ##   GetSinCombos
    ##   SecondDerivG
    ##   PlotOverlayTimeSeries
    ##   ExtraParamsFromFit
    ##   SvdOnComplex
    ##   FitCombinations
    ##   start
    ##   GetTopNModelsFromFits
    ##   AddDf
    ##   GetPhi.se
    ##   GetRSquaredFromFits
    ##   jmeth
    ##   GetTFs
    ##   PlotAmpPhase
    ##   RemoveP2Name
    ##   ChunkDatGenesToFile
    ##   NGenesByAmp
    ##   RemoveLowlyExpressedGenes
    ##   PlotHeatmapNconds
    ##   n
    ##   log1pExp
    ##   MakeCnamesLivKidWTKO
    ##   circular_phase24H_histogram
    ##   PlotGeneAcrossTissuesRnaseq
    ##   .N.sub.liv
    ##   GetHistCounts
    ##   PlotExprsHeatmap
    ##   DelPhi.Dela
    ##   PriorG
    ##   DelPhi.Delb
    ##   prot.long
    ##   RemoveCommasBraces
    ##   GetGenesFromMotifs
    ##   ConvertToSingleDay
    ##   MakeDesMatRunFit
    ##   PlotMeanExprsOfModel
    ##   CalculatePeriodogram
    ##   GetTopGenesFromSvd
    ##   PlotGeneByRhythmicParameters
    ##   ProjectOnVector
    ##   FitToMat
    ##   GetPhaseFromParams
    ##   FilterLiverGenes
    ##   PlotGeneAcrossTissues
    ##   PlotProteomics
    ##   NormalizeComplexMat
    ##   PlotPeriodogram
    ##   FindMaxFreqs
    ##   PlotSvdFromGeneList
    ##   multiplot
    ##   LoadFitsFromModels
    ##   ConcatenateFits
    ##   tissue.order
    ##   fits.bytiss
    ##   fits.long
    ##   .ls.objects
    ##   .Last
    ##   NGenesByAmp.long
    ##   CoefToModelName
    ##   LoadChunkRunNconds
    ##   dat.wtko.collapsed
    ##   RemoveLowExprsPseudoShortGenes
    ##   GetInterval
    ##   ModeG
    ##   DiffPhase
    ##   PlotTpmAcrossTissuesWTKO
    ##   NrhythToStr
    ##   AddColname
    ##   Transform
    ##   Center
    ##   MakeDesMatFromModelName
    ##   OrderDecreasing
    ##   TemporalToFrequency2
    ##   DelFA.DelI
    ##   remove.wfat
    ##   GetBestModel
    ##   FilterGenesByTissue
    ##   dat.wtko
    ##   dotsize
    ##   PlotComplexLong
    ##   GetMaxPhaseDiffFromParams
    ##   PlotTpmAcrossTissues
    ##   LoadActivitiesLong
    ##   SaveDatGenes
    ##   GetFourierEntropy
    ##   SubsetByMaxBicWeight
    ##   GetSelectionCriterion
    ##   dat.complex
    ##   PlotActivitiesWithSE.dhs
    ##   fits.long.filt
    ##   LoadPhosphoData
    ##   MakeModelName
    ##   PlotFitDiagnostics
    ##   PlotRnaseqAcrossTissues
    ##   AdjustPhase
    ##   SHOWTHIS
    ##   .rdataPath
    ##   ProjectToFrequency2
    ##   PcToMotif
    ##   PlotAmpPhaseAllTissues
    ##   GetRhythModel.cos
    ##   TissueFromCname
    ##   PlotActivitiesWithSE.rnaseq
    ##   PlotActivitiesWithSE.wtko
    ##   InvGamma
    ##   GetBayesFactor
    ##   GetFlatModel
    ##   is_top_N
    ##   DoFourier
    ##   IsRhythmic
    ##   LoadActivitiesLongKidneyLiver
    ##   TimeFromCname
    ##   ProjectToFrequency
    ##   SumWeights
    ##   PlotGeneNormalized
    ##   GetTimesTissuesGenoKL
    ##   AICFromLmFit
    ##   gg_color_hue
    ##   GetNrhythFromModel
    ##   ConvertArgToPhase
    ##   GetRhythmicFormula.Shared
    ##   FracVarByGeneList
    ##   FitRhythmicDesMat
    ##   MakeDesMatRunFitEnv
    ##   CoefToParams
    ##   tfs
    ##   PlotFirstNComponents
    ##   PlotDiagnostics
    ##   RemoveDashes
    ##   GetNoiseFactor
    ##   GetRhythModel.sin
    ##   make_circ_coord
    ##   PlotEigengene
    ##   GetGenoFromSamp
    ##   PlotmRNAActivityProtein
    ##   ModelLikelihood2
    ##   PlotMeanActivitiesWithSE
    ##   ggplotColours
    ##   PlotPeriodogramLong
    ##   PlotHeatmapGeneList
    ##   LaplaceApprox
    ##   PlotBeforeAfter
    ##   GetMotifFromGene
    ##   TemporalToFrequency
    ##   GetRhythmicFormula
    ##   NormalizeWeights
    ##   PlotArgsMatrix
    ##   PlotAmpPhaseTissue
    ##   PlotGeneTissuesWTKO
    ##   LoadLivKid
    ##   MakeHash
    ##   GetRhythModel
    ##   PlotGenesInModel
    ##   AtanPositive
    ##   plot.i
    ##   RemoveP2AndCommaBracesDashes
    ##   DelF.DelI
    ##   PlotComponentOriginal
    ##   PlotMeanActivitiesWithSE.singleintercept
    ##   GetAmpPhaseFromActivities
    ##   GetAllCombos
    ##   PlotComplex2.act
    ##   GetAmp.se
    ##   GetCosCombos
    ##   PlotTissuePcaFromGenelist
    ##   FitModel
    ##   ModelToTissue
    ##   SaveChunk
    ##   GetTimeFromSamp
    ##   NRhythmicFromString
    ##   LoadDesMatDatGeneRunFits
    ##   hr.shift
    ##   NRhythmicFromVector
    ##   PlotRelampHeatmap
    ##   GetAmp
    ##   PCbiplot
    ##   TemporalToFrequencyDatLong
    ##   MakeDesMatChunks
    ##   mrna.hl
    ##   MaxDiffPhaseVector
    ##   DelR.Dela
    ##   DelR.Delb
    ##   FilterModelByAmp
    ##   GetOmegas
    ##   GetTopMotifs
    ##   `.__T__$:base`
    ##   SetdiffTimepoints
    ##   FilterCombos
    ##   omega
    ##   LoadProteomicsData
    ##   GetDesignMatrices
    ##   FitModels
    ##   AddRhythmicColumns
    ##   PlotActivitiesWithSE
    ##   PlotGOTermsPhase
    ##   BICWeight
    ##   DelFB.DelI
    ##   GetSinusoid
    ##   SubsetFullDesignMatrix
    ##   outliers
    ##   FilterKidneyLiverGenes

Plot some genes
---------------

    jgene <- "Dbp"
    dat.sub <- subset(dat.long, gene == jgene)

**Expression of Dbp across tissues and time from Hogenesch data: **
`#+ fig.width=4, fig.height=4, dpi=50`

    PlotGeneAcrossTissues(dat.sub) + theme_bw()  + theme(aspect.ratio = 1, legend.position = "none")

![](README_files/figure-markdown_strict/unnamed-chunk-4-1.png)

we do model selection to identify rhythmic parameters shared or
differenet across tissues **Model selection output automatically groups
tissues by rhythmic parameters: **
`#+ fig.width=4, fig.height=4, dpi=50`

    jgene.byparam <- "Ube2u"
    PlotGeneByRhythmicParameters(fits.long, subset(dat.long, experiment == "array"),
                                      jgene.byparam, amp.filt = 0.2, jtitle=jgene.byparam, facet.rows = 1, jcex = 8,
                                      pointsize = 0)

    ## Warning: Detecting old grouped_df format, replacing `vars` attribute by
    ## `groups`

![](README_files/figure-markdown_strict/unnamed-chunk-5-1.png)

Here we show a liver-specific rhythmic gene. We plot just the microarray
data for visualization purposes, although the model is fit on both
microarray and RNA-seq together. \#\#\# Plot genes for WT and KO to see
whether a gene is clock-controlled or clock-independent (e.g., driven by
feeding rhythms)

Dbp is clock controlled, because it is flat in Bmal1 KO **WT vs KO DBP:
** `#+ fig.width=4, fig.height=4, dpi=50`

    PlotGeneTissuesWTKO(subset(dat.wtko, gene == jgene), jtitle = jgene, jsize = 10, single.day = TRUE)

    ## Warning: Detecting old grouped_df format, replacing `vars` attribute by
    ## `groups`

![](README_files/figure-markdown_strict/unnamed-chunk-6-1.png)

The model selection for dat.wtko object is called \`fits.long.filt\`\`,
here I show that Dbp is rhythmic in liver and kidney wild type, but not
KO

    print(subset(fits.long.filt, gene == jgene))

    ## Warning: Detecting old grouped_df format, replacing `vars` attribute by
    ## `groups`

    ## # A tibble: 1 x 12
    ## # Groups:   gene, method [1]
    ##   gene  model weight weight.raw param.list method n.params n.rhyth amp.avg
    ##   <fct> <fct>  <dbl>      <dbl> <list>     <chr>     <int>   <int>   <dbl>
    ## 1 Dbp   Live…  0.966      -48.7 <dbl [6]>  g=1001        1       2    2.83
    ## # … with 3 more variables: phase.sd <dbl>, phase.maxdiff <dbl>,
    ## #   phase.avg <dbl>

Complex-valued SVD
------------------

summarize tissue-wide genes using complex-valued SVD, a compact way to
visualize how modules of genes oscillate across tissues

    svdcomponent <- 1  # look at first singular values, which captures the most variance
    genes.tw <- as.character(subset(fits.long, n.rhyth >= 8)$gene)
    genes.tw.wtko <- as.character(subset(fits.long.filt, model %in% c("Liver_SV129,Kidney_SV129"))$gene)

    # on Hogenesch data
    s.tw <- SvdOnComplex(subset(dat.complex, gene %in% genes.tw), value.var = "exprs.transformed")

    ## Warning: Detecting old grouped_df format, replacing `vars` attribute by
    ## `groups`

    # on Naef data
    s.tw.wtko <- SvdOnComplex(subset(dat.freq, gene %in% genes.tw.wtko), value.var = "exprs.transformed")

    ## Warning: Detecting old grouped_df format, replacing `vars` attribute by
    ## `groups`

    eigens.tw <- GetEigens(s.tw, period = 24, comp = svdcomponent, add.arrow = TRUE, jsize = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))

    ## Warning in if (!is.na(label.gene)) {: the condition has length > 1 and only
    ## the first element will be used

    jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

**Tissue-wide SVD module: ** `#+ fig.width=4, fig.height=4, dpi=50`

    multiplot(eigens.tw$u.plot, eigens.tw$v.plot, layout = jlayout)

![](README_files/figure-markdown_strict/unnamed-chunk-9-1.png)

The gene loadings show how the phase and amplitude relates to each
other. For example, we find /Arntl/ to oscillate in phase with /Npas2/
across tissues, but antiphasic with /Dbp/.

The tissue loadings show how the oscillations of each gene relates
across tissues. Here we see in this tissue-wide module that these genes
oscillate in nearly all tissues, with coherent phases. However, we find
the amplitudes vary (we often see brain tissues have lower amplitudes
than other tissues like liver)

    eigens.tw.wtko <- GetEigens(s.tw.wtko, period = 24, add.arrow = TRUE, comp = svdcomponent, jsize = 2, label.n = 15, eigenval = TRUE, adj.mag = TRUE, constant.amp = dotsize, peak.to.trough = TRUE, label.gene = c("Dbp", "Arntl", "Per2", "Nr1d1"))

    ## Warning in if (!is.na(label.gene)) {: the condition has length > 1 and only
    ## the first element will be used

    jlayout <- matrix(c(1, 2), 1, 2, byrow = TRUE)

**Tissue-wide SVD module on WT and KO data **
`#+ fig.width=4, fig.height=4, dpi=50`

    multiplot(eigens.tw.wtko$u.plot, eigens.tw.wtko$v.plot, layout = jlayout)

![](README_files/figure-markdown_strict/unnamed-chunk-11-1.png)

Here we find the gene loadings of the clock-controlled genes to be very
comparable to the Hogenesch data (Hogenesch data will contain both
clock-controlled and clock-independent genes). The addition of WT and
/Bmal1/ KO data now shows that indeed these genes are clock-controlled:
we see the amplitudes of these genes in kidney and liver /Bmal1/ KO
data.