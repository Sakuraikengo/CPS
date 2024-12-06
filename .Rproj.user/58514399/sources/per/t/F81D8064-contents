##############################################
# Title : 3.2.CPS
# Author : Kengo Sakurai
# Date : 2023-07-21
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(RAINBOWR)
library(gaston)
library(lme4)
library(breedSimulatR)
library(lpSolve)
library(stringr)
library(doParallel)
library(tictoc)
library(abind)
library(qtl)
packages <- c("RAINBOWR", "gaston", "lme4", "breedSimulatR", "lpSolve", "stringr", "abind", "qtl")
source(file = "R/1.0.function.R")

# 1.2. Save dir
dirSave <- "results/3.2.CPS/"

if (!dir.exists(dirSave)) {
  dir.create(dirSave, recursive = T)
}

cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i = x1:300, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
  # i <- 1
  # 1.3. Read data
  simPath <- "results/2.1.simulationSetting150/"
  initPop <- readRDS(file = paste0(simPath, i, "_initPop.rds"))
  map <- read.csv(file = paste0(simPath, i, "_map.csv"),
                  row.names = 1)
  markerEffectTrueMat <- read.csv(file = paste0(simPath, i, "_markerEffects.csv"),
                                  row.names = 1)
  markerEffectTrue <- c(markerEffectTrueMat)$x
  names(markerEffectTrue) <- rownames(markerEffectTrueMat)
  qtn <- names(markerEffectTrue)[markerEffectTrue != 0]
  qtnMinus <- markerEffectTrue[qtn] < 0

  # read data for visualizing the location of QTNs
  cross <- read.cross(format = "csvs",
                      genfile = paste0(simPath, i, "_genoQTL.csv"),
                      phefile = paste0(simPath, i, "_phenoQTL.csv"))


  # calculating the T matrix for calculating progeny variance
  tMatList <- CalcT(ObjectMap = map, qtn = qtn, k = f - 1)
  tMat1 <- tMatList$tMat1
  tMat2 <- tMatList$tMat2
  tMatSelf <- tMatList$tSelf
  rm(tMatList)
  gc()

  ###### 2. Analysis ######
  generation <- 0
  genoMat <- initPop$genoMat

  # calculating the frequency of desirable qtn
  CalcFreqQTNs(genoMat = genoMat,
               qtn = qtn,
               qtnMinus = qtnMinus,
               dirSave = dirSave,
               i = i,
               generation = generation)

  trueU <- c(genoMat %*% markerEffectTrue)
  trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
  rownames(trueU) <- rownames(genoMat)

  resEachList <- list(genoMat = genoMat,
                      trueU = trueU,
                      generation = generation)
  saveRDS(resEachList, paste0(dirSave, i, "_resList_", generation, ".rds"))
  saveRDS(resEachList, paste0(dirSave, i, "_F", f, "_", generation, ".rds"))

  # checkt the old files
  midFileList0 <- list.files(path = dirSave,
                             pattern = paste0("^", i, "_nowPop_"))
  midFileListInd <- order(as.numeric(str_sub(midFileList0,
                                             start = nchar(i) + 9,
                                             end = -5)))
  midFileList <- paste0(dirSave, midFileList0)[midFileListInd]
  if (length(midFileList) > 0) {
    # read the old data
    nextPop <- readRDS(paste0(dirSave, i, "_nowPop_1.rds"))

  } else {
    # calculating the mean of genotypic value for each cross
    predProgMean <- CalcProgMean(trueU = trueU)

    ##### Calculating the variance without simulation ########
    indName <- rownames(trueU)
    comb <- combn(indName, 2)
    combName <- apply(comb, 2, function(x) {
      paste0(x[1], "_", x[2])
    })

    # extracting the gamet information
    gametArray <- ExtractGamet(nowPop = initPop, qtn = qtn)

    sigmaMat <- CalcVarMat(comb = comb,
                           gametArray = gametArray,
                           markerEffectTrue1 = markerEffectTrue,
                           markerEffectTrue2 = markerEffectTrue,
                           qtn = qtn,
                           tMat1 = tMat1,
                           tMat2 = tMat2,
                           f = f,
                           nTrait = nTrait)

    rownames(sigmaMat) <- rownames(predProgMean) <- combName

    predProgMat <- ExpectedFnValue(sigmaMat = sigmaMat,
                                   predProgMean = predProgMean,
                                   intensityCross = intensityCross)
    rm(predProgMean, sigmaMat)
    gc()

    crossTable <- LP(predProgMat = predProgMat,
                     indName = indName,
                     comb = comb,
                     capa = capa,
                     nCrosses = nCrosses,
                     nProg = nProg,
                     generation = generation)
    rm(predProgMat)
    gc()


    # creating the next generation
    nextPop <- population$new(name = "C1 offspring",
                              inds = makeCrosses(crosses = crossTable, pop = initPop))
    saveRDS(nextPop, file = paste0(dirSave, i, "_nowPop_", (generation + 1), ".rds"))
  }
  ###### creating the next pop ########
  for (generation in 1:nGeneration) {
    # generation <- 50
    nowPop <- nextPop
    genoMat <- nowPop$genoMat

    # calculating the frequency of desirable qtn
    CalcFreqQTNs(genoMat = genoMat,
                 qtn = qtn,
                 qtnMinus = qtnMinus,
                 dirSave = dirSave,
                 i = i,
                 generation = generation)

    # calculating the true genetic values
    trueU <- c(genoMat %*% markerEffectTrue)
    trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
    rownames(trueU) <- rownames(genoMat)

    # check the results file
    midFile0 <- paste0(dirSave, i, "_nowPop_", (generation + 1), ".rds")
    midFile1 <- paste0(dirSave, i, "_F8_", (generation - generation %% 2), ".rds")

    if (file.exists(midFile0) & file.exists(midFile1)) {
      nextPop <- readRDS(midFile0)
    } else {
      resEachList <- list(trueU = trueU,
                          generation = generation)
      saveRDS(resEachList, paste0(dirSave, i, "_resList_", generation, ".rds"))

      indName <- rownames(trueU)
      comb <- combn(indName, 2)
      combName <- apply(comb, 2, function(x) {
        paste0(x[1], "_", x[2])
      })

      # calculating the mean of genotypic value for each cross
      predProgMean <- CalcProgMean(trueU = trueU)

      # extracting the gamet information
      gametArray <- ExtractGamet(nowPop = nowPop, qtn = qtn)

      sigmaMat <- CalcVarMat(comb = comb,
                             gametArray = gametArray,
                             markerEffectTrue1 = markerEffectTrue,
                             markerEffectTrue2 = markerEffectTrue,
                             qtn = qtn,
                             tMat1 = tMat1,
                             tMat2 = tMat2,
                             f = f,
                             nTrait = nTrait)

      rownames(sigmaMat) <- rownames(predProgMean) <- combName

      predProgMat <- ExpectedFnValue(sigmaMat = sigmaMat,
                                     predProgMean = predProgMean,
                                     intensityCross = intensityCross)
      rm(predProgMean, sigmaMat)
      gc()

      crossTable <- LP(predProgMat = predProgMat,
                       indName = indName,
                       comb = comb,
                       capa = capa,
                       nCrosses = nCrosses,
                       nProg = nProg,
                       generation = generation)
      rm(predProgMat)
      gc()

      # creating the next generation
      nextPop <- population$new(name = paste0("C", generation, " offspring"),
                                inds = makeCrosses(crosses = crossTable, pop = nowPop))
      saveRDS(nextPop, file = paste0(dirSave, i, "_nowPop_", (generation + 1), ".rds"))

      # selecting the genotypes for releasing varieties
      if (generation %% 2 == 0) {
        selfSelected <- SelectForSelfing1(trueU = trueU,
                                          gametArray = gametArray,
                                          markerEffectTrue1 = markerEffectTrue,
                                          qtn = qtn,
                                          tMat1 = tMat1,
                                          tMatSelf = tMatSelf,
                                          f = f - 2,
                                          nTrait = nTrait,
                                          nSelf = nSelf,
                                          intensity = intensitySelf)

        # list contains "phenotypic data", "genome data",
        # "estimated genetic value", and "true genetic value"
        resF8List <- CreateF8(selfSelected = selfSelected,
                              nowPop = nowPop,
                              nProgSelf = nProgSelf,
                              markerEffectEstimated = markerEffectTrue,
                              markerEffectTrue = markerEffectTrue,
                              nTrait = nTrait,
                              generation = generation)
        # resEachList <- c(resEachList, list(resEachListNew))
        saveRDS(resF8List, file = paste0(dirSave, i, "_F", f, "_", generation, ".rds"))
      }
    }
  }
  # summary of the recurrent genomic selection
  rgsFilePathList <- paste0(dirSave, i, "_resList_", 0:nGeneration, ".rds")
  rgsRowNames <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))
  rgsFileName <- paste0(dirSave, i, "_geneticValMat.csv")
  SummarizeResults(dataPathList = rgsFilePathList,
                   rowNamesData = rgsRowNames,
                   plotTitle = "Genetic Value of Generation ",
                   saveFileName = rgsFileName)

  # summary of the Inbred4 population
  f8FilePathList <- paste0(dirSave, i, "_F8_", seq(0, nGeneration, 2), ".rds")
  f8RowNames <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))
  f8FileName <- paste0(dirSave, i, "_resultMat.csv")
  SummarizeResults(dataPathList = f8FilePathList,
                   rowNamesData = f8RowNames,
                   plotTitle = "Genetic Value of F8 from Generation ",
                   saveFileName = f8FileName)

  # summary of the desired allele frequency
  freqFilePathList <- paste0(dirSave, i, "_desiredFreqQTN_", 0:nGeneration, ".csv")
  SummarizeFreq(freqFilePathList = freqFilePathList,
                nGeneration = nGeneration,
                i = i,
                cross = cross)
}
stopCluster(cl)
