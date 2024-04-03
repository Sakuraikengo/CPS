##############################################
# Title : 3.0.GS
# Author : Kengo Sakurai
# Date : 2023-07-21
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(RAINBOWR)
library(gaston)
library(breedSimulatR)
library(stringr)
library(doParallel)
library(lme4)
library(abind)
library(qtl)
packages <- c("RAINBOWR", "gaston", "lme4", "breedSimulatR", "stringr", "abind", "qtl")
source(file = "R/1.0.function.R")

# 1.2. Save dir
dirSave <- "results/3.0.GS/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave, recursive = T)
}

cl <- makeCluster(40)
registerDoParallel(cl)
foreach(i = 1:300, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
  # i <- 1
  # 1.3. Read the breeding simulation data
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
  # Extracting the genome data
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

  # selecting the candidates
  cand <- names(sort(trueU[, 1], decreasing = T)[1:nCrosses])
  crossTable <- RandomMate(cand = cand,
                           nCrosses = nCrosses,
                           nProg = nProg,
                           generation = generation)

  # creating the next generation
  nextPop <- population$new(name = "C1 offspring",
                            inds = makeCrosses(crosses = crossTable, pop = initPop))

  ###### creating the next pop ########
  # do the GP
  for (generation in 1:nGeneration) {
    # generation <- 3
    nowPop <- nextPop
    genoMat <- nowPop$genoMat

    # calculating the frequency of desirable qtns and save
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

    resEachList <- list(genoMat = genoMat,
                        trueU = trueU,
                        generation = generation)
    saveRDS(resEachList, paste0(dirSave, i, "_resList_", generation, ".rds"))

    # extracting the gamet information
    gametArray <- ExtractGamet(nowPop = nowPop, qtn = qtn)
    cand <- names(sort(trueU[, 1], decreasing = T)[1:nCrosses])

    crossTable <- RandomMate(cand = cand,
                             nCrosses = nCrosses,
                             nProg = nProg,
                             generation = generation)

    # creating the next generation
    nextPop <- population$new(name = paste0("C", generation, " offspring"),
                              inds = makeCrosses(crosses = crossTable, pop = nowPop))
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

      # list contains, "genome data", and "true genetic value"
      resF8List <- CreateF8(selfSelected = selfSelected,
                            nowPop = nowPop,
                            nProgSelf = nProgSelf,
                            markerEffectEstimated = markerEffectTrue,
                            markerEffectTrue = markerEffectTrue,
                            nTrait = nTrait,
                            generation = generation)
      saveRDS(resF8List, file = paste0(dirSave, i, "_F", f, "_", generation, ".rds"))
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
