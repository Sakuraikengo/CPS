##############################################
# Title : 3.2.CPS
# Author : Kengo Sakurai
# Date : 2024-05-14
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

for (heritability in heritabilityVec) {
  # heritability <- heritabilityVec[1]

  cl <- makeCluster(60)
  registerDoParallel(cl)
  foreach(i = 1:nRep, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
    # i <- 1

    # 1.2. Save dir
    dirSave <- paste0("results/3.2.CPS/h2_", heritability, "/")
    if (!dir.exists(dirSave)) {
      dir.create(dirSave, recursive = T)
    }

    # 1.3. Read the breeding simulation data
    simPath <- "results/2.1.simulationSetting150/"
    initPop <- readRDS(file = paste0(simPath, "h2_", heritability, "/", i, "_initPop.rds"))
    initResList <- readRDS(file = paste0(simPath, "h2_", heritability, "/", i, "_resList_0.rds"))
    map <- read.csv(file = paste0(simPath, i, "_map.csv"),
                    row.names = 1)
    markerEffectTrueMat <- read.csv(file = paste0(simPath, i, "_markerEffects.csv"),
                                    row.names = 1)
    markerEffectTrue <- c(markerEffectTrueMat)$x
    names(markerEffectTrue) <- rownames(markerEffectTrueMat)

    markerEffectEstimatedMat <- initResList$markerEffectEstimated
    markerEffectEstimated <- markerEffectEstimatedMat[, 1]
    beta <- initResList$beta

    qtn <- names(markerEffectTrue)[markerEffectTrue != 0]
    qtnMinus <- markerEffectTrue[qtn] < 0

    # read data for visualizing the location of QTNs
    cross <- read.cross(format = "csvs",
                        genfile = paste0(simPath, i, "_genoQTL.csv"),
                        phefile = paste0(simPath, i, "_phenoQTL.csv"))

    # Calculating D2 matrix
    D2List <- CalcD2(ObjectMap = map, markerEffectMat = markerEffectEstimatedMat, nTrait = nTrait, k = f - 1)
    D2 <- D2List$D2
    D2Self <- D2List$D2Self


    ###### 2. Analysis ######
    generation <- 0
    # Extracting the genome data
    genoMat <- initResList$genoMat

    # calculating the frequency of desirable qtn
    CalcFreqQTNs(genoMat = genoMat,
                 qtn = qtn,
                 qtnMinus = qtnMinus,
                 dirSave = dirSave,
                 i = i,
                 generation = generation)

    predU <- initResList$predU
    trueU <- initResList$trueU
    saveRDS(initResList, paste0(dirSave, i, "_resList_", generation, ".rds"))

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
      predProgMean <- CalcProgMean(trueU = predU)

      ##### Calculating the variance without simulation ########
      indName <- rownames(trueU)
      comb <- combn(indName, 2)
      combInd <- combn(1:length(indName), 2)
      combName <- apply(comb, 2, function(x) {
        paste0(x[1], "_", x[2])
      })

      # extracting the gamet information
      gametArray <- ExtractGamet(nowPop = initPop, qtn = colnames(genoMat))
      dimnames(gametArray)[[2]] <- colnames(genoMat)

      # Calculating the progeny variance for each cross
      sigmaMat <- CalcVarMat(D2 = D2,
                              combInd = combInd,
                              gametArray = gametArray,
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
      # generation <- 59
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
      trueU <- c((genoMat - 1) %*% markerEffectTrue)
      trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
      rownames(trueU) <- rownames(genoMat)

      # calculating the true genetic values
      predU <- c((genoMat - 1) %*% markerEffectEstimated) + beta
      predU <- matrix(predU, nrow = nrow(genoMat), ncol = nTrait)
      rownames(predU) <- rownames(genoMat)

      # check the results file
      midFile0 <- paste0(dirSave, i, "_nowPop_", (generation + 1), ".rds")
      midFile1 <- paste0(dirSave, i, "_F8_", (generation - generation %% 2), ".rds")

      if (file.exists(midFile0) & file.exists(midFile1)) {
        nextPop <- readRDS(midFile0)
      } else {
        # save the result
        resEachList <- list(pheno = NULL,
                            genoMat = genoMat,
                            predU = predU,
                            trueU = trueU,
                            markerEffectEstimated = markerEffectEstimated,
                            beta = beta,
                            generation = generation)
        saveRDS(resEachList, paste0(dirSave, i, "_resList_", generation, ".rds"))

        # calculating the mean of genotypic value for each cross
        predProgMean <- CalcProgMean(trueU = predU)

        ##### Calculating the variance without simulation ########
        indName <- rownames(trueU)
        comb <- combn(indName, 2)
        combInd <- combn(1:length(indName), 2)
        combName <- apply(comb, 2, function(x) {
          paste0(x[1], "_", x[2])
        })

        # extracting the gamet information
        gametArray <- ExtractGamet(nowPop = nowPop, qtn = colnames(genoMat))
        dimnames(gametArray)[[2]] <- colnames(genoMat)

        # Calculating the progeny variance for each cross
        sigmaMat <- CalcVarMat(D2 = D2,
                                combInd = combInd,
                                gametArray = gametArray,
                                nTrait = nTrait)
        rownames(sigmaMat) <- rownames(predProgMean) <- combName
        sigmaMat[sigmaMat < 0] <- 0

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
          selfSelected <- SelectForSelfing(trueU = predU,
                                            gametArray = gametArray,
                                            D2Self = D2Self,
                                            nTrait = nTrait,
                                            nSelf = nSelf,
                                            intensity = intensitySelf)

          # list contains, "genome data", and "true genetic value"
          resF8List <- CreateF8(selfSelected = selfSelected,
                                nowPop = nowPop,
                                nProgSelf = nProgSelf,
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
    initFilePath <- paste0(simPath, i, "_resOriginList.rds")
    f8FilePathList <- paste0(dirSave, i, "_F8_", seq(0, nGeneration, 2), ".rds")
    f8FilePathList[1] <- initFilePath
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
}

