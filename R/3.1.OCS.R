##############################################
# Title : 3.1.OCS
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

for (targetGeneration in targetGenerationVec) {
  # targetGeneration <- targetGenerationVec[2]

  for (degree in degreeVec) {
    # degree <- degreeVec[3]

    for (HeStarRatio in HeStarRatioVec) {
      # HeStarRatio <- HeStarRatioVec[3]
      # 1.2. Save dir
      dirSave <- paste0("results/3.1.OCS/", targetGeneration, "_", degree, "_", HeStarRatio, "/")

      if (!dir.exists(dirSave)) {
        dir.create(dirSave, recursive = T)
      }

      cl <- makeCluster(60)
      registerDoParallel(cl)
      foreach(i = 1:300, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
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
        #
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

        # Calculating the He0 (the initial neutral diversity)
        genoMatSel <- genoMat[, qtn]
        He0 <- CalcHe0(genoMat = genoMatSel)[1]
        HeStar <- He0 * HeStarRatio

        if (length(midFileList) > 0) {
          # read the old data
          nextPop <- readRDS(paste0(dirSave, i, "_nowPop_1.rds"))

        } else {
          # Calculating the He_t (genetic diversity required in generation t)
          if ((generation + 1) > targetGeneration) {
            He_t <- HeStar
          } else {
            He_t <- He0 + ((generation + 1) / targetGeneration)^degree * (HeStar - He0)
          }

          # Calculating K (genetic relationship matrix)
          K <- ((((genoMatSel-1) %*% t(genoMatSel-1)) / ncol(genoMatSel)) + 1) / 2

          # calculating the mean of genotypic value for each cross
          predProgMean <- CalcProgMean(trueU = trueU)

          indName <- rownames(trueU)
          comb <- combn(indName, 2)
          combName <- apply(comb, 2, function(x) {
            paste0(x[1], "_", x[2])
          })
          rownames(predProgMean) <- combName

          # A design matrix linking the N potential parents to the first (respectively second) parent in the cross list
          Z1 <- t(design.Z(pheno.labels = comb[1, ], geno.names = indName))
          Z2 <- t(design.Z(pheno.labels = comb[2, ], geno.names = indName))

          # selecting the crosses partially randomlly
          selectedInd <- SearchInitOC(predProgMean = predProgMean,
                                      comb = comb,
                                      combName = combName,
                                      K = K,
                                      Z1 = Z1,
                                      Z2 = Z2,
                                      He = He_t,
                                      initInd = initInd,
                                      maxIter = maxIter,
                                      nCrosses = nCrosses,
                                      capa = capa)
          if (length(selectedInd) != 10) {
            print("We can't find optimal crosses")
            stop()
          }

          # optimizing the selected crosses
          selectedOptimalInd <- SearchOptimalOC(selectedInd = selectedInd,
                                                predProgMean = predProgMean,
                                                comb = comb,
                                                combName = combName,
                                                K = K,
                                                Z1 = Z1,
                                                Z2 = Z2,
                                                He = He_t,
                                                nCrosses = nCrosses,
                                                capa = capa)
          rm(Z1, Z2)
          gc()
          bestCrossMat <- t(comb[, selectedOptimalInd])

          crossTable <- data.frame(ind1 = bestCrossMat[, 1],
                                   ind2 = bestCrossMat[, 2],
                                   n = nProg,
                                   names = paste0("C", generation, "N",
                                                  formatC(1:nrow(bestCrossMat),
                                                          width = 3, flag = "0")))

          # creating the next generation
          nextPop <- population$new(name = "C1 offspring",
                                    inds = makeCrosses(crosses = crossTable, pop = initPop))
          saveRDS(nextPop, file = paste0(dirSave, i, "_nowPop_", (generation + 1), ".rds"))
        }
        ###### creating the next pop ########
        # do the GP
        for (generation in 1:nGeneration) {
          # generation <- 1
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
            resEachList <- list(genoMat = genoMat,
                                trueU = trueU,
                                generation = generation)
            saveRDS(resEachList, paste0(dirSave, i, "_resList_", generation, ".rds"))
            rm(resEachList)
            gc()

            # Calculating the He_t (genetic diversity required in generation t)
            # He_t <- He0 + ((HeStar - He0) * generation) / nGeneration
            if ((generation + 1) > targetGeneration) {
              He_t <- HeStar
            } else {
              He_t <- He0 + ((generation + 1) / targetGeneration)^degree * (HeStar - He0)
            }

            # Calculating K (genetic relationship matrix)
            genoMatSel <- genoMat[, qtn]
            K <- ((((genoMatSel-1) %*% t(genoMatSel-1)) / ncol(genoMatSel)) + 1) / 2

            # calculating the mean of genotypic value for each cross
            predProgMean <- CalcProgMean(trueU = trueU)

            indName <- rownames(trueU)
            comb <- combn(indName, 2)
            combName <- apply(comb, 2, function(x) {
              paste0(x[1], "_", x[2])
            })
            rownames(predProgMean) <- combName

            # A design matrix linking the N potential parents to the first (respectively second) parent in the cross list
            Z1 <- t(design.Z(pheno.labels = comb[1, ], geno.names = indName))
            Z2 <- t(design.Z(pheno.labels = comb[2, ], geno.names = indName))

            # extracting the gamet information
            gametArray <- ExtractGamet(nowPop = nowPop, qtn = qtn)

            selectedInd <- SearchInitOC(predProgMean = predProgMean,
                                        comb = comb,
                                        combName = combName,
                                        K = K,
                                        Z1 = Z1,
                                        Z2 = Z2,
                                        He = He_t,
                                        initInd = initInd,
                                        maxIter = maxIter,
                                        nCrosses = nCrosses,
                                        capa = capa)
            if (length(selectedInd) != 10) {
              print("We can't find optimal crosses")
              stop()
            }

            selectedOptimalInd <- SearchOptimalOC(selectedInd = selectedInd,
                                                  predProgMean = predProgMean,
                                                  comb = comb,
                                                  combName = combName,
                                                  K = K,
                                                  Z1 = Z1,
                                                  Z2 = Z2,
                                                  He = He_t,
                                                  nCrosses = nCrosses,
                                                  capa = capa)
            rm(Z1, Z2)
            gc()
            bestCrossMat <- t(comb[, selectedOptimalInd])
            crossTable <- data.frame(ind1 = bestCrossMat[, 1],
                                     ind2 = bestCrossMat[, 2],
                                     n = nProg,
                                     names = paste0("C", generation, "N",
                                                    formatC(1:nrow(bestCrossMat),
                                                            width = 3, flag = "0")))

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
    }
  }
}

