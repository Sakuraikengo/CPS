##############################################
# Title : 2.1.simulationSetting150
# Author : Kengo Sakurai
# Date : 2024-03-14
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(breedSimulatR)
library(stringr)
library(qtl)
library(doParallel)
library(RAINBOWR)
library(MASS)
packages <- c("breedSimulatR", "stringr", "qtl", "RAINBOWR", "MASS")
source(file = "R/1.0.function.R")

# 1.1. Parameters
nCluster <- 4
nPheno <- 1

# 1.2. Save dir
dirSave <- "results/2.1.simulationSetting150/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

seedIndCsv <- paste0("results/seedInd.csv")
seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
seedInd <- c(as.matrix(seedInd))

cl <- makeCluster(60)
registerDoParallel(cl)
foreach(i = 1:nRep, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
  # i <- 1
  # 1.3. Read data
  simPath <- "results/2.0.simulationSetting198/"
  initPop <- readRDS(file = paste0(simPath, i, "_initPop.rds"))
  markerEffectTrueMat <- read.csv(file = paste0(simPath, i, "_markerEffects.csv"),
                                  row.names = 1)
  markerEffectTrue <- c(markerEffectTrueMat)$x
  names(markerEffectTrue) <- rownames(markerEffectTrueMat)
  qtn <- names(markerEffectTrue)[markerEffectTrue != 0]

  ###### 2. Analysis ######
  # determine the selected individuals
  genoMat <- initPop$genoMat

  trueU <- c((genoMat - 1) %*% markerEffectTrue)
  trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
  rownames(trueU) <- rownames(genoMat)

  # Clustering into 4 groups
  set.seed(seedInd[i])
  cluster <- kmeans(x = genoMat, centers = nCluster)$cluster

  indSel <- sapply(1:nCluster, function(l) {
    if (length(cluster[cluster == l]) == 1) {
      return(cluster[cluster == l])
    } else {
      groupL <- names(cluster[cluster == l])
      top <- sort(trueU[groupL, ], decreasing = T)[1]
      return(names(top))
    }
  })

  # save the information of 4 original lines
  resOriginList <- list(pheno = NULL,
                        genoMat = genoMat[indSel, ],
                        predU = trueU[indSel, , drop = F],
                        trueU = trueU[indSel, , drop = F],
                        generation = -1)
  saveRDS(resOriginList, paste0(dirSave, i, "_resOriginList.rds"))


  set.seed(seedInd[i])
  indSelOrder <- sample(indSel)
  crossTable <- data.frame(ind1 = indSelOrder[1:2],
                           ind2 = indSelOrder[3:4],
                           n = 1,
                           names = paste0("F1N",
                                          formatC(1:2,
                                                  width = 3, flag = "0")))

  # creating the F1 population
  f1Pop <- population$new(name = "F1",
                          inds = makeCrosses(crosses = crossTable, pop = initPop))

  # creating the 4-way cross population
  genoMat <- f1Pop$genoMat
  cand <- rownames(genoMat)

  crossTable <- data.frame(ind1 = cand[1],
                           ind2 = cand[2],
                           n = 150,
                           names = paste0("F2N001"))

  f2Pop <- population$new(name = "F2",
                          inds = makeCrosses(crosses = crossTable, pop = f1Pop))

  genoMat <- f2Pop$genoMat

  trueU <- c((genoMat - 1) %*% markerEffectTrue)
  trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
  rownames(trueU) <- rownames(genoMat)

  for (heritability in heritabilityVec) {
    # heritability <- heritabilityVec[1]
    # preparing a save directory
    dirSaveEach <- paste0(dirSave, "h2_", heritability, "/")
    if (!dir.exists(dirSaveEach)) {
      dir.create(dirSaveEach)
    }


    # calculating the residual variance (sigmaE) based on the heritability
    sigmaG <- var(trueU) * (length(trueU) - 1) / length(trueU)
    sigmaE <- (sigmaG / heritability) - sigmaG


    # simulating the phenotypic values
    set.seed(seedInd[i])
    phenoMat <- GetPheno(uMat = trueU,
                         sigmaE = sigmaE,
                         nPheno = nPheno)

    # building a genomic prediction model
    GPres <- GP(phenoMat = phenoMat, genoMat = genoMat - 1)
    predU <- GPres$uPredMat
    markerEffectEstimated <- GPres$markerEffPredMat
    beta <- GPres$betaEstimated

    # save the f2 information
    resList <- list(pheno = phenoMat,
                    genoMat = genoMat,
                    predU = predU,
                    trueU = trueU,
                    markerEffectEstimated = markerEffectEstimated,
                    beta = beta)
    saveRDS(resList, paste0(dirSaveEach, i, "_f2Pop.rds"))

    # generating the initial population
    selfTable <- data.frame(ind1 = rownames(f2Pop$genoMat),
                            ind2 = rownames(f2Pop$genoMat),
                            n = 1,
                            names = paste0("F3",
                                           str_sub(rownames(f2Pop$genoMat), start = 3)))
    f3Pop <- population$new(name = "F3",
                            inds = makeCrosses(crosses = selfTable, pop = f2Pop))
    saveRDS(f3Pop, file = paste0(dirSaveEach, i, "_initPop.rds"))

    resEachList <- list(pheno = NULL,
                        genoMat = f3Pop$genoMat,
                        predU = (f3Pop$genoMat - 1) %*% markerEffectEstimated + beta,
                        trueU = (f3Pop$genoMat - 1) %*% markerEffectTrue,
                        markerEffectEstimated = markerEffectEstimated,
                        beta = beta,
                        generation = 0)
    saveRDS(resEachList, paste0(dirSaveEach, i, "_resList_0.rds"))
  }

  # save other information
  genoPath <- paste0(simPath, i, "_genoQTL.csv")
  genoQTL <- read.csv(genoPath)
  genoSortedQTL <- cbind(id = genoQTL[, 1],
                         genoQTL[, sort(colnames(genoQTL)[2:ncol(genoQTL)])])
  write.csv(genoSortedQTL,
            file = paste0(dirSave, i, "_genoQTL.csv"), row.names = F)

  fileNames <- paste0(simPath, i,
                      c("_markerEffects.csv",
                        "_map.csv",
                        "_phenoQTL.csv"))
  fileNamesNew <- gsub(pattern = simPath, replacement = dirSave, fileNames)
  lapply(1:length(fileNames), function(l) {
    file.copy(from = fileNames[l], to = fileNamesNew[l])
  })
}
stopCluster(cl)
