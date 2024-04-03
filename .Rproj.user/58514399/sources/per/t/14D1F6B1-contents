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
packages <- c("breedSimulatR", "stringr", "qtl")
source(file = "R/1.0.function.R")

# 1.2. Save dir
dirSave <- "results/2.1.simulationSetting150/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

seedIndCsv <- paste0("results/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:50000, 1000, replace = F)
  write.csv(x = seedInd, file = paste0(seedIndCsv))
}

cl <- makeCluster(60)
registerDoParallel(cl)
foreach(i = 1:300, .export = ls(envir = parent.frame()), .packages = packages) %dopar% {
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

  trueU <- c(genoMat %*% markerEffectTrue)
  trueU <- matrix(trueU, nrow = nrow(genoMat), ncol = nTrait)
  rownames(trueU) <- rownames(genoMat)


  # Clustering into 4 groups
  set.seed(seedInd[i])
  cluster <- kmeans(x = genoMat, centers = 4)$cluster

  indSel <- sapply(1:4, function(l) {
    if (length(cluster[cluster == l]) == 1) {
      return(cluster[cluster == l])
    } else {
      groupL <- names(cluster[cluster == l])
      top <- sort(trueU[groupL, ], decreasing = T)[1]
      return(names(top))
    }
  })

  set.seed(seedInd[i])
  indSelOrder <- sample(indSel)

  crossTable <- data.frame(ind1 = indSelOrder[1:2],
                           ind2 = indSelOrder[3:4],
                           n = 1,
                           names = paste0("F1N",
                                          formatC(1:2,
                                                  width = 3, flag = "0")))

  # creating the F1 population
  f1Pop <- population$new(name = "F1 offspring",
                          inds = makeCrosses(crosses = crossTable, pop = initPop))
  saveRDS(f1Pop, file = paste0(dirSave, i, "_F1Pop.rds"))

  # creating the 4-way cross population
  genoMat <- f1Pop$genoMat
  cand <- rownames(genoMat)

  crossTable <- data.frame(ind1 = cand[1],
                           ind2 = cand[2],
                           n = 150,
                           names = paste0("F2N001"))

  g2Pop <- population$new(name = "F2 offspring",
                          inds = makeCrosses(crosses = crossTable, pop = f1Pop))
  saveRDS(g2Pop, file = paste0(dirSave, i, "_initPop.rds"))

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
