##############################################
# Title : 4.0.ComparisonOCS
# Author : Kengo Sakurai
# Date : 2024-03-06
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(plotrix)
library(stringr)
library(ggplot2)
library(abind)
source(file = "R/1.0.function.R")

# 1.1. Parameters
nGeneration <- 20 # the number of generations
nRep <- 100 # the number of replication
tRep <- rep(targetGenerationVec, each = length(HeStarRatioVec) * length(degreeVec))
sRep <- rep(rep(degreeVec, each = length(HeStarRatioVec)), length(targetGenerationVec))
heRep <- rep(HeStarRatioVec, length(targetGenerationVec) * length(degreeVec))

methodFactor <- paste0(tRep, "_", sRep, "_", heRep)
methodColor <- c("gray20", "#FF0000", "green4")


target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c("top", "top10Mean")
generationInd1 <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))

target2 <- "RGS"
baseFileName2 <- "geneticValMat.csv"
index2 <- c("top", "top10Mean", "variance")
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

# 1.2. Save dir
dirSave <- "results/4.0.ComparisonOCS/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# 1.3. Read data
dirPath <- paste0("results/3.1.OCS/", methodFactor, "/")
pathList <- list()
for (i in 1:length(dirPath)) {
  pathList[[i]] <- list(t = tRep[i],
                        s = sRep[i],
                        he = heRep[i],
                        path = dirPath[i])
}

###### 2. Analysis ######
SimSumOCS(pathList = pathList,
          target = target1,
          baseFileName = baseFileName1,
          index = index1,
          generationInd = generationInd1,
          methodFactor = methodFactor,
          methodColor = methodColor)

SimSumOCS(pathList = pathList,
          target = target2,
          baseFileName = baseFileName2,
          index = index2,
          generationInd = generationInd2,
          methodFactor = methodFactor,
          methodColor = methodColor)
