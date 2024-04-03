##############################################
# Title : 4.1.Comparison
# Author : Kengo Sakurai
# Date : 2023-08-16
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
nRep <- 300 # the number of replication
methodFactor <- c("GS", "OCS1", "OCS2", "CPS")
methodColor <- c("gray20", "green4", "#4169E1", "#FF0000")

target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c("top", "top10Mean")
generationInd1 <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))

target2 <- "RGS"
baseFileName2 <- "geneticValMat.csv"
index2 <- c("top", "top10Mean", "variance")
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

index3 <- c("WeightedFixed", "WeightedLost", "WeightedNotFixed")

# 1.2. Save dir
dirSave <- "results/4.1.Comparison/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# 1.3. Read data
gsPath <- "results/3.0.GS/"
ocsPath1 <- "results/3.1.OCS/20_1_0.01/"
ocsPath2 <- "results/3.1.OCS/20_1_0.3/"
cpsPath <- "results/3.2.CPS/"

pathList <- list(list(method = "GS", path = gsPath),
                 list(method = "OCS1", path = ocsPath1),
                 list(method = "OCS2", path = ocsPath2),
                 list(method = "CPS", path = cpsPath))


###### 2. Analysis ######
SimSum(pathList = pathList,
       target = target1,
       baseFileName = baseFileName1,
       index = index1,
       generationInd = generationInd1,
       methodFactor = methodFactor,
       methodColor = methodColor)

SimSum(pathList = pathList,
       target = target2,
       baseFileName = baseFileName2,
       index = index2,
       generationInd = generationInd2,
       methodFactor = methodFactor,
       methodColor = methodColor)

SimSumAllele(pathList = pathList,
             generationInd = generationInd2,
             index = index3,
             methodFactor = methodFactor,
             methodColor = methodColor)

