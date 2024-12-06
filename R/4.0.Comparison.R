##############################################
# Title : 13.0.Comparison
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
methodFactor <- c("GS", "OCS", "CPS")
h2List <- c(0.3, 0.6)
methodColor <- c("gray20", "#4169E1", "#FF0000")

target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c("top", "variance")
generationInd1 <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))

target2 <- "RGS"
baseFileName2 <- "geneticValMat.csv"
index2 <- c("top", "variance")
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

index3 <- c("WeightedFixed", "WeightedLost", "WeightedNotFixed")

# 1.2. Save dir
for (h2 in h2List) {
  # h2 <- h2List[2]
  dirSave <- paste0("results/4.0.Comparison/h2_", h2, "/")
  if (!dir.exists(dirSave)) {
    dir.create(dirSave, recursive = T)
  }

  # 1.3. Read data
  gsPath <- paste0("results/3.0.GS/h2_", h2, "/")
  ocsPath <- paste0("results/3.1.OCS/h2_", h2, "/60_1_0.01/")
  cpsPath <- paste0("results/3.2.CPS/h2_", h2, "/")

  pathList <- list(list(method = "GS", path = gsPath),
                   list(method = "OCS", path = ocsPath),
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

  index3 <- c("WeightedFixed", "WeightedLost", "WeightedNotFixed")
  SimSumAllele(pathList = pathList,
               generationInd = generationInd2,
               index = index3,
               methodFactor = methodFactor,
               methodColor = methodColor)

  index3 <- c("FixedSum", "NotFixedSum", "Mean")
  index <- index3
  generationInd <- generationInd2
  resList <- lapply(pathList, function(pathEach) {
    # pathEach <- pathList[[1]]
    method <- pathEach$method
    path <- pathEach$path

    # sort the file based on file name
    fileName <- list.files(path, pattern = "desiredFreqQTN.csv")
    fileNum <- as.numeric(str_sub(fileName, start = 1, end = -nchar(fileName[1])))
    fileOrder <- order(fileNum)

    # gathering the results files of "F4" and "Recurrent genomic selection"
    resList <- list.files(path, pattern = "desiredFreqQTN.csv", full.names = T)[fileOrder]
    resArray <- array(data = NA, dim = c(length(generationInd), length(index), nRep))
    dimnames(resArray) <- list(generationInd,
                               index,
                               1:nRep)

    for (i in 1:nRep) {
      # i <- 1
      resMat0 <- as.matrix(read.csv(resList[[i]], row.names = 1))
      markerEffectAll <- read.csv(paste0("results/2.0.simulationSetting150/", i, "_markerEffects.csv"),
                                  row.names = 1)
      markerEffectSel <- markerEffectAll[rownames(resMat0), ]
      markerEffectAbs <- abs(markerEffectSel)

      resMat <- apply(resMat0, 2, function(x) {
        # x <- resMat0[, 10]
        marker <- x * 2
        marker[markerEffectSel < 0] <- 2 - marker[markerEffectSel < 0]
        fix <- sum((marker * markerEffectSel)[(marker == 2)])
        mid <- sum((marker * markerEffectSel)[(marker != 2)])
        mu <- sum(marker * markerEffectSel)
        return(c(Fixed = fix, NotFixed = mid, Mean = mu))
      })
      resArray[, , i] <- t(resMat[, 1:length(generationInd)])
    }

    resMeanMat <- apply(resArray, c(1, 2), mean, na.rm = T)
    resSeMat <- apply(resArray, c(1, 2), std.error, na.rm = T)

    generation <- as.numeric(str_sub(rownames(resMeanMat), start = 2))
    generationRep <- rep(generation, ncol(resMeanMat))
    methodRep <- rep(method, length(generationRep))
    indexRep <- rep(colnames(resMeanMat), each = nrow(resMeanMat))
    dfEach <- data.frame(Generation = generationRep,
                         Method = methodRep,
                         Index = indexRep,
                         Value = c(resMeanMat),
                         SE = c(resSeMat))
    return(resDf = dfEach)
  })

  resDf <- do.call(rbind, resList)
  resDf$Method <- factor(resDf$Method, levels = methodFactor)
  col <- methodColor

  for (indexEach in index) {
    # indexEach <- index[1]
    resDfEach <- resDf[resDf$Index == indexEach, ]
    comparisonList <- tapply(resDfEach[, "Value"], resDfEach$Method, function(x) {
      resDfEach[resDfEach$Method == "CPS", "Value"] / x
    })
    comparisonMat <- do.call(rbind, comparisonList)
    colnames(comparisonMat) <- unique(resDfEach$Generation)
    write.csv(comparisonMat, paste0(dirSave, indexEach, "_compare.csv"))

    g <- ggplot(resDfEach, aes(x = Generation, y = Value, color = Method)) +
      geom_line(size = 1.2) +
      theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 24),
            axis.title.y = element_text(size = 6)) +
      ylim(c(-3, 60)) +
      scale_color_manual(values = col)

    png(paste0(dirSave, indexEach, ".png"),
        width = 1440, height = 1440, res = 214)
    print(g)
    dev.off()

  }
}

