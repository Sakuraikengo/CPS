# 1.1. Parameters
heritabilityVec <- c(0.3, 0.6) # heritability in each simulation
nProg <- 15 # the number of progeny from each cross
nProgSelf <- 50 # the number of progeny for selfing
probSelf <- 1 - (1 / nProgSelf) # probability (which will decide the selection intensity for selfing)
intensitySelf <- qnorm(probSelf, 0, 1) # selection intensity (selfing)
nProgAll <- nProg * nProgSelf # the number of progeny from each cross (considering selfing)
probCross <- 1 - (1 / nProgAll) # probability (which will decide the selection intensity for cross)
intensityCross <- qnorm(probCross, 0, 1) # selection intensity (cross)
probNext <- 1 - (1 / nProg) # probability (which will decide the selection intensity for next generation)
intensityNext <- qnorm(probNext, 0, 1) #selection intensity (next generation)
capa <- 2 # how many times can we use each genotype for cross
nCrosses <- 10 # the number of crosses
nSelf <- 2 # the number of selfing genotypes for releasing varieties
nGeneration <- 60 # the number of generations
nTrait <- 1 # the number of trait
nRep <- 300 # the number of breeding simulations
f <- 8 # the generation for evaluation (F8)
initInd <- 4 # initial number for searching init OC
maxIter = 1000 # max iteration for searching init OC
HeStarRatioVec <- c(0.01) # final genetic diversity ratio
degreeVec <- c(1) # mathmatical degree for the penalty term for OCS
targetGenerationVec <- c(60) # target generation in OCS

# collecting phenotypic data
GetPheno <- function(uMat, sigmaE, nPheno) {
  phenoList <- lapply(1:length(sigmaE), function(x) {
    residual <- mvrnorm(n = nPheno, mu = rep(0, nrow(uMat)), Sigma = diag(sigmaE[x], nrow = nrow(uMat)))
    phenoVec <- rep(uMat[, x], each = nPheno) + c(residual)
    pheno <- matrix(phenoVec, ncol = 1)
    rownames(pheno) <- rep(rownames(uMat), each = nPheno)
    return(pheno)
  })
  phenoMat <- do.call(cbind, phenoList)
  return(phenoMat)
}

GP <- function(phenoMat, genoMat) {
  amat <- calcGRM(genoMat = genoMat, methodGRM = "A.mat")
  Z <- design.Z(pheno.labels = rownames(phenoMat), geno.names = rownames(amat))
  ZETA <- list(A = list(Z = Z, K = amat))

  markerEffPredMat <- matrix(0, nrow = ncol(genoMat), ncol = ncol(phenoMat))
  uPredMat <- matrix(0, nrow = nrow(genoMat), ncol = ncol(phenoMat))
  rownames(markerEffPredMat) <- colnames(genoMat)
  rownames(uPredMat) <- rownames(genoMat)
  colnames(markerEffPredMat) <- colnames(uPredMat) <- colnames(phenoMat)

  for (traitInd in 1:ncol(phenoMat)) {
    # traitInd <- 1
    model <- EMM.cpp(y = phenoMat[, traitInd], ZETA = ZETA, n.core = 1)
    # h2 <- model$Vu / (model$Vu + model$Ve)

    # estimate marker effects from estimated genotypic values
    if (min(eigen(amat)$values) < 1e-08) {
      diag(amat) <- diag(amat) + 1e-06
    }
    mrkEffects <- t(genoMat) %*% solve(genoMat %*% t(genoMat)) %*% model$u
    mrkEffects <- c(mrkEffects)

    # input the result into matrix
    markerEffPredMat[, traitInd] <- mrkEffects
    uPredMat[, traitInd] <- model$u + model$beta
  }
  return(list(uPredMat = uPredMat,
              markerEffPredMat = markerEffPredMat,
              betaEstimated = model$beta))
}

# calculating the frequency of desirable qtns and save
CalcFreqQTNs <- function(genoMat, qtn, qtnMinus, dirSave, i, generation) {
  genoMatSel <- genoMat[, qtn]
  freq <- apply(genoMatSel, 2, function(x) {
    sum(x / 2) / length(x)
  })
  desiredFreq <- abs(qtnMinus - freq)
  write.csv(desiredFreq, paste0(dirSave, i, "_desiredFreqQTN_", generation, ".csv"))
}

# summarize the results of each simulation
SummarizeResults <- function(dataPathList, rowNamesData, plotTitle, saveFileName) {
  resList <- lapply(dataPathList, function(dataPathEach) {
    # dataPathEach <- dataPathList[[2]]
    resGenereationEach <- readRDS(dataPathEach)
    generation <- resGenereationEach$generation
    geneticValueTrue <- resGenereationEach$trueU
    plotPath0 <- gsub(dataPathEach, pattern = "resList", replacement = "u")
    plotPath <- gsub(plotPath0, pattern = ".rds", replacement = ".png")
    png(plotPath)
    hist(c(geneticValueTrue), breaks = 10, freq = F,
         xlab = "Trait1", ylab = "Frequency",
         xlim = c(-5, 100), ylim = c(0, 1),
         main = paste0(plotTitle, generation))
    dev.off()

    top10 <- sort(geneticValueTrue[, 1], decreasing = T)[1:10]
    top <- sort(geneticValueTrue[, 1], decreasing = T)[1]
    top10Mean <- mean(top10)

    gVariance <- var(geneticValueTrue) * (length(geneticValueTrue) - 1) / length(geneticValueTrue)
    res <- c(top, top10Mean, gVariance)
    names(res) <- c("top", "top10Mean", "variance")
    return(res)
  })
  resMat <- do.call(rbind, resList)
  rownames(resMat) <- rowNamesData
  write.csv(resMat, file = saveFileName)
}

# summary of the desired allele frequency
SummarizeFreq <- function(freqFilePathList, nGeneration, i, cross) {
  freqPath <- paste0(dirSave, i, "_desiredFreqQTN.csv")
  if (!file.exists(freqPath)) {
  freqList <- lapply(freqFilePathList, function(freqFilePathEach) {
    # freqFilePathEach <- freqFilePathList[[1]]
    freqEach <- read.csv(freqFilePathEach)
    freqVec <- freqEach[, 2]
    names(freqVec) <- freqEach$X
    return(freqVec)
  })
  freqMat <- do.call(cbind, freqList)
  write.csv(freqMat, file = paste0(dirSave, i, "_desiredFreqQTN.csv"))
  file.remove(freqFilePathList)
  } else {
    freqMat <- read.csv(freqPath, row.names = 1)
  }

  chr <- str_sub(rownames(freqMat), 4, 5)
  png(filename = paste0(dirSave, i, "_desiredFreqQTN.png"))
  matplot(t(freqMat), type = "l", col = chr)
  dev.off()

  notFixedInd <- freqMat[, nGeneration] > 0 & freqMat[, nGeneration] < 1
  notFixedQTN <- rownames(freqMat)[notFixedInd]
  for (eachChrInd in 1:length(cross$geno)) {
    # eachChrInd <- 1
    names(cross$geno[[eachChrInd]]$map)[!(names(cross$geno[[eachChrInd]]$map) %in% notFixedQTN)] <- NA
  }
  # visualize the not fixed QTNs
  png(filename = paste0(dirSave, i, "_notFixedQTN.png"))
  plotMap(cross, show.marker.names = T)
  dev.off()
}

# calculating the D(2) matrix based on recombination frequency
CalcD2 <- function(ObjectMap, markerEffectMat, nTrait, k) {
  # getting the number of chr
  nChr <- length(unique(ObjectMap$chr))

  # separating to each chr
  D2List <- lapply(unique(ObjectMap$chr), function(eachChr) {
    # eachChr <- unique(ObjectMap$chr)[2]
    mapEach <- ObjectMap[ObjectMap$chr == eachChr, ]
    effectEach <- markerEffectMat[ObjectMap$chr == eachChr, ]
    markerName <- rownames(mapEach)

    # calculating the distance of each marker set
    myDist <- sapply(1:nrow(mapEach),
                     function(x) abs(mapEach$linkMapPos[x] - mapEach$linkMapPos))

    # convert the distance to the recombination frequency (Haldane)
    c1 <- 0.5 * (1 - exp(-2 * (myDist / 100)))

    # r for the future F"k+1" generations
    ck <- 2*c1*(1 - ((0.5)^k)*(1 - 2*c1)^k) / (1 + 2*c1)

    # r for the future F"k" generations (for selfing)
    l <- k - 1
    cSelf <- 2*c1*(1 - ((0.5)^l)*(1 - 2*c1)^l) / (1 + 2*c1)

    D1_1 <- (1 - ck) * (1 - 2*c1) / 4
    D1_2 <- (1 - 2*ck - (0.5*(1 - 2*c1))^(k)) / 4
    D1Self_1 <- (1 - cSelf) * (1 - 2*c1) / 4
    D1Self_2 <- (1 - 2*cSelf - (0.5*(1 - 2*c1))^(l)) / 4

    # calculating the D2
    D2_1 <- diag(effectEach) %*% D1_1 %*% diag(effectEach)
    D2_2 <- diag(effectEach) %*% D1_2 %*% diag(effectEach)
    D2Self_1 <- diag(effectEach) %*% D1Self_1 %*% diag(effectEach)
    D2Self_2 <- diag(effectEach) %*% D1Self_2 %*% diag(effectEach)

    rownames(D2_1) <- rownames(D2_2) <- rownames(D2Self_1) <- rownames(D2Self_2) <- markerName
    colnames(D2_1) <- colnames(D2_2) <- colnames(D2Self_1) <- colnames(D2Self_2) <- markerName

    return(list(D2 = list(D2_1, D2_2),
                D2Self = list(D2Self_1, D2Self_2)))
  })
  D2 <- lapply(D2List, function(D2ListEach) {
    D2ListEach$D2
  })
  D2Self <- lapply(D2List, function(D2ListEach) {
    D2ListEach$D2Self
  })
  return(list(D2 = D2, D2Self = D2Self))
}


# Calculating the progeny variance for each cross
CalcVarMat <- function(D2, combInd, gametArray, nTrait) {
  sigmaMatEachChr <- sapply(D2, function(D2each) {
    # D2each <- D2[[1]]
    D2_1_each <- D2each[[1]]
    D2_2_each <- D2each[[2]]
    markerName <- rownames(D2_1_each)

    markerMat_1 <- gametArray[, markerName, 1] * 2
    markerMat_2 <- gametArray[, markerName, 2] * 2

    mu1_1 <- diag(markerMat_1 %*% D2_1_each %*% t(markerMat_1)) / 4
    mu1_2 <- diag(markerMat_2 %*% D2_1_each %*% t(markerMat_2)) / 4
    gamma1 <- diag(markerMat_1 %*% D2_1_each %*% t(markerMat_2)) / 4
    gamma2 <- diag(markerMat_2 %*% D2_1_each %*% t(markerMat_1)) / 4

    gammaMat2_1 <- markerMat_1 %*% D2_2_each %*% t(markerMat_1) / 4
    gammaMat2_2 <- markerMat_2 %*% D2_2_each %*% t(markerMat_2) / 4
    mu2_1 <- diag(gammaMat2_1)
    mu2_2 <- diag(gammaMat2_2)
    gammaMat2_3 <- markerMat_1 %*% D2_2_each %*% t(markerMat_2) / 4
    gammaMat2_4 <- markerMat_2 %*% D2_2_each %*% t(markerMat_1) / 4

    p1 <- combInd[1, ]
    p2 <- combInd[2, ]

    eq1 <- mu1_1[p1] + mu1_2[p1] + mu1_1[p2] + mu1_2[p2]
    eq2 <- gamma1[p1] + gamma1[p2] + gamma2[p1] + gamma2[p2]
    eq3 <- mu2_1[p1] + mu2_2[p1] + mu2_1[p2] + mu2_2[p2]
    eq4 <- sapply(1:length(p1), function(x) {
      p1Each <- p1[x]
      p2Each <- p2[x]

      x1_1 <- gammaMat2_1[p1Each, p2Each]
      x1_2 <- gammaMat2_2[p1Each, p2Each]
      x1_3 <- gammaMat2_3[p1Each, p2Each]
      x1_4 <- gammaMat2_4[p1Each, p2Each]
      x2_1 <- gammaMat2_1[p2Each, p1Each]
      x2_2 <- gammaMat2_2[p2Each, p1Each]
      x2_3 <- gammaMat2_3[p2Each, p1Each]
      x2_4 <- gammaMat2_4[p2Each, p1Each]
      return(x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4)
    })
    sigmaEach <- eq1 - eq2 + 2*eq3 - eq4
    return(sigmaEach)
  })
  sigmaVec <- apply(sigmaMatEachChr, 1, sum)
  sigmaMat <- matrix(sigmaVec, ncol = nTrait)
  return(sigmaMat)
}

# Extracting the best crosses using Linear Programing(LP)
LP <- function(predProgMat, # the potential of each cross
               indName, # the name of individual
               comb, # the matrix of each cross
               capa, # capacity for each edge (individual)
               nCrosses, # the number of crosses
               nProg, # the number of progenies for each cross
               generation # count of the generation
) {
  # matrix to vector
  dataVec <- c(-predProgMat)
  crossNameSel <- rownames(predProgMat)

  # edges
  ed1 <- paste0("L_", indName)
  ed2 <- paste0("R_", indName)
  ed <- c(ed1, ed2)

  # nodes
  s_es <- paste0("S-", ed1)
  t_es <- paste0(ed2, "-T")
  c_es <- paste0("S-L_", comb[1, ],
                 "-R_",
                 comb[2, ], "-T")
  es <- c(s_es, t_es, c_es)

  # cost
  s_cs <- rep(0, length(s_es))
  t_cs <- rep(0, length(t_es))
  c_cs <- dataVec
  cs <- c(s_cs, t_cs, c_cs)

  # limitations
  # Normally, we make a design matrix, but it will be so large.
  # we made a M times 3 matrix (M depends on the limitation)
  # first column means the location of row
  # second column means the location of column
  # third column means the value of limitation
  stc_lim <- matrix(c(rep(1:length(cs), 2),
                      rep(1, length(cs))), ncol = 3)
  s_lim <- matrix(c(rep(length(cs) + 1, length(s_cs)),
                    1:length(s_cs),
                    rep(1, length(s_cs))), ncol = 3)
  t_lim <- matrix(c(rep(length(cs) + 2, length(s_cs)),
                    (length(s_cs)+1):(length(s_cs) + length(t_cs)),
                    rep(1, length(s_cs))), ncol = 3)
  t_lim_end <- t_lim[1, 1]

  # the amount of flow of input and output from each edge is must be the same
  node_num <- length(s_es) + length(t_es) + 2*length(c_es)
  ed_lim <- matrix(NA,
                   nrow = node_num,
                   ncol = 3)
  for (i in 1:length(ed)) {
    # i <- 1
    # extracting each edge
    ed_each <- ed[i]
    ed_lim_ind <- t_lim_end + i

    # index for searching the nodes (edge to edge) which come from "ed_each"
    c_ind <- paste0("-", ed_each, "-")

    # start to edge (edge to end) (-1)
    start <- min(which(is.na(ed_lim[, 1])))
    ed_lim[start, ] <- c(ed_lim_ind, i, -1)

    # edge to edge (1)
    ed_lim_each <- which(grepl(pattern = c_ind, x = es))
    if (length(ed_lim_each) < 1) {
      next
    }

    # add the limitation (the amount of flow is the same)
    ed_lim_mat <- matrix(c(rep(ed_lim_ind, length(ed_lim_each)),
                           ed_lim_each,
                           rep(1, length(ed_lim_each))), ncol = 3)
    start <- min(which(is.na(ed_lim[, 1])))
    ed_lim[start:(start + length(ed_lim_each) - 1), ] <- ed_lim_mat
  }

  ed_lim_end <- ed_lim[nrow(ed_lim), 1] + 1

  # the total amount of flow of specific edge
  # (start to the edge & the edge to end)
  # must be limited to "capa"
  sym_lim_row_ind <- ed_lim_end:(ed_lim_end + length(s_cs) - 1)
  sym_lim_start <- matrix(c(sym_lim_row_ind,
                            1:length(s_cs),
                            rep(1, length(s_cs))), ncol = 3)
  sym_lim_end_ind <- (1+length(s_cs)):(length(s_cs)+length(t_cs))
  sym_lim_end <- matrix(c(sym_lim_row_ind,
                          sym_lim_end_ind,
                          rep(1, length(t_cs))), ncol = 3)
  sym_lim <- rbind(sym_lim_start, sym_lim_end)
  lim <- rbind(stc_lim, s_lim, t_lim, ed_lim, sym_lim)

  # lim <- rbind(stc_lim, s_lim, t_lim, ed_lim, sym_lim)
  eq <- c(rep("<=", length(c(s_cs, t_cs))),
          rep("<=", length(c_cs)),
          rep("==", 2 + length(ed)),
          rep("<=", length(s_cs)))

  # capacity
  s_ws <- rep(capa, length(s_es))
  t_ws <- rep(capa, length(t_es))
  c_ws <- rep(1, length(c_es))
  flow <- rep(nCrosses, 2)
  ed_ws <- rep(0, length(ed))
  sym_ws <- rep(capa, length(s_cs))
  ws <- c(s_ws, t_ws, c_ws, flow, ed_ws, sym_ws)
  res <- lp(direction = "min",
            objective.in = cs,
            dense.const = lim,
            const.dir = eq,
            const.rhs = ws,
            all.int = T)
  ind <- as.character(round(res$solution))
  # sum(dataVec[ind[(2 * n + 1):length(cs)] == 1])
  selCrosses0 <- es[ind == 1]
  selCrosses0 <- selCrosses0[grepl(pattern = "^S.{1,}-T$", x = selCrosses0)]
  selCrosses <- str_sub(selCrosses0, 5, -3)
  selCrossesList <- str_split(selCrosses, pattern = "-R_")
  bestCrossMat <- do.call(rbind, selCrossesList)
  if (max(table(c(bestCrossMat))) > 2) {
    print("Each genotype must not be used more than nCrosses times")
  } else {
    crossTable <- data.frame(ind1 = bestCrossMat[, 1],
                             ind2 = bestCrossMat[, 2],
                             n = nProg,
                             names = paste0("C", generation, "N",
                                            formatC(1:nrow(bestCrossMat), width = 3, flag = "0")))
    return(crossTable)
  }
}

CalcProgMean <- function(trueU) {
  trueProgMeanMat <- apply(trueU, 2, function(trueUeach) {
    trueProgMeanAll <- (rep(trueUeach, each = length(trueUeach)) + rep(trueUeach, length(trueUeach))) / 2
    trueProgMeanMat <- matrix(trueProgMeanAll, nrow = length(trueUeach), ncol = length(trueUeach))
    rownames(trueProgMeanMat) <- colnames(trueProgMeanMat) <- names(trueUeach)
    trueProgMean <- trueProgMeanMat[lower.tri(trueProgMeanMat)]
    return(trueProgMean)
  })
  return(trueProgMeanMat)
}

ExtractGamet <- function(nowPop, qtn) {
  gametArrayList <- lapply(nowPop$inds, function(x) {
    # x <- initPop$inds[[1]]
    gamet <- do.call(cbind, x$haplo$values)
    gamet <- gamet[, qtn]
    array(data = c(t(gamet)), dim = c(1, ncol(gamet), 2))
  })
  gametArray <- do.call(abind, c(gametArrayList, along = 1))
  return(gametArray)
}

ExpectedFnValue <- function(sigmaMat, predProgMean, intensityCross) {
  ExpectedProgFn <- predProgMean + intensityCross*sqrt(sigmaMat)
  return(ExpectedProgFn)
}


CreateF8 <- function(selfSelected, nowPop, nProgSelf, markerEffectTrue, nTrait, generation) {
  f8GenerationInd <- str_split(selfSelected[1], pattern = "N")[[1]][1]
  f8Generation <- as.numeric(str_sub(f8GenerationInd, 2)) + 8
  f8Generation <- formatC(f8Generation, width = 2, flag = "0")
  selfTable <- data.frame(ind1 = selfSelected,
                          ind2 = selfSelected,
                          n = nProgSelf,
                          names = paste0("C", f8Generation, "F2N",
                                         formatC(1:length(selfSelected), width = 3, flag = "0")))
  f2Pop <- population$new(name = "F2 offspring",
                          inds = makeCrosses(crosses = selfTable, pop = nowPop))
  selfNowPop <- f2Pop
  for (l in 2:7) {
    selfTable <- data.frame(ind1 = rownames(selfNowPop$genoMat),
                            ind2 = rownames(selfNowPop$genoMat),
                            n = 1,
                            names = paste0("C", f8Generation, "F", l+1, "N",
                                           formatC(1:nrow(selfNowPop$genoMat), width = 3, flag = "0")))
    selfNextPop <- population$new(name = paste0("F", l+1, " offspring"),
                                  inds = makeCrosses(crosses = selfTable, pop = selfNowPop))
    selfNowPop <- selfNextPop
  }
  genoMatNew <- selfNowPop$genoMat
  geneticValTrueNew <- (genoMatNew - 1) %*% markerEffectTrue
  return(list(genoMat = genoMatNew,
              trueU = geneticValTrueNew,
              generation = generation))
}


# Calculating the He0 (the initial neutral diversity)
CalcHe0 <- function(genoMat) {
  alleleFreq <- apply(genoMat, 2, function(x) {
    sum(x/2) / length(x)
  })
  He0 <- (2 * t(alleleFreq) %*% (1 - alleleFreq)) / ncol(genoMat)
}

SearchInitOC <- function(predProgMean, comb, combName, K, Z1, Z2, He, initInd, maxIter, nCrosses, capa) {
  # the index for the selected crosses
  selectedInd <- sample(1:ncol(comb), initInd)
  while (initInd <= nCrosses & length(selectedInd) < nCrosses) {
    initInd <- initInd + 1
    iter <- 1 # the number of iteration
    print(initInd)
    # the index for the selected crosses
    selectedInd <- sample(1:ncol(comb), initInd)
    while (iter < maxIter & length(selectedInd) < nCrosses) {
      selectedInd <- sample(1:ncol(comb), initInd)
      selectedCross <- comb[, selectedInd]
      selectedName <- c(selectedCross)
      # check "the same cross is not selected",
      # "each individual should not be used more than 2 times"
      if (max(table(selectedName)) > 2) {
        iter <- iter + 1
        next
      }

      # set the contribution of Parent1 and Parent2
      contribution1 <- contribution2 <- rep(0, length(combName))
      names(contribution1) <- names(contribution2) <- combName
      contribution1[selectedInd] <- 1 / 2
      contribution2[selectedInd] <- 1 / 2
      contribution <- (Z1 %*% contribution1 + Z2 %*% contribution2) / length(selectedInd)

      # genetic diversity
      DA <- 1 - (t(contribution) %*% K %*% contribution)

      # adding a cross which meets the criterion one by one
      for (cycleInd in (initInd+1):nCrosses) {
        # cycleInd <- initInd+1
        # calculating DB for all crosses at the same time
        n <- length(selectedInd)
        contribution1B <- rep(1 / 2, ncol(Z1))
        contribution2B <- rep(1 / 2, ncol(Z1))
        Y <- Z1 * contribution1B + Z2 * contribution2B

        DB_1 <- (1 - DA)
        DB_2_tmp <- ((1 / n)^2) * t(Y) %*% K * t(Y)
        DB_2 <- apply(DB_2_tmp, 1, sum)
        DB_3 <- (2 / n) * t(contribution) %*% K %*% Y
        DB <- 1 - (((n / (n + 1))^2) * (DB_1 + c(DB_2) + c(DB_3)))

        if (any(DB > He)) {
          # extracting the crosses which meet the criterion
          crossNewCandInd <- which(DB > He)
          crossOrder <- order(predProgMean[crossNewCandInd, ], decreasing = T)
          cand <- 1
          # selecting the crosses based on the predProgMean value
          while (cand <= length(crossOrder)) {
            # cand <- 1
            crossNewInd <- crossNewCandInd[crossOrder][cand]
            selectedCandInd <- c(selectedInd, crossNewInd)
            selectedCandCross <- comb[, selectedCandInd]
            selectedCandName <- c(selectedCandCross)
            # check "the same cross is not selected",
            # "each individual should not be used more than 2 times"
            if (max(table(selectedCandInd)) == 1 & max(table(selectedCandName)) <= capa) {
              # add the new cross to the selected individuals
              selectedInd <- selectedCandInd
              selectedCross <- selectedCandCross
              selectedName <- selectedCandName
              DA <- DB[crossNewCandInd[crossOrder][cand]]
              contribution1 <- contribution2 <- rep(0, length(combName))
              names(contribution1) <- names(contribution2) <- combName
              contribution1[selectedInd] <- 1 / 2
              contribution2[selectedInd] <- 1 / 2
              contribution <- (Z1 %*% contribution1 + Z2 %*% contribution2) / length(selectedInd)
              print("update!!")
              # stop the while iteration
              cand <- 1e10
              # stop()
            } else {
              # try next cross
              cand <- cand + 1
            }
          }
        } else {
          break
        }
        if (cand != 1e10) {
          # if the crosses don't meet the criterion, go to the next iter
          break
        }
      }
      iter <- iter + 1
    }
  }
  if (max(table(selectedCandInd)) == 1 & max(table(selectedCandName)) <= capa) {
    return(selectedInd)
  } else {
    return(NA)
  }
}

SearchOptimalOC <- function(selectedInd, predProgMean, comb, combName, K, Z1, Z2, He, nCrosses, capa) {
  selectedOptimalInd <- selectedInd
  valueOrderInd <- 1
  while (valueOrderInd <= nCrosses) {
    print(valueOrderInd)
    valueOrder <- order(predProgMean[selectedOptimalInd, ])
    valueSum <- sum(predProgMean[selectedOptimalInd, ])
    print(valueSum)

    selectedHighInd <- selectedOptimalInd[-valueOrderInd]

    selectedHighCross <- comb[, selectedHighInd]
    selectedHighName <- c(selectedHighCross)

    # set the contribution of Parent1 and Parent2
    contribution1 <- contribution2 <- rep(0, length(combName))
    names(contribution1) <- names(contribution2) <- combName
    contribution1[selectedHighInd] <- 1 / 2
    contribution2[selectedHighInd] <- 1 / 2
    contribution <- (Z1 %*% contribution1 + Z2 %*% contribution2) / length(selectedHighInd)

    # genetic diversity
    DA <- 1 - (t(contribution) %*% K %*% contribution)

    # calculating DB for all crosses at the same time
    n <- length(selectedHighInd)
    contribution1B <- rep(1 / 2, ncol(Z1))
    contribution2B <- rep(1 / 2, ncol(Z1))
    Y <- Z1 * contribution1B + Z2 * contribution2B

    DB_1 <- (1 - DA)
    DB_2_tmp <- ((1 / n)^2) * t(Y) %*% K * t(Y)
    DB_2 <- apply(DB_2_tmp, 1, sum)
    DB_3 <- (2 / n) * t(contribution) %*% K %*% Y
    DB <- 1 - (((n / (n + 1))^2) * (DB_1 + c(DB_2) + c(DB_3)))

    if (any(DB > He)) {
      # extracting the crosses which meet the criterion
      crossNewCandInd <- which(DB > He)
      crossOrder <- order(predProgMean[crossNewCandInd, ], decreasing = T)
      cand <- 1
      # selecting the crosses based on the predProgMean value
      while (cand <= length(crossOrder)) {
        # cand <- 1
        crossNewInd <- crossNewCandInd[crossOrder][cand]
        selectedCandInd <- c(selectedHighInd, crossNewInd)
        selectedCandCross <- comb[, selectedCandInd]
        selectedCandName <- c(selectedCandCross)
        check1 <- max(table(selectedCandInd)) == 1
        check2 <- max(table(selectedCandName)) <= capa
        check3 <- sum(predProgMean[selectedCandInd, ]) > valueSum
        check <- c(check1, check2, check3)
        # check "the same cross is not selected",
        # "each individual should not be used more than 2 times"
        if (all(check)) {
          # add the new cross to the selected individuals
          selectedOptimalInd <- selectedCandInd
          selectedCross <- selectedCandCross
          selectedName <- selectedCandName
          DA <- DB[crossNewCandInd[cand]]
          contribution1 <- contribution2 <- rep(0, length(combName))
          names(contribution1) <- names(contribution2) <- combName
          contribution1[selectedOptimalInd] <- 1 / 2
          contribution2[selectedOptimalInd] <- 1 / 2
          contribution <- (Z1 %*% contribution1 + Z2 %*% contribution2) / length(selectedOptimalInd)
          print("update!!")
          # stop the while iteration
          cand <- 1e10
          # stop()
        } else {
          # try next cross
          if (cand == length(crossOrder)) {
            valueOrderInd <- valueOrderInd + 1
          }
          cand <- cand + 1
        }
      }
    } else {
      valueOrderInd <- valueOrderInd + 1
    }
  }
  return(selectedOptimalInd)
}


SelectForSelfing <- function(trueU,
                              gametArray,
                              D2Self,
                              nTrait,
                              nSelf,
                              intensity) {
  selfInd <- rbind(1:nrow(trueU), 1:nrow(trueU))

  varSelfMat <- CalcVarMat(D2 = D2Self,
                            combInd = selfInd,
                            gametArray = gametArray,
                            nTrait = nTrait)
  varSelfMat[varSelfMat < 0] <- 0
  rownames(varSelfMat) <- rownames(trueU)
  fnValue <- trueU + intensity*sqrt(varSelfMat)
  selfSelected <- names(sort(fnValue[, 1], decreasing = T))[1:nSelf]
  return(selfSelected)
}



RandomMate <- function(cand, nCrosses, nProg, generation) {
  candRandom <- sample(x = cand, size = nCrosses, replace = F)
  candRandomRep <- rep(candRandom, each = 2)
  crossMat <- matrix(c(candRandomRep[2:length(candRandomRep)], candRandomRep[1]),
                     ncol = 2, byrow = T)
  crossTable <- data.frame(ind1 = crossMat[, 1],
                           ind2 = crossMat[, 2],
                           n = nProg,
                           names = paste0("C", generation, "N",
                                          formatC(1:nrow(crossMat), width = 3, flag = "0")))

}

SimSum <- function(pathList, target, baseFileName, index, generationInd, methodFactor, methodColor) {
  resList <- lapply(pathList, function(pathEach) {
    # pathEach <- pathList[[1]]
    method <- pathEach$method
    path <- pathEach$path

    # sort the file based on file name
    fileName <- list.files(path, pattern = baseFileName)
    fileNum <- as.numeric(str_sub(fileName, start = 1, end = -nchar(fileName[1])))
    fileOrder <- order(fileNum)

    # gathering the results files of "F4" and "Recurrent genomic selection"
    resList <- list.files(path, pattern = baseFileName, full.names = T)[fileOrder]
    resArray <- array(data = NA, dim = c(length(generationInd), length(index), nRep))
    dimnames(resArray) <- list(generationInd,
                               index,
                               1:nRep)

    for (i in 1:nRep) {
      # i <- 1
      resMat0 <- as.matrix(read.csv(resList[[i]], row.names = 1))
      resMat0 <- resMat0[1:dim(resArray)[1], index]

      # read the results of initial population
      resInit <- readRDS(paste0(path, i, "_resList_0.rds"))
      trueU <- resInit$trueU
      predU <- resInit$predU
      maxU <- max(trueU)
      accuracy <- cor(trueU, predU)
      gVar <- var(trueU) * (length(trueU) - 1) / length(trueU)
      if (length(index) == 2) {
        resMat0[1, ] <- c(maxU, gVar)
      } else if (length(index) == 3) {
        resMat0[1, ] <- c(maxU, gVar, accuracy)
      }

      for (l in 1:length(index)) {
        if (grepl(pattern = "^top", index[l])) {
          resArray[, l, i] <- (resMat0[, l] - maxU) / sqrt(gVar)
        } else {
          resArray[, l, i] <- resMat0[, l]
        }
      }
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
    return(list(method = method,
                resArray = resArray, resDf = dfEach))
  })

  resAllList <- lapply(resList, function(resListEach) {
    # resListEach <- resList[[1]]
    method <- resListEach$method
    resArrayEach <- resListEach$resArray
    resFin <- resArrayEach[dim(resArrayEach)[1], , ]
    rownames(resFin) <- paste0(method, dimnames(resArrayEach)[[2]])
    return(resFin)
  })

  resAllMat <- do.call(rbind, resAllList)
  for (i in 1:length(index)) {
    # i <- 2
    resIndEach0 <- index[i]
    resIndEach <- paste0(resIndEach0, "$")

    resDfList <- lapply(resList, function(resListEach) {
      # resListEach <- resList[[1]]
      return(resListEach$resDf)
    })
    resDf <- do.call(rbind, resDfList)

    resDfEach <- resDf[resDf$Index == resIndEach0, ]
    resDfEach$Method <- factor(resDfEach$Method, levels = methodFactor)
    col <- methodColor

    comparisonList <- tapply(resDfEach[, "Value"], resDfEach$Method, function(x) {
      resDfEach[resDfEach$Method == "CPS", "Value"] / x
    })
    comparisonMat <- do.call(rbind, comparisonList)
    colnames(comparisonMat) <- unique(resDfEach$Generation)
    write.csv(comparisonMat, paste0(dirSave, target, "_", resIndEach0, "_compare.csv"))

    g <- ggplot(resDfEach, aes(x = Generation, y = Value, color = Method)) +
      geom_line(size = 1.2) +
      theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 24),
            axis.title.y = element_text(size = 6)) +
      scale_color_manual(values = col)

    if (resIndEach0 == "top") {
      g <- g + ylim(c(0, 5.1))
    }

    png(paste0(dirSave, target, "_", resIndEach0, ".png"),
        width = 1440, height = 1440, res = 214)
    print(g)
    dev.off()

    if (resIndEach0 != "variance") {
      simList <- lapply(resList, function(eachResList) {
        # eachResList <- resList[[1]]
        eachResMat <- eachResList$resArray[, resIndEach0, ]
        return(eachResMat)
      })
      simArray <- do.call(abind, c(simList, along = 3))
      topMat <- apply(simArray, c(1, 2), which.max)
      ratioMat <- apply(topMat, 1, function(l) {
        # l <- topMat[2, ]
        ratio <- table(factor(l, levels = 1:length(methodFactor))) / 300
        return(ratio)
      })
      ratio <- c(ratioMat[, 2:ncol(ratioMat)])
      ratioVec <- c(rep(0, length(methodFactor)), ratio)

      generation <- as.numeric(str_sub(rownames(topMat), start = 2))
      generationRep <- rep(generation, each = length(methodFactor))
      methodRep <- rep(methodFactor, length(generation))
      indexRep <- rep(resIndEach0, length(generationRep))
      dfEach <- data.frame(Generation = generationRep,
                           Method = methodRep,
                           Index = indexRep,
                           Value = ratioVec)
      dfEach$Method <- factor(dfEach$Method, levels = methodFactor)
      write.csv(dfEach, paste0(dirSave, target, "_", resIndEach0, "_count.csv"))


      dfRatio <- dfEach[dfEach$Generation != 0, ]
      g <- ggplot(dfRatio, aes(x = Generation, y = Value, fill = Method)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              axis.title.x = element_text(size = 24),
              axis.title.y = element_text(size = 6)) +
        scale_fill_manual(values = col)

      png(paste0(dirSave, target, "_", resIndEach0, "_count.png"),
          width = 1440, height = 1440, res = 214)
      print(g)
      dev.off()
    }
  }
}


SimSumAllele <- function(pathList, generationInd, index, methodFactor, methodColor) {

  resList <- lapply(pathList, function(pathEach) {
    # pathEach <- pathList[[3]]
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
      markerEffectAll <- read.csv(paste0("results/2.1.simulationSetting150/", i, "_markerEffects.csv"),
                                  row.names = 1)
      markerEffectSel <- markerEffectAll[rownames(resMat0), ]
      markerEffectAbs <- abs(markerEffectSel)

      resMat <- apply(resMat0, 2, function(x) {
        fav <- sum(markerEffectAbs[x == 1]) / sum(markerEffectAbs)
        dis <- sum(markerEffectAbs[x == 0]) / sum(markerEffectAbs)
        mid <- 1 - fav - dis
        return(c(Fixed = fav, Lost = dis, NotFixed = mid))
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
      ylim(c(0, 1)) +
      scale_color_manual(values = col)

    png(paste0(dirSave, indexEach, ".png"),
        width = 1440, height = 1440, res = 214)
    print(g)
    dev.off()

  }
}


plotMap <- function (x, map2, chr, horizontal = FALSE, shift = TRUE, show.marker.names = FALSE,
                     alternate.chrid = FALSE, ...)
{
  dots <- list(...)
  if ("main" %in% names(dots)) {
    themain <- dots$main
    usemaindefault <- FALSE
  }
  else usemaindefault <- TRUE
  if ("xlim" %in% names(dots)) {
    xlim <- dots$xlim
    usexlimdefault <- FALSE
  }
  else usexlimdefault <- TRUE
  if ("ylim" %in% names(dots)) {
    ylim <- dots$ylim
    useylimdefault <- FALSE
  }
  else useylimdefault <- TRUE
  if ("xlab" %in% names(dots))
    xlab <- dots$xlab
  else {
    if (horizontal)
      xlab <- "Location (cM)"
    else xlab <- "Chromosome"
  }
  if ("ylab" %in% names(dots))
    ylab <- dots$ylab
  else {
    if (horizontal)
      ylab <- "Chromosome"
    else ylab <- "Location (cM)"
  }
  map <- x
  if (inherits(map, "cross"))
    map <- pull.map(map)
  if (!missing(map2) && inherits(map2, "cross"))
    map2 <- pull.map(map2)
  if (!inherits(map, "map") || (!missing(map2) && !inherits(map2,
                                                            "map")))
    warning("Input should have class \"cross\" or \"map\".")
  if (!missing(map2) && is.matrix(map[[1]]) != is.matrix(map2[[1]]))
    stop("Maps must be both sex-specific or neither sex-specific.")
  if (!missing(chr)) {
    map <- map[matchchr(chr, names(map))]
    if (!missing(map2))
      map2 <- map2[matchchr(chr, names(map2))]
  }
  sex.sp <- FALSE
  if (is.matrix(map[[1]])) {
    one.map <- FALSE
    sex.sp <- TRUE
    if (!missing(map2)) {
      if (is.logical(map2)) {
        horizontal <- map2
        map2 <- lapply(map, function(a) a[2, ])
        map <- lapply(map, function(a) a[1, ])
      }
      else {
        Map1 <- lapply(map, function(a) a[1, , drop = TRUE])
        Map2 <- lapply(map, function(a) a[2, , drop = TRUE])
        Map3 <- lapply(map2, function(a) a[1, , drop = TRUE])
        Map4 <- lapply(map2, function(a) a[2, , drop = TRUE])
        old.mfrow <- par("mfrow")
        on.exit(par(mfrow = old.mfrow))
        par(mfrow = c(2, 1))
        class(Map1) <- class(Map2) <- class(Map3) <- class(Map4) <- "map"
        plotMap(Map1, Map3, horizontal = horizontal,
                shift = shift, show.marker.names = show.marker.names,
                alternate.chrid = alternate.chrid)
        plotMap(Map2, Map4, horizontal = horizontal,
                shift = shift, show.marker.names = show.marker.names,
                alternate.chrid = alternate.chrid)
        return(invisible(NULL))
      }
    }
    else {
      map2 <- lapply(map, function(a) a[2, ])
      map <- lapply(map, function(a) a[1, ])
    }
  }
  else {
    if (!missing(map2))
      one.map <- FALSE
    else one.map <- TRUE
  }
  if (one.map) {
    n.chr <- length(map)
    if (!show.marker.names) {
      chrpos <- 1:n.chr
      thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr * 2, by = 2)
      thelim <- range(chrpos) + c(-0.35, 2.35)
    }
    if (shift)
      map <- lapply(map, function(a) a - a[1])
    maxlen <- max(unlist(lapply(map, max)))
    if (horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- c(0, maxlen)
      if (useylimdefault)
        ylim <- rev(thelim)
      plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
           yaxs = "i", xlab = xlab, ylab = ylab, yaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        segments(min(map[[i]]), chrpos[i], max(map[[i]]),
                 chrpos[i])
        segments(map[[i]], chrpos[i] - 0.25, map[[i]],
                 chrpos[i] + 0.25)
        if (show.marker.names)
          text(map[[i]], chrpos[i] + 0.35, names(map[[i]]),
               srt = 90, adj = c(1, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 2,
                                            at = chrpos[i], labels = names(map)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- thelim
      if (useylimdefault)
        ylim <- c(maxlen, 0)
      plot(0, 0, type = "n", ylim = ylim, xlim = xlim,
           xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        segments(chrpos[i], min(map[[i]]), chrpos[i],
                 max(map[[i]]))
        segments(chrpos[i] - 0.25, map[[i]], chrpos[i] +
                   0.25, map[[i]])
        if (show.marker.names)
          if (any(!is.na(names(map[[i]]))))

            segments(chrpos[i] - 0.25, map[[i]][!(is.na(names(map[[i]])))], chrpos[i] +
                       0.25, map[[i]][!(is.na(names(map[[i]])))], col = "green")
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 1,
                                            at = chrpos[i], labels = names(map)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    if (usemaindefault)
      title(main = "Genetic map")
    else if (themain != "")
      title(main = themain)
  }
  else {
    map1 <- map
    if (is.matrix(map2[[1]]))
      stop("Second map appears to be a sex-specific map.")
    if (length(map1) != length(map2))
      stop("Maps have different numbers of chromosomes.")
    if (any(names(map1) != names(map2))) {
      cat("Map1: ", names(map1), "\n")
      cat("Map2: ", names(map2), "\n")
      stop("Maps have different chromosome names.")
    }
    if (shift) {
      map1 <- lapply(map1, function(a) a - a[1])
      map2 <- lapply(map2, function(a) a - a[1])
    }
    n.mar1 <- sapply(map1, length)
    n.mar2 <- sapply(map2, length)
    markernames1 <- lapply(map1, names)
    markernames2 <- lapply(map2, names)
    if (any(n.mar1 != n.mar2)) {
      if (show.marker.names) {
        warning("Can't show marker names because of different numbers of markers.")
        show.marker.names <- FALSE
      }
    }
    else if (any(unlist(markernames1) != unlist(markernames2))) {
      if (show.marker.names) {
        warning("Can't show marker names because markers in different orders.")
        show.marker.names <- FALSE
      }
    }
    n.chr <- length(map1)
    maxloc <- max(c(unlist(lapply(map1, max)), unlist(lapply(map2,
                                                             max))))
    if (!show.marker.names) {
      chrpos <- 1:n.chr
      thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr * 2, by = 2)
      thelim <- range(chrpos) + c(-0.4, 2.4)
    }
    if (!horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- thelim
      if (useylimdefault)
        ylim <- c(maxloc, 0)
      plot(0, 0, type = "n", ylim = ylim, xlim = xlim,
           xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        if (max(map2[[i]]) < max(map1[[i]]))
          map2[[i]] <- map2[[i]] + (max(map1[[i]]) -
                                      max(map2[[i]]))/2
        else map1[[i]] <- map1[[i]] + (max(map2[[i]]) -
                                         max(map1[[i]]))/2
        segments(chrpos[i] - 0.3, min(map1[[i]]), chrpos[i] -
                   0.3, max(map1[[i]]))
        segments(chrpos[i] + 0.3, min(map2[[i]]), chrpos[i] +
                   0.3, max(map2[[i]]))
        wh <- match(markernames1[[i]], markernames2[[i]])
        for (j in which(!is.na(wh))) segments(chrpos[i] -
                                                0.3, map1[[i]][j], chrpos[i] + 0.3, map2[[i]][wh[j]])
        if (any(is.na(wh)))
          segments(chrpos[i] - 0.4, map1[[i]][is.na(wh)],
                   chrpos[i] - 0.2, map1[[i]][is.na(wh)])
        wh <- match(markernames2[[i]], markernames1[[i]])
        if (any(is.na(wh)))
          segments(chrpos[i] + 0.4, map2[[i]][is.na(wh)],
                   chrpos[i] + 0.2, map2[[i]][is.na(wh)])
        if (show.marker.names)
          text(chrpos[i] + 0.35, map2[[i]], names(map2[[i]]),
               adj = c(0, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 1,
                                            at = chrpos[i], labels = names(map1)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map1)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map1)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- c(0, maxloc)
      if (useylimdefault)
        ylim <- rev(thelim)
      plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab, yaxt = "n", yaxs = "i")
      a <- par("usr")
      for (i in 1:n.chr) {
        if (max(map2[[i]]) < max(map1[[i]]))
          map2[[i]] <- map2[[i]] + (max(map1[[i]]) -
                                      max(map2[[i]]))/2
        else map1[[i]] <- map1[[i]] + (max(map2[[i]]) -
                                         max(map1[[i]]))/2
        segments(min(map1[[i]]), chrpos[i] - 0.3, max(map1[[i]]),
                 chrpos[[i]] - 0.3)
        segments(min(map2[[i]]), chrpos[i] + 0.3, max(map2[[i]]),
                 chrpos[[i]] + 0.3)
        wh <- match(markernames1[[i]], markernames2[[i]])
        for (j in which(!is.na(wh))) segments(map1[[i]][j],
                                              chrpos[i] - 0.3, map2[[i]][wh[j]], chrpos[i] +
                                                0.3)
        if (any(is.na(wh)))
          segments(map1[[i]][is.na(wh)], chrpos[i] -
                     0.4, map1[[i]][is.na(wh)], chrpos[i] - 0.2)
        wh <- match(markernames2[[i]], markernames1[[i]])
        if (any(is.na(wh)))
          segments(map2[[i]][is.na(wh)], chrpos[i] +
                     0.4, map2[[i]][is.na(wh)], chrpos[i] + 0.2)
        if (show.marker.names)
          text(map2[[i]], chrpos[i] + 0.35, names(map2[[i]]),
               srt = 90, adj = c(1, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 2,
                                            at = chrpos[i], labels = names(map1)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map1)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map1)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    if (usemaindefault) {
      if (!sex.sp)
        title(main = "Comparison of genetic maps")
      else title(main = "Genetic map")
    }
    else if (themain != "")
      title(main = themain)
  }
  invisible()
}
