##############################################
# Title : 2.0.simulationSetting198
# Author : Kengo Sakurai
# Date : 2023-07-21
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(breedSimulatR)
library(stringr)
library(gaston)
library(lme4)
library(RAINBOWR)
library(qtl)
source(file = "R/1.0.function.R")

# 1.1. Parameters
threshold <- 0.6 # for maximum LD
nSNP <- 200 # the number of SNPs in each chr
nChr <- 20 # the number of chromosomes
lchr <- 1e9 # physical length of each chr (no meaning)
lchrCm <- 100 # map length of each chr
nQtn <- 1000 # the number of QTN
nTrait <- 1

# 1.2. Save dir
dirSave <- "results/2.0.simulationSetting198/"
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

# 1.3. Read data
genoInfoPath <- "raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz"
genoInfoRaw0 <- read.vcf(genoInfoPath)

genoInfoRaw0@snps$id <- paste0("Chr",
                               formatC(genoInfoRaw0@snps$chr, width = 2, flag = "0"),
                               "_",
                               formatC(genoInfoRaw0@snps$pos, width = 8, flag = "0"))
genomeFil <- select.snps(genoInfoRaw0, condition = maf >= 0.1)
genomeFil <- LD.thin(genomeFil, threshold = threshold)

genomeMat <- as.matrix(genomeFil)
snpInfoList <- str_split(colnames(genomeMat), pattern = "_")
snpInfoMat <- do.call(rbind, snpInfoList)
snpInfoDf0 <- data.frame(chr = as.character(snpInfoMat[, 1]),
                         SNPid = colnames(genomeMat),
                         physPos = as.numeric(snpInfoMat[, 2]))
snpInfoDf0 <- snpInfoDf0[order(snpInfoDf0$SNPid), ]

# calculating map position
mapPosList <- tapply(X = snpInfoDf0$physPos, INDEX = snpInfoDf0$chr, function(eachChr) {
  # eachChr <- snpInfoDf0$PysPos[snpInfoDf0$Chr == snpInfoDf0$Chr[1]]
  pysPosMax <- max(eachChr)
  pysPosMin <- min(eachChr)
  mapPos <- (eachChr -  pysPosMin) * 100 / (pysPosMax - pysPosMin)
  return(mapPos)
})
mapPos <- unlist(mapPosList)
snpInfoDf <- data.frame(snpInfoDf0,
                        linkMapPos = mapPos)

for (i in 1:300) {
  # i <- 1
  # sampling the SNPs (totally nSNPs)
  set.seed(seedInd[i])
  snpSel <- unlist(tapply(snpInfoDf$SNPid, snpInfoDf$chr, sample, nSNP))
  snpSel <- sort(snpSel)
  snpCoord <- snpInfoDf[snpInfoDf$SNPid %in% snpSel, ]
  rownames(snpCoord) <- snpCoord$SNPid
  genomeMatSel <- genomeMat[, snpSel]

  # #### Creating the initial population ####
  # create specie object
  # We set "lchr", but it has no meaning
  specie_statEx <- specie$new(specName = "Soybean",
                              nChr = nChr,
                              lchr = lchr,
                              lchrCm = lchrCm)

  # # create SNPinfo object
  SNPs <- SNPinfo$new(SNPcoord = snpCoord,
                      specie = specie_statEx)
  map <- SNPs$SNPcoord
  map <- map[order(map$SNPid), ]
  write.csv(map, file = paste0(dirSave, i, "_map.csv"))

  # create population object
  initPop <- createPop(geno = genomeMatSel,
                       SNPinfo = SNPs,
                       popName = "Initial population")
  saveRDS(initPop, file = paste0(dirSave, i, "_initPop.rds"))

  # set the marker effect
  set.seed(seedInd[i])
  qtn <- sample(names(initPop$maf > 0.1), nQtn)
  set.seed(seedInd[i])
  yield <- trait$new(name = "Yield",
                     qtn = qtn,
                     qtnEff = rnorm(nQtn, 0, 0.35))

  # save the marker effect
  markerEff <- rep(0, length(map$SNPid))
  names(markerEff) <- map$SNPid
  markerEff[yield$qtn] <- yield$qtnEff
  write.csv(markerEff, file = paste0(dirSave, i, "_markerEffects.csv"))

  ##### heritability = 1 #####
  phenolab <- phenotyper$new(name = "Pheno lab",
                             traits = yield,
                             plotCost = 150,
                             mu = 0,
                             he = 1,
                             pop = initPop)
  saveRDS(phenolab, file = paste0(dirSave, i, "_phenolab.rds"))

  # save a genome data and phenotype data for visualizing the qtn data
  genoMat <- initPop$genoMat
  genoMatSel <- genoMat[, qtn]
  genoMatSel[genoMatSel == 0] <- "H"
  genoMatSel[genoMatSel == 2] <- "A"
  mapSel <- map[qtn, ]
  chr <- as.numeric(str_sub(mapSel$chr, 4, 5))
  genoQTL <- data.frame(chr, mapSel$linkMapPos, t(genoMatSel))
  genoDf <- data.frame(id = c("", "", rownames(genoMatSel)), t(genoQTL))
  write.csv(genoDf,
            file = paste0(dirSave, i, "_genoQTL.csv"), row.names = F)

  pheno0 <- phenolab$trial(pop = initPop, rep = 1)$data
  pheno <- pheno0[, 1:2, drop = F]
  colnames(pheno) <- c("id", "trait1")
  rownames(pheno) <- pheno0$ind
  write.csv(pheno, file = paste0(dirSave, i, "_phenoQTL.csv"))

  # read data for visualizing the location of QTNs
  cross <- read.cross(format = "csvs",
                      genfile = paste0(dirSave, i, "_genoQTL.csv"),
                      phefile = paste0(dirSave, i, "_phenoQTL.csv"))

  # visualize the location of QTNs
  png(filename = paste0(dirSave, i, "_QTNs.png"))
  plotMap(cross)
  dev.off()

}
