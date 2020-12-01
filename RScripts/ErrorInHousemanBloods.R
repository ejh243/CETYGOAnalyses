## Script to test error metric using Houseman data set
## Dorothea Seiler Vellame
## started 18-11-2020

### load data and assign training and testing #########
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

compositeCellType = "Blood"
platform<-"450k"
referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
data(list = referencePkg)
referenceRGset <- get(referencePkg)
pheno = pData(referenceRGset)$CellType

## only keep the 6 wanted cell types
index = which(pheno == "Bcell" | 
                pheno == "CD4T" | 
                pheno == "CD8T" | 
                pheno == "Gran" | 
                pheno == "Mono" | 
                pheno == "NK")
pheno = as.factor(pheno[index]) 
betas = referenceRGset[,index]

## subset betas for those in EPIC only


# betas = getBeta(betas)

x = read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7, header = T)
EPIC450kLoci = getAnnotation(betas)
y = subsetByLoci(betas, includeLoci = rownames(EPIC450kLoci)[rownames(EPIC450kLoci) %in% x$IlmnID])
betas = y



## split into training and testing and add to phenotype cols
phenoTrain = c(rep("Test", length(pheno)))
for(i in 1:length(levels(pheno))){
  cellTypeIndex = which(pheno == levels(pheno)[i])
  set.seed(123)
  phenoTrain[sample(cellTypeIndex, ceiling(length(cellTypeIndex)/2))] = "Train"
}


## save as RGset
betasTrain = betas[,phenoTrain == "Train"]
betasTest = betas[,phenoTrain == "Test"]
pheno = data.frame(celltype = pheno, trainTest = phenoTrain)

save(betas, betasTrain, betasTest, pheno, file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestRGSet.Rdata")

## save as matrix
betas = getBeta(preprocessRaw(betas))
betasTrain = getBeta(preprocessRaw(betasTrain))
betasTest = getBeta(preprocessRaw(betasTest))
save(betas, betasTrain, betasTest, pheno, file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## plot PCA of betas per cell type coloured by train test
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")
library(ggfortify)
library(ggplot2)
library(cowplot)
library(scales)

plotDat = list()
colours = hue_pal()(6)
shapes = c(16, 17, 15, 3, 7, 8)
for (i in 1:length(levels(pheno$celltype))){
  dat = betas[, pheno$celltype == levels(pheno$celltype)[i]]
  datPheno = pheno[pheno$celltype == levels(pheno$celltype)[i],]
  betaVar = apply(dat, 1, var, na.rm = T)
  topBeta = dat[order(betaVar, decreasing = T)[1:1000],]
  plot = prcomp(t(topBeta)) 
  plotDat[[i]] = autoplot(plot, col = colours[i], data = datPheno, shape = "trainTest", size = 2) + 
    theme_cowplot(18) + 
    ggtitle(levels(pheno$celltype)[i]) +
    theme(legend.position = "none")
}

# betaVar = apply(betas, 1, var, na.rm = T)
# topBeta = betas[order(betaVar, decreasing = T)[1:1000],]
# plotDat[[7]] = autoplot(prcomp(t(topBeta)), data = pheno, col = "celltype", shape = "trainTest", size = 2) +
#   theme_cowplot(18) 

cellPlots = plot_grid(plotDat[[1]],
          plotDat[[2]],
          plotDat[[3]],
          plotDat[[4]],
          plotDat[[5]],
          plotDat[[6]],
          ncol = 3, labels = "AUTO")

# fullCellPCA = plotDat[[7]] +
#   labs(shape = "", col = "Cell type")

## quantile normalise together and apart
quantileBetasTrain = getBeta(preprocessQuantile(betasTrain, fixOutliers = TRUE,
                                                removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                                quantileNormalize = TRUE, stratified = TRUE, 
                                                mergeManifest = FALSE, sex = NULL))
quantileBetasTest = getBeta(preprocessQuantile(betasTest, fixOutliers = TRUE,
                                               removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                               quantileNormalize = TRUE, stratified = TRUE, 
                                               mergeManifest = FALSE, sex = NULL))
quantileBetas = getBeta(preprocessQuantile(betas, fixOutliers = TRUE,
                                               removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                               quantileNormalize = TRUE, stratified = TRUE, 
                                               mergeManifest = FALSE, sex = NULL))


phenoTrain = pheno[pheno$trainTest == "Train",]
phenoTest = pheno[pheno$trainTest == "Test",]

## save data
save(quantileBetasTrain, quantileBetasTest, quantileBetas, phenoTest, phenoTrain,
     file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")



### Optimise number of CpGs ###########################
## load normalised data 
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")

## source model functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

# ## make models
# modelListCpG = list()
# cpg = c(1:25,seq(30,55,5), seq(60, 200, 10))
# for (i in 1:length(cpg)) {
#   modelListCpG[[i]] = pickCompProbes(rawbetas = as.matrix(quantileBetasTrain),
#                                      cellTypes = levels(as.factor(phenoTrain$celltype)),
#                                      cellInd = as.factor(phenoTrain$celltype),
#                                      numProbes =  cpg[i],
#                                      probeSelect = "auto")
#   names(modelListCpG)[i] = paste("model", cpg[i], "CG", sep = "")
# }
# 
# ## save model list
# save(modelListCpG, file = "/mnt/data1/Thea/ErrorMetric/data/nCpGModels/nCpGModels.Rdata")

load("/mnt/data1/Thea/ErrorMetric/data/nCpGModels/nCpGModels.Rdata")
library(reshape2)

## predict value of test data for each model
cpg = c(1:25,seq(30,55,5), seq(60, 200, 10))
trueProp = singleCellProportionMatrix(phenoTest[,1])
rmseTvP = c()

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

for (i in 1:length(cpg)) {
  pred = projectCellTypeWithError(quantileBetasTest,  modelType = "ownModel", ownModelData = modelListCpG[[i]])
  rmseTvP = c(rmseTvP, RMSE(melt(trueProp)[,3], melt(pred[,-which(colnames(pred) %in% c("nCGmissing", "error"))])[,3]))
}


corPlot = data.frame(rmseTvP, cpg)

ggplot(corPlot, aes(x = cpg, y = rmseTvP)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Number of CpGs", y = "RMSE")


## save model that used 150 CpGs
HousemanBlood150CpGModel = modelListCpG[[which(cpg == 150)]]
save(HousemanBlood150CpGModel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")



### Plot heirarchical cluster of training cpgs used in model #####
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")

## load data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")

## load functions
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

modelBetas = GetModelCG(quantileBetasTrain, list(HousemanBlood150CpGModel))
library(gplots)
library(viridis)
library(scales)
library(ComplexHeatmap)

# colours6 = hue_pal()(6)

col = list(Celltype = c("Bcell" = "#F8766D", "CD4T" = "#B79F00",
                        "CD8T" = "#00BA38",  "Gran" = "#00BFC4",
                        "Mono" = "#619CFF",  "NK" = "#F564E3"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(Celltype = phenoTrain$celltype,
  col = col)

Heatmap(modelBetas, name = "DNAm",
          top_annotation = ha, show_row_names = F, show_column_names = F)


### check effect of normalisation #####################
## load normalised data 
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## source model functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## create model in separately normalised train data
sepNormalisedModel = pickCompProbes(rawbetas = quantileBetasTrain,
                                    cellTypes = levels(as.factor(phenoTrain$celltype)),
                                    cellInd = as.factor(phenoTrain$celltype),
                                    numProbes =  150,
                                    probeSelect = "auto")

## create model in combined normalised train data
combNormalisedModel = pickCompProbes(rawbetas = quantileBetas[, pheno$trainTest == "Train"],
                                     cellTypes = levels(as.factor(phenoTrain$celltype)),
                                     cellInd = as.factor(phenoTrain$celltype),
                                     numProbes =  150,
                                     probeSelect = "auto")

## apply both models to their respective test data sets
sepNormalisedPrediction = projectCellTypeWithError(quantileBetasTest, model = "ownModel", ownModelData = sepNormalisedModel)
combNormalisedPrediction = projectCellTypeWithError(quantileBetas[, pheno$trainTest == "Test"], model = "ownModel", ownModelData = combNormalisedModel)
trueProp = singleCellProportionMatrix(phenoTest[,1])

library(ggplot2)
library(cowplot)
library(reshape2)

# get absolute difference in predictions
sepMelt = melt(sepNormalisedPrediction[,-which(colnames(sepNormalisedPrediction) %in% c("nCGmissing", "error"))])
combMelt = melt(combNormalisedPrediction[,-which(colnames(combNormalisedPrediction) %in% c("nCGmissing", "error"))])
absDiffSep = abs(sepMelt[,3] - melt(trueProp)[,3])
absDiffComb = abs(combMelt[,3] - melt(trueProp)[,3])

plotDat = data.frame(diff = c(absDiffComb, absDiffSep), 
                     norm = rep(c("Seperate", "Combined"), each = length(absDiffComb)),
                     cell = as.factor(c(as.character(sepMelt[,2]), as.character(combMelt[,2]))))

## use a paired one-sided t-test to compare the groups
ttest = t.test(absDiffComb, absDiffSep, paired = T, alternative = "less")

ggplot(plotDat, aes(x = norm, y = diff)) +
  geom_violin() +
  geom_jitter(aes(col = cell, shape = cell)) +
  theme_cowplot(18) +
  labs(x = "Normalization", y = "Absolute difference between\ntrue and predicted") +
  theme(legend.title = element_blank()) +
  ylim(c(-0.001,0.3)) +
  annotate("text", x = 1.5, y = 0.3, label = paste("p =", signif(ttest$p.value, 3)), size = 5)



### create models using 3:6 cell types
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")

## make design matrix of booleans for which cell type will be present
designMatrix = expand.grid(c(T,F), c(T,F), c(T,F), c(T,F), c(T,F), c(T,F))
designMatrix = designMatrix[apply(designMatrix, 1, sum) >= 3,]

cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")
cellTypeShorthand = c("B", "C4", "C8", "G", "M", "NK")

modelList = list()
for (i in 1:nrow(designMatrix)){
  modellingDat = CellTypeSubsetBetasAndPheno(cellTypes[unlist(designMatrix[i,])],
                                             quantileBetasTrain, phenoTrain, phenoColName = "celltype", justBetas = F)
  modelList[[i]] = pickCompProbes(rawbetas = as.matrix(modellingDat[[1]]),
                                  cellTypes = levels(as.factor(as.character(modellingDat[[2]]$celltype))),
                                  cellInd = as.factor(as.character(modellingDat[[2]]$celltype)),
                                  numProbes =  150,
                                  probeSelect = "auto")
  names(modelList)[i] = paste("model", paste(cellTypeShorthand[unlist(designMatrix[i,])], sep = "", collapse = ""), sep = "_")
}

save(modelList, file = "/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")





bulk = CellTypeProportionSimulator(GetModelCG(quantileBetasTest, modelList),
                            phenoTest,
                            phenoColName = "celltype",
                            nBulk = 15,
                            proportionsMatrixType = "random")

stackedPlots = ModelCompareStackedBar(bulk[[1]],
                           modelList,
                           nCpGPlot = F,
                           sampleNamesOnPlots = F,
                           trueComparison = T,
                           trueProportions = bulk[[2]])

save(modelList, bulk, stackedPlots, file = "/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")









































## check effect of age, sex in EXTEND and US #########
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")
model = sepNormalisedModel

## load functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")

## load data
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
us = dat
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batches_1_2_normalised_together.rdat")
ex = betas

rm(dat, betas, pheno, sepNormalisedModel)

## get error for US
usPred = projectCellTypeWithError(us, modelType = "ownModel", ownModelData = model)

## get error for EX
exPred = projectCellTypeWithError(ex, modelType = "ownModel", ownModelData = model)