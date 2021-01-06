## Script to test error metric using Houseman data set
## Dorothea Seiler Vellame
## started 18-11-2020

### Load data and assign training and testing #########
library(minfi)
library(wateRmelon)
library(FlowSorted.Blood.450k)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

compositeCellType = "Blood"
platform<-"450k"
referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
data(list = referencePkg)
referenceRGset <- get(referencePkg)
phenoDat = pData(referenceRGset)$CellType

## only keep the 6 wanted cell types
index = which(phenoDat == "Bcell" | 
                phenoDat == "CD4T" | 
                phenoDat == "CD8T" | 
                phenoDat == "Gran" | 
                phenoDat == "Mono" | 
                phenoDat == "NK")
phenoDat = as.factor(phenoDat[index]) 
betas = referenceRGset[,index]

## subset betas for those in EPIC only
x = read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7, header = T)
EPIC450kLoci = getAnnotation(betas)
y = subsetByLoci(betas, includeLoci = rownames(EPIC450kLoci)[rownames(EPIC450kLoci) %in% x$IlmnID])


## subset for those in EXTEND and Understanding Society
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
us = dat
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batches_1_2_normalised_together.rdat")
ex = betas
rm(pheno)

USEXrows = rownames(us)[rownames(us) %in% rownames(ex)]
EPIC450kLoci = getAnnotation(y)
y = subsetByLoci(y, includeLoci = rownames(EPIC450kLoci)[rownames(EPIC450kLoci) %in% USEXrows])
betas = y

## split into training and testing and add to phenotype cols
phenoTrain = c(rep("Test", length(phenoDat)))
for(i in 1:length(levels(phenoDat))){
  cellTypeIndex = which(phenoDat == levels(phenoDat)[i])
  set.seed(123)
  phenoTrain[sample(cellTypeIndex, ceiling(length(cellTypeIndex)/2))] = "Train"
}


## save as RGset
betasTrain = betas[,phenoTrain == "Train"]
betasTest = betas[,phenoTrain == "Test"]
pheno = data.frame(celltype = phenoDat, trainTest = phenoTrain)

phenoTrain = pheno[pheno$trainTest == "Train",]
phenoTest = pheno[pheno$trainTest == "Test",]

save(betas, betasTrain, betasTest, pheno, file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestRGSet.Rdata")

## save as matrix
betas = getBeta(preprocessRaw(betas))
betasTrain = getBeta(preprocessRaw(betasTrain))
betasTest = getBeta(preprocessRaw(betasTest))
save(betas, betasTrain, betasTest, pheno, phenoTrain, phenoTest,
     file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")


## quantile normalise together and apart
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestRGSet.Rdata")
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




## save data
save(quantileBetasTrain, quantileBetasTest, quantileBetas, phenoTest, phenoTrain, pheno,
     file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")

## plot PCA of betas per cell type coloured by train test
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")
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


pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/cellTrainTestPCA.pdf", height = 7, width = 10)
plot_grid(plotDat[[1]],
          plotDat[[2]],
          plotDat[[3]],
          plotDat[[4]],
          plotDat[[5]],
          plotDat[[6]],
          ncol = 3, labels = "AUTO")
dev.off()

## full PCA to show celltype similarity
betaVar = apply(betas, 1, var, na.rm = T)
topBeta = betas[order(betaVar, decreasing = T)[1:1000],]
plotDatPCAFULL = autoplot(prcomp(t(topBeta)), data = pheno, col = "celltype", shape = "celltype", size = 2) +
  theme_cowplot(18)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/cellPCA.pdf",height = 7, width = 7)
plotDatPCAFULL +
  labs(col = "Cell type", shape = "Cell type")
dev.off()


### Check effect of normalisation #####################
## load normalised data 
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## source model functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## create model in separately normalised train data
unnormalisedModel = pickCompProbes(rawbetas = betasTrain,
                                   cellTypes = levels(as.factor(phenoTrain$celltype)),
                                   cellInd = as.factor(phenoTrain$celltype),
                                   numProbes =  150,
                                   probeSelect = "auto")

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

## apply models to their respective test data sets
unnormalisedPrediction = projectCellTypeWithError(betasTest, model = "ownModel", ownModelData = unnormalisedModel)
sepNormalisedPrediction = projectCellTypeWithError(quantileBetasTest, model = "ownModel", ownModelData = sepNormalisedModel)
combNormalisedPrediction = projectCellTypeWithError(quantileBetas[, pheno$trainTest == "Test"], model = "ownModel", ownModelData = combNormalisedModel)
unModNormDatPrediction = projectCellTypeWithError(quantileBetasTest, model = "ownModel", ownModelData = unnormalisedModel)
normModUnDatPrediction = projectCellTypeWithError(betasTest, model = "ownModel", ownModelData = sepNormalisedModel)
trueProp = singleCellProportionMatrix(phenoTest[,1])

library(ggplot2)
library(cowplot)
library(reshape2)

# get absolute difference in predictions
unnormMelt = melt(unnormalisedPrediction[,-which(colnames(unnormalisedPrediction) %in% c("nCGmissing", "error"))])
sepMelt = melt(sepNormalisedPrediction[,-which(colnames(sepNormalisedPrediction) %in% c("nCGmissing", "error"))])
combMelt = melt(combNormalisedPrediction[,-which(colnames(combNormalisedPrediction) %in% c("nCGmissing", "error"))])
unModNormDatMelt = melt(unModNormDatPrediction[,-which(colnames(unModNormDatPrediction) %in% c("nCGmissing", "error"))])
normModUnDatMelt = melt(normModUnDatPrediction[,-which(colnames(normModUnDatPrediction) %in% c("nCGmissing", "error"))])

absDiffUn = abs(unnormMelt[,3] - melt(trueProp)[,3])
absDiffSep = abs(sepMelt[,3] - melt(trueProp)[,3])
absDiffComb = abs(combMelt[,3] - melt(trueProp)[,3])
absDiffunModNormDat = abs(unModNormDatMelt[,3] - melt(trueProp)[,3])
absDiffnormModUnDat = abs(normModUnDatMelt[,3] - melt(trueProp)[,3])


plotDat = data.frame(diff = c(absDiffUn, absDiffComb, absDiffSep, absDiffunModNormDat, absDiffnormModUnDat), 
                     norm = rep(c("Unnormalised", "Seperate", "Combined", "Unnormalised Train\nNormalised Test", "Normalised Train\nUnnormalised Test"), each = length(absDiffComb)),
                     cell = as.factor(c(as.character(unnormMelt[,2]), as.character(sepMelt[,2]), as.character(combMelt[,2]), as.character(sepMelt[,2]), as.character(sepMelt[,2]))))

## use a paired t-test to compare each group to comb, the default
dat = cbind(absDiffUn, absDiffComb, absDiffSep, absDiffunModNormDat, absDiffnormModUnDat)
library(xtable)
xtable(rbind.data.frame(t.testMatrix(dat), mean = signif(colMeans(dat),3), 
                        SD = signif(apply(dat, 2, sd),3)))

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/sepCombNormalisationComparison.pdf", 
    height = 7, width = 11)
ggplot(plotDat, aes(x = norm, y = diff)) +
  geom_violin() +
  geom_jitter(aes(col = cell, shape = cell)) +
  theme_cowplot(18) +
  labs(x = "Normalisation", y = "Absolute difference between\ntrue and predicted") +
  theme(legend.title = element_blank())
dev.off()



### Optimise number of CpGs ###########################
## load normalised data 
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## source model functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

# ## make models
# modelListCpG = list()
# cpg = c(1:25,seq(30,55,5), seq(60, 200, 10))
# for (i in 1:length(cpg)) {
#   modelListCpG[[i]] = pickCompProbes(rawbetas = as.matrix(betasTrain),
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
  pred = projectCellTypeWithError(betasTest,  modelType = "ownModel", ownModelData = modelListCpG[[i]])
  rmseTvP = c(rmseTvP, RMSE(melt(trueProp)[,3], melt(pred[,-which(colnames(pred) %in% c("nCGmissing", "error"))])[,3]))
}


corPlot = data.frame(rmseTvP, cpg)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/nCpGNeededForModel.pdf", height = 7, width = 7)
ggplot(corPlot, aes(x = cpg, y = rmseTvP)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Number of CpGs", y = "RMSE")
dev.off()

## save model that used 150 CpGs
HousemanBlood150CpGModel = modelListCpG[[which(cpg == 150)]]
save(HousemanBlood150CpGModel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")



### Plot heirarchical cluster of training cpgs used in model #####
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")

## load data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## load functions
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

modelBetas = GetModelCG(betasTrain, list(HousemanBlood150CpGModel))
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

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/heatmapForModelCpGs.pdf", height = 8, width = 7)
Heatmap(modelBetas, name = "DNAm",
        top_annotation = ha, show_row_names = F, show_column_names = F, show_row_dend = F)
dev.off()



### Create models using 3:6 cell types ################
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## make design matrix of booleans for which cell type will be present
designMatrix = expand.grid(c(T,F), c(T,F), c(T,F), c(T,F), c(T,F), c(T,F))
designMatrix = designMatrix[apply(designMatrix, 1, sum) >= 3,]

cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")
cellTypeShorthand = c("B", "C4", "C8", "G", "M", "NK")

# modelList = list()
# for (i in 1:nrow(designMatrix)){
#   modellingDat = CellTypeSubsetBetasAndPheno(cellTypes[unlist(designMatrix[i,])],
#                                              betasTrain, phenoTrain, phenoColName = "celltype", justBetas = F)
#   modelList[[i]] = pickCompProbes(rawbetas = as.matrix(modellingDat[[1]]),
#                                   cellTypes = levels(as.factor(as.character(modellingDat[[2]]$celltype))),
#                                   cellInd = as.factor(as.character(modellingDat[[2]]$celltype)),
#                                   numProbes =  150,
#                                   probeSelect = "auto")
#   names(modelList)[i] = paste("model", paste(cellTypeShorthand[unlist(designMatrix[i,])], sep = "", collapse = ""), sep = "_")
# }
# 
# save(modelList, file = "/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")
# 
# load("/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")
# 
# ## mean proportions of each cell type in whole blood (from Reinius2012)
# meanBloodProp = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43)
# sdBloodProp = c(1.44, 3.12, 3.13, 9.19, 3.17, 1.5)
# 
# # each cell type will be simulated with mean, mean +- sd, mean +- 2sd
# 
# bulk = CellTypeProportionSimulator(GetModelCG(betasTest, modelList),
#                             phenoTest,
#                             phenoColName = "celltype",
#                             nBulk = 30,
#                             proportionsMatrixType = "own",
#                             proportionsMatrix = simPropMaker(meanBloodProp, sdBloodProp))
# 
# stackedPlots = ModelCompareStackedBar(bulk[[1]],
#                            modelList,
#                            nCpGPlot = F,
#                            sampleNamesOnPlots = F,
#                            trueComparison = T,
#                            trueProportions = bulk[[2]])
# 
# save(modelList, bulk, stackedPlots, file = "/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")


load("/mnt/data1/Thea/ErrorMetric/data/nCellTypeModels/VaryNCellsData.Rdata")

x = stackedPlots[[1]]$data

## add columns for which cell types in each data set and then compare specific models
modelPresent = matrix(ncol = length(cellTypes), nrow = nrow(x), data = 0)
colnames(modelPresent) = paste(cellTypeShorthand, "cell")


for(i in 1:6){
  modelPresent[grep(cellTypeShorthand[i],x$model),i] = 1
}

y = rowSums(modelPresent)
plotDat = data.frame(x,modelPresent, sums = y)

## colour by number of cells in model
plotDatBox = spread(plotDat, key = c(cellType), value = c(proportion_pred))

## no stats comparison for any with simulated data as n is decided by me

pdf("/mnt/data1/Thea/ErrorMetric/plots/badModels/violinnCelltypeModels.pdf", height = 6, width = 6)
ggplot(plotDatBox, aes(x = as.factor(sums), y = error, fill = as.factor(sums))) +
  geom_violin() +
  theme_cowplot(18) +
  labs(x = "Number of cell types in the model", y = "DSRMSE") +
  ylim(c(0, max(plotDatBox$error))) +
  theme(legend.position = "none") 
dev.off() 


## general Q: 
## why do some still predict well? low proportion of that cell type? cell type less important?

## compare those with only 5 to simplify the question
plotDat5 = plotDatBox[plotDatBox$sums == 5, ]

model5Index = which(rowSums(designMatrix) == 5)
models5 = modelList[model5Index]

wantedModelNames = sapply(strsplit(names(models5), "_"), function(x){return(x[[2]])})
names(models5) = wantedModelNames

## for each cell type, plot stacked bar of true and actual and error for their own simulated data 
meanBloodProp = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43)
models5Compared = list()
for (i in 1:6){
  celltypePreBool = 1:6
  bulkProp = simPropMaker2(meanBloodProp, celltypeBool = celltypePreBool == i, cellNames = cellTypes)
  bulk = CellTypeProportionSimulator(GetModelCG(betasTest, list(models5[[i]])),
                                     phenoTest,
                                     phenoColName = "celltype",
                                     nBulk = nrow(bulkProp),
                                     proportionsMatrixType = "own",
                                     proportionsMatrix = bulkProp)
  model = list(models5[[i]])
  names(model) = wantedModelNames[i]
  models5Compared[[i]] = ModelCompareStackedBar(bulk[[1]], 
                                                modelList = model, 
                                                trueComparison = T,
                                                noise = F,
                                                trueProportions = bulk[[2]],
                                                nCpGPlot = F,
                                                sampleNamesOnPlots = F)
}

for(i in 1:6){
  leg = get_legend(models5Compared[[i]][[3]]
                   + theme(legend.position=c(0.05,0.8),legend.direction = "horizontal", legend.title = element_blank())
                   + guides(fill = guide_legend(nrow = 1)))
  plots = plot_grid(models5Compared[[i]][[1]] + theme(legend.position = "none"),
                    models5Compared[[i]][[2]] + theme(legend.position = "none"),
                    models5Compared[[i]][[3]] + theme(legend.position = "none"), ncol = 1,
                    rel_heights = c(0.6,1,1), labels = "AUTO", axis = "rl", align = "v" )
  pdf(paste("/mnt/data1/Thea/ErrorMetric/plots/badModels/stackedBar5cell",cellTypes[i],"missing.pdf", sep = ""), height = 9, width = 7)     
  print(plot_grid(plots, leg, ncol = 1, rel_heights = c(1,0.08)))
  dev.off()
}


## plot error metrics for each 10 samples together to highlight the difference in magnitude
erPlotDat = models5Compared[[1]][[1]]$data
for (i in 2:6){
  erPlotDat = rbind.data.frame(erPlotDat, models5Compared[[i]][[1]]$data)
}

pdf("/mnt/data1/Thea/ErrorMetric/plots/badModels/errorOnly5Celltypes.pdf", height = 4, width = 7)
ggplot(erPlotDat, aes(x = sample, y = error, col = model)) +
  geom_point() +
  theme_cowplot(18) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "Sample", y = "DSRMSE", col = "Model") +
  ylim(c(0, max(erPlotDat$error)))
dev.off()





### Compare error when noise is introduced ############

## load functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")

## load testing data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## mean proportions of each cell type in whole blood (from Reinius2012)
meanBloodProp = matrix(nrow = 1, byrow = T, data = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43))
colnames(meanBloodProp) = levels(phenoTest$celltype)

noise = seq(0,0.95,0.05)

## create simulated samples with increasing noise
testData = CellTypeProportionSimulator(betas = betasTest, 
                                       pheno = phenoTest, 
                                       phenoColName = "celltype", 
                                       nBulk = length(noise), 
                                       proportionsMatrixType = "own",
                                       proportionsMatrix = meanBloodProp,
                                       noiseIn = T,
                                       proportionNoise = noise)


stackedWithNoise = ModelCompareStackedBar(testBetas = testData[[1]], 
                                          modelList = list(Predicted = HousemanBlood150CpGModel), 
                                          trueComparison = T,
                                          noise = T,
                                          trueProportions = testData[[2]],
                                          nCpGPlot = F,
                                          sampleNamesOnPlots = F)

leg = get_legend(stackedWithNoise[[3]]
                 + theme(legend.position=c(0.03,0.8),legend.direction = "horizontal", legend.title = element_blank())
                 + guides(fill = guide_legend(nrow = 1)))
plots = plot_grid(stackedWithNoise[[1]] + theme(legend.position = "none"),
                  stackedWithNoise[[2]] + theme(legend.position = "none"),
                  stackedWithNoise[[3]] + theme(legend.position = "none"), ncol = 1,
                  rel_heights = c(0.6,1,1), labels = "AUTO", axis = "rl", align = "v" )

pdf("/mnt/data1/Thea/ErrorMetric/plots/badData/simWithNoise.pdf", height = 9, width = 7.5)     
print(plot_grid(plots, leg, ncol = 1, rel_heights = c(1,0.08)))
dev.off()



### Check increasing missingness of CpGs ##############
## load functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")

## load testing data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## mean proportions of each cell type in whole blood (from Reinius2012)
meanBloodProp = matrix(nrow = 1, byrow = T, data = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43))/100
colnames(meanBloodProp) = levels(phenoTest$celltype)

## create a single representative sample
testData = CellTypeProportionSimulator(betas = GetModelCG(betasTest, list(HousemanBlood150CpGModel)), 
                                       pheno = phenoTest, 
                                       phenoColName = "celltype", 
                                       nBulk = 1, 
                                       proportionsMatrixType = "own",
                                       proportionsMatrix = meanBloodProp,
                                       noiseIn = F)


## create function to add x NAs to betas
MakeNAsInBetas = function(proportionNA, betas){
  nToBeNA = floor(nrow(betas)*proportionNA)
  NAIndex = sample(nrow(betas), nToBeNA)
  betas[NAIndex,] = NA
  return(betas)
}

propMissing = seq(0,0.9,0.05)
x = lapply(propMissing, MakeNAsInBetas, testData[[1]])

plotDat = ErrorAcrossDataSets(x, HousemanBlood150CpGModel)
plotDat = cbind.data.frame(plotDat, propMissing = seq(0,0.9,0.05))

pdf("/mnt/data1/Thea/ErrorMetric/plots/badData/simWithMissingCpGs.pdf", height = 4, width = 7) 
ggplot(plotDat, aes(x = propMissing, y = error)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Proportion of CpGs missing", y = "DSRMSE")
dev.off()



### Check effect of age, sex in EXTEND and US #########
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")
model =  HousemanBlood150CpGModel

## load functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load libraries
library(ggplot2)
library(cowplot)

## load data
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
us = dat
usPheno = pheno
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batches_1_2_normalised_together.rdat")
ex = betas
exPheno = pheno
rm(dat, betas, pheno, HousemanBlood150CpGModel)

## get error for US
usPred = projectCellTypeWithError(us, modelType = "ownModel", ownModelData = model)

## get error for EX
exPred = projectCellTypeWithError(ex, modelType = "ownModel", ownModelData = model)

sexUS = usPheno$nsex
sexUS = ifelse(sexUS =="1", "Male", "Female")
allPheno = rbind.data.frame(cbind.data.frame(usPred, data = "US", age = usPheno$confage, sex = sexUS),
                            cbind.data.frame(exPred, data = "EX", age = exPheno$Age, sex = exPheno$Sex))

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/sexAcrossEXandUS.pdf", height = 7, width = 7) 
ggplot(allPheno, aes(x = data, y = error, fill = sex)) +
  geom_violin() +
  theme_cowplot(18) +
  labs(y = "DSRMSE", x = "Data set", fill = "Sex")
dev.off()

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/ageAcrossEXandUS.pdf", height = 7, width = 7) 
ggplot(allPheno, aes(x = as.numeric(age), y = error, col = data)) +
  geom_point() +
  theme_cowplot(18) +
  labs(y = "DSRMSE", x = "Age", col = "Data")
dev.off()

t.test(allPheno[allPheno$sex == "Female","error"], allPheno[allPheno$sex == "Male","error"], alternative = "greater")

summary(lm(error ~ as.numeric(age), data = allPheno))





### Plot outputs from Essex data
## open gds and make into matrix
library(gdsfmt)
gfile = openfn.gds("/mnt/data1/Thea/ErrorMetric/data/EssexOutput/sub.gds")

dat = cbind.data.frame(read.gdsn(index.gdsn(gfile$root, "Pred")),
                       Age = read.gdsn(index.gdsn(gfile$root, "Age")),
                       Sex = read.gdsn(index.gdsn(gfile$root, "Sex")),
                       Tissue = read.gdsn(index.gdsn(gfile$root, "Tissue")),
                       SubTissue = read.gdsn(index.gdsn(gfile$root, "SubTissue")),
                       DatasetOrigin = read.gdsn(index.gdsn(gfile$root, "DatasetOrigin")))
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")
colnames(dat)[1:8] = c(colnames(HousemanBlood150CpGModel$coefEsts), "error", "nCGmissing")

## remove the .gds from DatasetOrigin
dat$DatasetOrigin = as.factor(unlist(strsplit(as.character(dat$DatasetOrigin), ".g"))[seq(1,nrow(dat)*2,2)])

## remove samples withmissing CpGs
dat = dat[dat$nCGmissing ==0,]

## create a Blood Bool to colour by
dat$blood = rep(0, nrow(dat))
dat$blood[dat$Tissue == "Blood" |
            dat$Tissue == "B Cells" |
            dat$Tissue == "Granulocyes" |
            dat$Tissue == "Neutrophils" |
            dat$Tissue == "NK" |
            dat$Tissue == "Lymph Node" |
            dat$Tissue == "T Cells"] = 1

## plot
library(ggplot2)
library(cowplot)
library(forcats)

pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexsAllTissueBoxplot.pdf", height = 9, width = 13)
ggplot(dat, aes(x = fct_reorder(Tissue, blood, .fun = median, .desc =TRUE), y = error, fill = as.factor(blood))) +
  geom_boxplot() +
  theme_cowplot(18) +
  scale_fill_manual(values = c("#0A8ABA", "#BA3A0A"), name = "Blood?", labels = c("No", "Yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = element_blank(), y = "DSRMSE")
dev.off()

## t test between blood and non blood samples
t.test(dat$error[dat$blood ==0], dat$error[dat$blood ==1])

# ## plot the same for only blood
# datB = dat[dat$blood == 1,]
# ggplot(datB, aes(x = fct_reorder(DatasetOrigin, error, .fun = median, .desc =F), y = error)) +
#   geom_boxplot() +
#   theme_cowplot(18) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Tissue", y = "DSRMSE") 
# 
# 
# ## check relationship with age
# ggplot(datB, aes(x = as.numeric(as.character(Age)), y = error, col = Sex)) +
#   geom_point() +
#   theme_cowplot(18) +
#   labs(x = "Age", y = "DSRMSE") 

## select purified blood cell types
datB = dat[dat$blood == 1,]
datB = datB[!(datB$Tissue == "Blood" | datB$Tissue == "Lymph Node" | datB$Tissue == "Neutrophils"), ]

datB$trueProp = rep(NA, nrow(datB))

datB$trueProp[datB$Tissue == "B Cells"] = datB$Bcell[datB$Tissue == "B Cells"]
datB$trueProp[datB$Tissue == "Granulocyes"] = datB$Gran[datB$Tissue == "Granulocyes"]
datB$trueProp[datB$Tissue == "NK"] = datB$NK[datB$Tissue == "NK"]
datB$trueProp[datB$Tissue == "T Cells"] = datB$CD4T[datB$Tissue == "T Cells"] + datB$CD8T[datB$Tissue == "T Cells"]

ggplot(datB, aes(x = trueProp, y = error, col = Tissue)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "True proportion", y = "DSRMSE")
  
ggplot(datB, aes(x = trueProp, y = error, shape = Tissue, col = DatasetOrigin)) +
  geom_point(size = 2) +
  theme_cowplot(18) +
  labs(x = "True proportion", y = "DSRMSE")

##  plot a grid of one cell type at a time
plot = list()
for (i in 1:4){
  plot[[i]] = ggplot(datB[datB$Tissue == levels(as.factor(as.character(datB$Tissue)))[i],],
                     aes(x = trueProp, y = error, col = DatasetOrigin)) +
    geom_point(size = 2) +
    theme_cowplot(18) +
    labs(x = "True proportion", y = "DSRMSE", coll = "Data") +
    xlim(c(min(datB$trueProp),max(datB$trueProp))) +
    ylim(c(0, max(datB$error))) +
    ggtitle(levels(as.factor(as.character(datB$Tissue)))[i])
}

pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodCellTypevsError.pdf", height = 10, width = 10)
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], labels = "AUTO")
dev.off()