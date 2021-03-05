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

### plot PCA of betas per cell type coloured by train test ######
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")
library(ggfortify)
library(ggplot2)
library(cowplot)
library(scales)

plotDat = list()
colours = hue_pal()(6)
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

plotsnoLeg = plot_grid(plotDat[[1]],
                       plotDat[[2]],
                       plotDat[[3]],
                       plotDat[[4]],
                       plotDat[[5]],
                       plotDat[[6]],
                       ncol = 3, labels = "AUTO")
legplot = ggplot(data.frame(one = betas[1, 1:10], two = betas[2, 1:10], tt = rep(c("Train", "Test"), 5)),
                 aes(x = one, y = two, shape = tt)) + 
  geom_point(size = 2) +
  labs(shape = "") +
  theme_cowplot(18)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/cellTrainTestPCA.pdf", height = 9, width = 12)
plot_grid(plotsnoLeg, get_legend(legplot), 
          ncol = 2,
          rel_widths = c(3, 0.3))
dev.off()

# ## full PCA to show celltype similarity
# betaVar = apply(betas, 1, var, na.rm = T)
# topBeta = betas[order(betaVar, decreasing = T)[1:1000],]
# plotDatPCAFULL = autoplot(prcomp(t(topBeta)), data = pheno, col = "celltype", shape = "celltype", size = 2) +
#   theme_cowplot(18)
# 
# pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/cellPCA.pdf",height = 7, width = 7)
# plotDatPCAFULL +
#   labs(col = "Cell type", shape = "Cell type")
# dev.off()


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
                                   numProbes =  50,
                                   probeSelect = "auto")

## create model in separately normalised train data
sepNormalisedModel = pickCompProbes(rawbetas = quantileBetasTrain,
                                    cellTypes = levels(as.factor(phenoTrain$celltype)),
                                    cellInd = as.factor(phenoTrain$celltype),
                                    numProbes =  50,
                                    probeSelect = "auto")

## create model in combined normalised train data
combNormalisedModel = pickCompProbes(rawbetas = quantileBetas[, pheno$trainTest == "Train"],
                                     cellTypes = levels(as.factor(phenoTrain$celltype)),
                                     cellInd = as.factor(phenoTrain$celltype),
                                     numProbes =  50,
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
                     norm = rep(c("UU", "N", "NN", "UN", "NU"), each = length(absDiffComb)),
                     cell = as.factor(c(as.character(unnormMelt[,2]), as.character(sepMelt[,2]), as.character(combMelt[,2]), as.character(sepMelt[,2]), as.character(sepMelt[,2]))))

## use a paired t-test to compare each group to comb, the default
dat = cbind(absDiffUn, absDiffComb, absDiffSep, absDiffunModNormDat, absDiffnormModUnDat)
library(xtable)
xtable(rbind.data.frame(t.testMatrix(dat), mean = signif(colMeans(dat),3), 
                        SD = signif(apply(dat, 2, sd),3)))

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/sepCombNormalisationComparison.pdf", 
    height = 5, width = 8)
ggplot(plotDat, aes(x = norm, y = diff)) +
  geom_violin() +
  geom_jitter(aes(col = cell, shape = cell)) +
  theme_cowplot(18) +
  labs(x = "Normalisation", y = "Absolute difference between\ntrue and predicted\n(Proportion of methylation)") +
  theme(legend.title = element_blank())
dev.off()



### Optimise number of CpGs ###########################
## load unnormalised data 
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## source model functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

# ## make models
# modelListCpG = list()
# cpg = c(1:25,seq(30,55,5), seq(60, 100, 10))
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
cpg = c(1:25,seq(30,55,5), seq(60, 100, 10))
trueProp = singleCellProportionMatrix(phenoTest[,1])
rmseTvP = c()

# RMSE = function(m, o){
#   sqrt(mean((m - o)^2))
# }


for (i in 1:length(cpg)) {
  pred = projectCellTypeWithError(betasTest,  modelType = "ownModel", ownModelData = modelListCpG[[i]])
  rmseTvP = c(rmseTvP, mean(abs(melt(trueProp)[,3] - melt(pred[,-which(colnames(pred) %in% c("nCGmissing", "error"))])[,3])))
}

nCG = c()
for(i in 1:length(modelListCpG)){
  nCG = c(nCG, nrow(modelListCpG[[i]]$coefEsts))
}


corPlot = data.frame(rmseTvP, cpg, nCG)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/nCpGNeededForModel.pdf", height = 7, width = 7)
ggplot(corPlot, aes(x = cpg, y = rmseTvP)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = expression(paste(italic("numProbes"))), 
                      y = "Absolute difference between\ntrue and predicted\n(Proportion of methylation)")
dev.off()

## save model that used 150 CpGs
HousemanBlood50CpGModel = modelListCpG[[which(cpg == 50)]]
save(HousemanBlood50CpGModel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")



### Plot heirarchical cluster of training cpgs used in model #####
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

## load data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## load functions
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

modelBetas = GetModelCG(betasTrain, list(HousemanBlood50CpGModel))
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

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/heatmapForModelCpGs.pdf", height = 6, width = 5)
Heatmap(modelBetas, name = "DNAm",
        top_annotation = ha, show_row_names = F, show_column_names = F, show_row_dend = F)
dev.off()

### plot PCA of model CpGs ############################
## full PCA to show celltype similarity
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

library(ggfortify)
library(ggplot2)
library(cowplot)
library(scales)

pheno = data.frame(celltype = colnames(HousemanBlood50CpGModel$coefEsts))
plotDatPCAFULL = autoplot(prcomp(t(HousemanBlood50CpGModel$coefEsts)), data = pheno, col = "celltype", shape = "celltype", size = 2) +
  theme_cowplot(18)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/modelCpGPCA.pdf",height = 6, width = 6)
plotDatPCAFULL +
  labs(col = "Cell type", shape = "Cell type")
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
#                                   numProbes = 50,
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
plotDat5$mod = factor(unlist(strsplit(as.character(plotDat5$model), "l_"))[seq(2,nrow(plotDat5)*2,2)],
                         levels = c("C4C8GMNK", "BC8GMNK", "BC4GMNK", "BC4C8MNK", "BC4C8GNK", "BC4C8GM"))


## violin plot of 5 cell types only
pdf("/mnt/data1/Thea/ErrorMetric/plots/badModels/violin5CelltypeModels.pdf", height = 6, width = 6)
ggplot(plotDat5, aes(x = mod, y = error, fill = mod)) +
  geom_violin() +
  theme_cowplot(18) +
  labs(x = "Model", y = "DSRMSE") +
  ylim(c(0, max(plotDat5$error))) +
  scale_x_discrete(limits = c("C4C8GMNK", "BC8GMNK", "BC4GMNK", "BC4C8MNK", "BC4C8GNK", "BC4C8GM")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()


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
  # theme(axis.text.x=element_blank(),
  # axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Sa" = "0.1",
                              "Sb" = "0.2",
                              "Sc" = "0.3",
                              "Sd" = "0.4",
                              "Se" = "0.5",
                              "Sf" = "0.6",
                              "Sg" = "0.7",
                              "Sh" = "0.8",
                              "Si" = "0.9",
                              "Sj" = "1.0")) +
  labs(x = "Proportion of missing cell type", y = "DSRMSE", col = "Model") +
  ylim(c(0, max(erPlotDat$error)))
dev.off()





### Compare error when noise is introduced ############

## load functions
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

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
                                          modelList = list(Predicted = HousemanBlood50CpGModel), 
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
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

## load testing data
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")

## mean proportions of each cell type in whole blood (from Reinius2012)
meanBloodProp = matrix(nrow = 1, byrow = T, data = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43))/100
colnames(meanBloodProp) = levels(phenoTest$celltype)

## create a single representative sample
testData = CellTypeProportionSimulator(betas = GetModelCG(betasTest, list(HousemanBlood50CpGModel)), 
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

plotDat = ErrorAcrossDataSets(x, HousemanBlood50CpGModel)
plotDat = cbind.data.frame(plotDat, propMissing = seq(0,0.9,0.05))

pdf("/mnt/data1/Thea/ErrorMetric/plots/badData/simWithMissingCpGs.pdf", height = 4, width = 7) 
ggplot(plotDat, aes(x = propMissing, y = error)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Proportion of CpGs missing", y = "DSRMSE")
dev.off()



### Check effect of age, sex in EXTEND and US #########
## load model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")
model =  HousemanBlood50CpGModel

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
rm(dat, betas, pheno, HousemanBlood50CpGModel)

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
  labs(y = "DSRMSE", x = "Dataset", fill = "Sex")
dev.off()

## create age from 38 phenotype
allPheno$ageDiff = abs(as.numeric(as.character(allPheno$age))-38)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/ageAcrossEXandUS.pdf", height = 7, width = 7) 
ggplot(allPheno, aes(x = ageDiff, y = error, col = data)) +
  geom_point() +
  theme_cowplot(18) +
  labs(y = "DSRMSE", x = "Absolute age difference", col = "Dataset")
dev.off()

t.test(allPheno[allPheno$sex == "Female","error"], allPheno[allPheno$sex == "Male","error"], alternative = "greater")

summary(lm(error ~ ageDiff, data = allPheno))


### check annotation of model CpGs to chromosomes #####
x = read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7, header = T)
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

x = x[,c("IlmnID", "CHR")]
cpgInMod = x[x$IlmnID %in% rownames(HousemanBlood50CpGModel$coefEsts),]
table(cpgInMod$CHR)

save(cpgInMod, file = "/mnt/data1/Thea/ErrorMetric/data/cpgInModel.Rdata")

### compare X chromosome in model between male and female in EX and US ####
load("/mnt/data1/Thea/ErrorMetric/data/cpgInModel.Rdata")
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
us = dat
usPheno = pheno
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batches_1_2_normalised_together.rdat")
ex = betas
exPheno = pheno
rm(dat, betas, pheno)

exT = ex[rownames(ex) %in% cpgInMod$IlmnID,]
usT = us[rownames(us) %in% cpgInMod$IlmnID,]

exT = exT[match(cpgInMod$IlmnID, rownames(exT)),]
usT = usT[match(cpgInMod$IlmnID, rownames(usT)),]

all(rownames(exT) == cpgInMod$IlmnID)
all(rownames(usT) == cpgInMod$IlmnID)

## subset for those in the x chromosome
exT = exT[cpgInMod$CHR == "X",]
usT = usT[cpgInMod$CHR == "X",]

sexUS = usPheno$nsex
sexUS = ifelse(sexUS =="1", "Male", "Female")

plotDat = rbind.data.frame(data.frame(t(exT), sampleID = colnames(exT), study = "EX", sex = exPheno$Sex),
                           data.frame(t(usT), sampleID = colnames(usT), study = "US", sex = sexUS))

statDat = data.frame(matrix(nrow = ncol(plotDat)-3, ncol = 5))
colnames(statDat) = c("CpG", "pValue", "Sig", "meanF", "meanM")
for (i in 1:(ncol(plotDat)-3)){
  x = t.test(plotDat[,i]~plotDat[,"sex"])
  statDat[i,1] = colnames(plotDat)[i]
  statDat[i,2] = signif(x$p.value, 2)
  statDat[i,4:5] = signif(x$estimate, 3)
}
statDat$Sig[statDat$pValue >= 0.05] = ""
statDat$Sig[statDat$pValue < 0.05] = "."
statDat$Sig[statDat$pValue < 0.01] = "*"
statDat$Sig[statDat$pValue < 0.001] = "**"
statDat$Sig[statDat$pValue < 0.0001] = "***"


library(xtable)
print(xtable(statDat), include.rownames=FALSE)

library(reshape2)
plotDat = melt(plotDat, id.vars = c("study", "sex", "sampleID"))

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/ErrorXchrCpGsInEXUS.pdf", height = 8, width = 14)
ggplot(plotDat, aes(x = variable, y = value, fill = sex)) +
  theme_cowplot(18) +
  geom_violin(position=position_dodge(0.5)) +
  facet_wrap(~study, ncol = 1) +
  labs(y = "Proportion of methylation", fill = "Sex", x = "X chromosome CpGs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15))
dev.off()




### Plot outputs from Essex data ######################

# ## Essex server code:
# ## open gds file
# library(gdsfmt)
# x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)
# 
# 
# betas = read.gdsn(index.gdsn(x, "rawbetas"))
# #dim(betas)
# 
# load("~/DSRMSE/models/HousemanBloodModel50CpG.Rdata")
# 
# pID = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
# betaIndex = pID %in% rownames(HousemanBlood50CpGModel$coefEsts)
# 
# rownames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
# colnames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "FullBarCode"))
# 
# betaMod = betas[betaIndex,]
# betaMod = betaMod[match(rownames(HousemanBlood50CpGModel$coefEsts), rownames(betaMod)),]
# 
# ## create column index, removing samples that have "", or unsorted
# tDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Tissue"))
# aDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Age"))
# sDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Sex"))
# dDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "DatasetOrigin"))
# stDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "SubTissue"))
# tInd = which(tDat == "" | 
#                tDat == "Unsorted Tissues" | 
#                tDat == "Unsorted Cell Line" | 
#                tDat == "Unsorted Tumours")
# 
# betaMod = betaMod[,-tInd]
# 
# gfile = createfn.gds("sub.gds")
# 
# add.gdsn(gfile, "Age", aDat[-tInd])
# add.gdsn(gfile, "Sex", sDat[-tInd])
# add.gdsn(gfile, "Tissue", tDat[-tInd])
# add.gdsn(gfile, "DatasetOrigin", dDat[-tInd])
# add.gdsn(gfile, "SubTissue", stDat[-tInd])
# add.gdsn(gfile, "rownames", rownames(betaMod))
# add.gdsn(gfile, "colnames", colnames(betaMod))
# 
# #source("~/DSRMSE/FunctionsForErrorTesting.R")
# source("~/DSRMSE/projectCellTypeWithError.R")
# model = HousemanBlood50CpGModel
# 
# ind = c(seq(1000, 20960, 1000))
# errPred = lapply(ind, function(ind){
#   pred = projectCellTypeWithError(betaMod[,(ind - 999):ind], model = "ownModel", ownModelData = model)
#   return(pred)})
# 
# errPred[[length(errPred)+1]] = projectCellTypeWithError(betaMod[,20001:20960], model = "ownModel", ownModelData = model)
# 
# ## extract all from list and add to gds file
# pred = errPred[[1]]
# for (i in 2:length(errPred)){
#   pred = rbind(pred, errPred[[i]])
# }
# 
# add.gdsn(gfile, "Pred", pred)
# closefn.gds(gfile)
# closefn.gds(x)
# q()
# 
# scp sub.gds dSeiler@knight.ex.ac.uk:/mnt/data1/Thea/ErrorMetric/data/EssexOutput/
#   

## open gds and make into matrix
library(gdsfmt)
gfile = openfn.gds("/mnt/data1/Thea/ErrorMetric/data/EssexOutput/sub.gds")

dat = cbind.data.frame(read.gdsn(index.gdsn(gfile$root, "Pred")),
                       Age = read.gdsn(index.gdsn(gfile$root, "Age")),
                       Sex = read.gdsn(index.gdsn(gfile$root, "Sex")),
                       Tissue = read.gdsn(index.gdsn(gfile$root, "Tissue")),
                       SubTissue = read.gdsn(index.gdsn(gfile$root, "SubTissue")),
                       Sample = read.gdsn(index.gdsn(gfile$root, "colnames")),
                       DatasetOrigin = read.gdsn(index.gdsn(gfile$root, "DatasetOrigin")))
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")
colnames(dat)[1:8] = c(colnames(HousemanBlood50CpGModel$coefEsts), "error", "nCGmissing")

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

## merge bloods for plot
dat$TissueBlood = dat$Tissue
dat$TissueBlood[dat$Tissue == "Blood" |
            dat$Tissue == "B Cells" |
            dat$Tissue == "Granulocyes" |
            dat$Tissue == "Neutrophils" |
            dat$Tissue == "NK" |
            dat$Tissue == "Lymph Node" |
            dat$Tissue == "T Cells"] = "Blood"

## close gds 
closefn.gds(gfile)

## plot
library(ggplot2)
library(cowplot)
library(forcats)
library(dplyr)

dat_summary = dat %>%
  group_by(TissueBlood) %>%
  tally()

dat = merge(dat, dat_summary,  by = "TissueBlood")

dat$TissueBlood = as.factor(as.character(dat$TissueBlood))
pos = c()
for(i in 1:length(levels(dat$TissueBlood))){
  pos = c(pos, max(dat$error[dat$Tissue == levels(dat$TissueBlood)[i]]))
}

dat.pos = data.frame(TissueBlood = levels(dat$TissueBlood), pos, n = dat_summary$n)

pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexsAllTissueBoxplot.pdf", height = 9, width = 14)
ggplot(dat, aes(x = fct_reorder(TissueBlood, blood, .fun = median, .desc =TRUE))) +
  geom_boxplot(aes(y = error, fill = as.factor(blood))) +
  theme_cowplot(18) +
  scale_fill_manual(values = c("#0A8ABA", "#BA3A0A"), name = "Blood?", labels = c("No", "Yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = element_blank(), y = "DSRMSE") +
  geom_text(data = dat.pos, aes(TissueBlood, label = n, y = pos+0.02))
dev.off()

## t test between blood and non blood samples
t.test(dat$error[dat$blood ==0], dat$error[dat$blood ==1])

# ## plot the same for only blood
# datB = dat[dat$blood == 1,]
# pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodBoxplot.pdf", height = 10, width = 16)
# ggplot(datB, aes(x = fct_reorder(DatasetOrigin, error, .fun = median, .desc =F), y = error, fill = "#BA3A0A")) +
#   geom_boxplot() +
#   theme_cowplot(18) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.x.bottom = element_text(size=12)) +
#   labs(x = "Study", y = "DSRMSE")
# dev.off()
# 
# 
# ## check relationship with age
# ## remove those without age or sex
# datBAS = datB[!datB$Sex == "",]
# pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodSexCheck.pdf", height = 10, width = 16)
# ggplot(datBAS, aes(x = fct_reorder(DatasetOrigin, error, .fun = median, .desc =F), y = error, fill = as.factor(as.character(Sex)))) +
#   geom_boxplot() +
#   theme_cowplot(18) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.x.bottom = element_text(size=12)) +
#   labs(x = "Study", y = "DSRMSE", fill = "Sex")
# dev.off()
# 
# pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodAgeCheck.pdf", height = 13, width = 9)
# ggplot(datBAS, aes(x = as.numeric(as.character(Age)), y = error, col = as.factor(as.character(DatasetOrigin)))) +
#   geom_point(size = 1.4) +
#   theme_cowplot(18) +
#   theme(legend.position = "bottom", legend.text = element_text(size=12)) +
#   labs(x = "Age", y = "DSRMSE", col = "Study") 
# dev.off()

## select purified blood cell types
datB = dat[dat$blood == 1,]
datB = datB[!(datB$Tissue == "Blood" | datB$Tissue == "Lymph Node" | datB$Tissue == "Neutrophils"), ]

datB$trueProp = rep(NA, nrow(datB))

datB$trueProp[datB$Tissue == "B Cells"] = datB$Bcell[datB$Tissue == "B Cells"]
datB$trueProp[datB$Tissue == "Granulocyes"] = datB$Gran[datB$Tissue == "Granulocyes"]
datB$trueProp[datB$Tissue == "NK"] = datB$NK[datB$Tissue == "NK"]
datB$trueProp[datB$Tissue == "T Cells"] = datB$CD4T[datB$Tissue == "T Cells"] + datB$CD8T[datB$Tissue == "T Cells"]

datB$Tissue = as.character(datB$Tissue)
datB$Tissue[datB$Tissue == "Granulocyes"] = "Granulocytes"
datB$Tissue = as.factor(as.character(datB$Tissue))

##  plot a grid of one cell type at a time
plot = list()
for (i in 1:4){
  plot[[i]] = ggplot(datB[datB$Tissue == levels(as.factor(as.character(datB$Tissue)))[i],],
                     aes(x = trueProp, y = error, col = DatasetOrigin)) +
    geom_point(size = 2) +
    theme_cowplot(18) +
    labs(x = "Predicted proportion", y = "DSRMSE", coll = "Data") +
    xlim(c(min(datB$trueProp),max(datB$trueProp))) +
    ylim(c(0, max(datB$error))) +
    ggtitle(levels(as.factor(as.character(datB$Tissue)))[i])
}

pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodCellTypevsError.pdf", height = 10, width = 15)
plot_grid(plot[[1]], plot[[2]], plot[[3]], plot[[4]], labels = "AUTO")
dev.off()


## stacked bar charts for each study in each cell type
## B cells
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

# predictions = datB[datB$Tissue == "T Cells",]
# 
# temp = cellTypeCompareStackedBar(predictions)

datB = datB[,-which(colnames(datB) =="Sample")]
datB = datB[,-which(colnames(datB) =="n")]
datB = datB[,-which(colnames(datB) =="TissueBlood")]

plotList = list()
datB$Tissue = as.factor(as.character(datB$Tissue))
for(i in 1:length(levels(datB$Tissue))){
  plotList = c(plotList, 
               cellTypeCompareStackedBar(
                 datB[datB$Tissue == levels(datB$Tissue)[i],]))
}

# pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodStackedBar.pdf")
# for(i in 1:length(plotList)){
#   print(plotList[[i]])
# }
# dev.off()

pdf("/mnt/data1/Thea/ErrorMetric/plots/EssexDataPlots/ErrorEssexBloodStackedBarExample.pdf", height = 9, width = 7)
print(plotList[[3]])
dev.off()



## plot trueProp against 1-sum(predicted)
datB$totalPred = rowSums(datB[,c("Bcell", "CD4T", "CD8T","Gran", "Mono","NK")])
datB$absDiffPred = abs(rowSums(datB[,c("Bcell", "CD4T", "CD8T","Gran", "Mono","NK")])-1)


ggplot(datB, aes(x = trueProp, y = absDiffPred, col = DatasetOrigin, shape = Tissue)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "True proportion", y = "|Sum of predictions - 1|")

ggplot(datB, aes(x = trueProp, y = absDiffPred, col = error, shape = Tissue)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "True proportion", y = "|Sum of predictions - 1|")


## doesn't show anything! Now look up the worst datasets in GEO in the hope that they're all cancer or something...
datBad = datB[datB$error <0.1 & datB$trueProp<0.75,]
levels(as.factor(as.character(datBad$DatasetOrigin)))

## get sample IDs for samples in GSE89251 with low error, bad pred and those with higher error
datBad[datBad$DatasetOrigin == "GSE89251" & datBad$trueProp < 0.5,"Sample"]
datB[datB$DatasetOrigin == "GSE89251" & datB$error > 0.2,"Sample"]


### comparison to estimateCellCounts.wmln #############
# library(gdsfmt)
# x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)
# 
# ### comparison to estimateCellCounts.wmln #############
# library(gdsfmt)
# x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)
# 
# library(wateRmelon)
# library(bigmelon)
# ## function for creating subgds
# ##  INPUT newGdsName - name of new file to be saved
# ##        gfileOld - file to be subset
# ##        pDataSplit - the pData column the data will be split using
# ##        pDataVal - the pData value wanted
# 
# setwd("~/gdsTemp/")
# 
# pDataVal = c("B Cells", "Granulocyes", "NK", "T Cells")
# newGdsName = "gsub"
# pDataSplit = "Tissue"
# gfileOld = x
# 
# 
# createSubGDS = function(pDataVal, newGdsName, pDataSplit, gfileOld){
#   gfile = createfn.gds(newGdsName)
#   
#   ## add fData
#   add.gdsn(gfile, "fData", fData(x))
#   
#   ## get index for col/rows to keep
#   pCol = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), pDataSplit))
#   index = which(pCol %in% pDataVal)
#   
#   ## add subset pData
#   add.gdsn(gfile, "pData", pData(x)[index,])
#   
#   ## add the rest of the data, subset
#   add.gdsn(gfile, "rawbetas", read.gdsn(index.gdsn(x, "rawbetas"))[,index])
#   add.gdsn(gfile, "methylated", read.gdsn(index.gdsn(x, "methylated"))[,index])
#   add.gdsn(gfile, "unmethylated", read.gdsn(index.gdsn(x, "unmethylated"))[,index])
#   add.gdsn(gfile, "pvals", read.gdsn(index.gdsn(x, "pvals"))[,index])
#   add.gdsn(gfile, "NBeads", read.gdsn(index.gdsn(x, "NBeads"))[,index])
#   
#   return(gfile)
# }





### Predict proportions in E risk purified samples (just those that overlap) ####
library(minfi)
load("/mnt/data1/Eilis/Projects/Asthma/QC/CombinedQC_WithMarkdown/AllRGSet.rdat")
betas = getBeta(preprocessRaw(RGSet))

pheno = read.csv("/mnt/data1/Eilis/Projects/Asthma/QC/CombinedQC_WithMarkdown/CellSortedAsthmaERisk_SamplesPassedQC.csv")
pheno = pheno[,c("Basename", "PatientID", "Sample.Type", "sex", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")]

## only keep cell types in prediction 
pheno = pheno[which(pheno$Sample.Type == "B-cells" |
                      pheno$Sample.Type == "CD4 T-cells" |
                      pheno$Sample.Type == "CD8 T-cells" |
                      pheno$Sample.Type == "Granulocytes" |
                      pheno$Sample.Type == "Monocytes"),]

## subset betas for those in pheno
betas = betas[,colnames(betas) %in% pheno$Basename]
# all(colnames(betas) == pheno$Basename)  # T

## Use model to predict cell types
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

pred = projectCellTypeWithError(betas, modelType = "ownModel", ownModelData = HousemanBlood50CpGModel)

## summarise data for plotting
predNames =  c("CD4T", "CD8T", "Bcell", "Gran", "Mono")
eRiskNames = c("CD4 T-cells","CD8 T-cells", "B-cells", "Granulocytes", "Monocytes")
cellP = cellM = Celltype = error = c()

for (i in 1:length(predNames)){
  temp = pred[pheno$Sample.Type == eRiskNames[i],]
  cellP = c(cellP, temp[,colnames(temp) == predNames[i]])
  Celltype = c(Celltype, rep(predNames[i], nrow(temp)))
  error = c(error, temp[,"error"])
  pTemp = pheno[pheno$Sample.Type == eRiskNames[i],]
  cellM = c(cellM, pTemp[,colnames(pTemp) == predNames[i]])
}

plotDat = data.frame(cellP, cellM, Celltype, error)

library(ggplot2)  
library(cowplot)  
library(scales)

colSub = hue_pal()(6)

## get correlation for each cell type
core = c()
for (i in 1:5){ 
  temp = plotDat[plotDat$Celltype == levels(plotDat$Celltype)[i],]
  core = c(core, cor(temp$cellP, temp$error))
}
core = signif(core,2)

labin = paste(levels(plotDat$Celltype), " (Cor = ", core, ")", sep = "")
names(labin) = levels(plotDat$Celltype)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/EriskCelltypeError.pdf", height = 11, width = 8)
ggplot(plotDat, aes(x = cellP, y = error, col = Celltype, shape = Celltype)) +
  geom_point(size = 2) +
  scale_color_manual(values = colSub[1:5]) +
  scale_shape_manual(values = c(20, 17, 15, 3, 7)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_cowplot(18) +
  facet_wrap(~Celltype, nrow = 5, labeller = labeller(Celltype = labin)) +
  labs(x = "Predicted proportion", y = "DSRMSE", col = "Cell type", shape = "Cell type") +
  theme(legend.position = "none")
dev.off()


