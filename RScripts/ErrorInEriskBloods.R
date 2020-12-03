## Script to test error metric using E-risk data set
## Dorothea Seiler Vellame
## started 09-10-2020

## knight file path 
path = "/mnt/data1/Thea/ErrorMetric/"

## ISCA file path
path = "/gpfs/ts0/projects/Research_Project-191406/ErrorMetric/"

## get stats on samples #####
load("/mnt/data1/Eilis/Projects/Asthma/CellTypeComparisons/Correlations_New/BetasSortedByCellType_NoY.rdat")

tableOfSamples = data.frame(matrix(ncol = 11, nrow = 30, data = 0))
colnames(tableOfSamples) = c("Sample", types, "Train", "Test")
tableOfSamples$Sample = individuals

for(i in 1:8){
  tableOfSamples[,i+1] = 1
  tableOfSamples[which(is.na(allphenos[[i]][,"PatientID"])),i+1] = 0
}
 
load(paste(path, "data/bloodEriskDataFull.Rdata", sep = ""))

trainTestInfo = unique(bloodPheno[,c("PatientID", "traintest")])
tableOfSamples$Train[which(tableOfSamples$Sample %in% trainTestInfo[which(trainTestInfo$traintest == "train"),1])] = 1
tableOfSamples$Test[which(tableOfSamples$Sample %in% trainTestInfo[which(trainTestInfo$traintest == "test"),1])] = 1
  
tableOfSamples = rbind.data.frame(tableOfSamples, c("Total", colSums(tableOfSamples[,2:ncol(tableOfSamples)])))
write.csv(tableOfSamples, file = paste(path, "data/tableOfERiskSamples.csv", sep = ""))
  
## QC original data ###################################
## make PCAs, remove non blood data and NA samples
## assign train test status and show randomness
## save train test separately


## load the data
load("/mnt/data1/Eilis/Projects/Asthma/CellTypeComparisons/Correlations_New/BetasSortedByCellType_NoY.rdat")

## create bloodbeta list containing only the blood realted data INCLUDES WHOLE BLOOD
bloodbeta = allbetas[2:7]
bloodpheno = allphenos[2:7]

## are any samples missing consistently? If so remove them
for(i in 2:7){
  print(which(is.na(allphenos[[i]][,"PatientID"])))
}
## remove samples 2 and 27 from pheno and beta

for (i in 1:6){
  bloodbeta[[i]] = bloodbeta[[i]][,-c(2, 27)]
  bloodpheno[[i]] = bloodpheno[[i]][-c(2, 27),]
}

## create blood pheno with sample ID , cell type, predicted age, sex
bloodphenomerge = bloodpheno[[1]][,c("PatientID", "Sample.Type", "sex", "PredictedAge")]
for (i in 2:6){
  bloodphenomerge = rbind.data.frame(bloodphenomerge, bloodpheno[[i]][,c("PatientID", "Sample.Type", "sex", "PredictedAge")])
}

## check predicted ages per group
library(ggplot2)
library(cowplot)

# pdf(paste(path, "plots/createTrainTest/predictedAgeAcrossAllCellTypes.pdf", sep = ""), height = 5, width = 8)
# ggplot(bloodphenomerge, aes(x = PatientID, y = PredictedAge, col = Sample.Type)) +
#   geom_point() +
#   theme_cowplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Patient ID", y = "Predicted age", col = "Cell type")
# dev.off()

## non whole blood will vary as the age predictor probably not built for cell populations
## Don't use age, ages are very similar

wholebloodphenomerge = bloodphenomerge[which(bloodphenomerge$Sample.Type == "whole blood"),]

## split by sex
set.seed(1234)
testMIndex = sample(which(wholebloodphenomerge$sex =="M"), 5)
trainMIndex = which(wholebloodphenomerge$sex =="M")[which(!which(wholebloodphenomerge$sex =="M") %in% testMIndex)]
set.seed(1234)
testFIndex = sample(which(wholebloodphenomerge$sex =="F"), 5)
trainFIndex = which(wholebloodphenomerge$sex =="F")[which(!which(wholebloodphenomerge$sex =="F") %in% testFIndex)]

## Add train test column to pheno file
wholebloodphenomerge[c(trainMIndex, trainFIndex),5] = "train" 
wholebloodphenomerge[c(testMIndex, testFIndex),5] = "test"

# pdf(paste(path, "plots/createTrainTest/trainTestSamplesAcrossPredictedAgeAndSex.pdf", sep = ""), height = 5, width = 5)
# ggplot(wholebloodphenomerge, aes(x = sex, y = PredictedAge, col = V5)) +
#   geom_jitter(width = 0.3, size = 2) +
#   theme_cowplot() +
#   labs(x = "Sex", y = "Predicted age", col = "")
# dev.off()


## make train test PCA across each cell type to show randomness
getMostVarN = function(data, n){
  sd = apply(data, 1, sd, na.rm = T)
  top = sort(sd, decreasing = T)[1:n]  
  return(data[rownames(data) %in% names(top),])
}

## make PCA of 100 most variable sites across all cells to check no mislabeling
topbetas = list()
for (i in 1:length(bloodbeta)){
  topbetas[[i]] = getMostVarN(bloodbeta[[i]], 100)
}

plotDat = list()
for (i in 1:length(bloodbeta)){
  pcaDat = prcomp(t(topbetas[[i]]))
  plotDat[[i]] = data.frame(pcaDat$x[,1:2], wholebloodphenomerge, bloodpheno[[i]]$Sample.Type)
  colnames(plotDat[[i]])[8] = "celltype"
}

plotL = list()
for  (i in 1:length(bloodbeta)){
  plotL[[i]] = ggplot(plotDat[[i]], aes(x = PC1, y = PC2, col = V5, shape = sex)) +
    geom_point(size = 2) +
    theme_cowplot(18) +
    ggtitle(plotDat[[i]]$celltype[1]) +
    labs(col = "")
}

pdf(paste(path, "plots/createTrainTest/PCAPerCellTypeTrainTest.pdf", sep = ""), height = 6, width = 15)
plot_grid(plotL[[1]], 
          plotL[[2]], 
          plotL[[3]], 
          plotL[[4]], 
          plotL[[5]], 
          plotL[[6]], 
          nrow = 2)
dev.off()

## add train test to each pheno
for  (i in 1:length(bloodbeta)){
  bloodpheno[[i]] = cbind.data.frame(bloodpheno[[i]], wholebloodphenomerge$V5)
  colnames(bloodpheno[[i]])[30] = "traintest"
}

bloodBeta = cbind.data.frame(bloodbeta[[1]],
                             bloodbeta[[2]],
                             bloodbeta[[3]],
                             bloodbeta[[4]],
                             bloodbeta[[5]],
                             bloodbeta[[6]])

bloodPheno = rbind.data.frame(bloodpheno[[1]],
                              bloodpheno[[2]],
                              bloodpheno[[3]],
                              bloodpheno[[4]],
                              bloodpheno[[5]],
                              bloodpheno[[6]])

## remove spaces from cell type phenotype
bloodPheno[which(bloodPheno$Sample.Type == "CD4 T-cells"), "Sample.Type"] = "CD4.T_cells"
bloodPheno[which(bloodPheno$Sample.Type == "CD8 T-cells"), "Sample.Type"] = "CD8.T_cells"
bloodPheno[which(bloodPheno$Sample.Type == "B-cells"), "Sample.Type"] = "B_cells"

## save blood data
save(bloodBeta, bloodPheno,
     file = paste(path, "data/bloodEriskDataFull.Rdata", sep = ""))

## save train and test data
betas = bloodBeta[,intersect(which(bloodPheno$traintest == "train"), which(bloodPheno$Sample.Type != "whole blood"))]
betastest = bloodBeta[,intersect(which(bloodPheno$traintest == "test"), which(bloodPheno$Sample.Type != "whole blood"))]

pheno = bloodPheno[intersect(which(bloodPheno$traintest == "train"), which(bloodPheno$Sample.Type != "whole blood")),]
phenotest = bloodPheno[intersect(which(bloodPheno$traintest == "test"), which(bloodPheno$Sample.Type != "whole blood")),]

save(betas, pheno,
     file = paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))
save(betastest, phenotest,
     file = paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))









## Verify that 150 CpG per cell type is enough ########
## Don't include plots with error as that's not in the narative yet! This should just show that 150 is more than enough
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))
load(paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))

## Plot number of CpGs vs correlation between actual & pred 
## make models 
modelListCpG = list()
cpg = c(1:25,seq(30,55,5), seq(60, 150, 10))
for (i in 1:length(cpg)) {
  modelListCpG[[i]] = pickCompProbes(rawbetas = as.matrix(betas),
                                     cellTypes = levels(as.factor(pheno$Sample.Type)),
                                     cellInd = as.factor(pheno$Sample.Type),
                                     numProbes =  cpg[i],
                                     probeSelect = "auto")
  names(modelListCpG)[i] = paste("model", cpg[i], "CG", sep = "")
}

## save model list
save(modelListCpG, file = paste(path, "data/CpGModels.Rdata", sep = ""))

load (paste(path, "data/CpGModels.Rdata", sep = ""))

## simulate 40 bulk samples
bulk = CellTypeProportionSimulator(betastest, phenotest, "Sample.Type", 40)
bulkBetas = bulk[[1]]
bulkPheno = bulk[[2]]


## get predicted value for each model
truevpred = list()
corDat = c()
for (i in 1:length(cpg)) {
  pred = projectCellTypeWithError(GetModelCG(bulkBetas, list(modelCG = modelListCpG[[i]])) , modelListCpG[[i]]$coefEsts)
  truevpred[[i]] = TrueVSPredictedPlot(pred, bulkPheno)
  corDat = c(corDat, cor(truevpred[[i]]$data$proportion_pred, truevpred[[i]]$data$proportion_true))
}

corPlot = data.frame(corDat, cpg)

p1 = ggplot(corPlot, aes(x = cpg, y = corDat)) +
  geom_point() +
  theme_cowplot(18) +
  labs(x = "Number of CpGs", y = "Correlation") +
  ylim(c(0.9965, 0.999))


## Use model with 150 CpGs
p2 = truevpred[[which(cpg == 150)]]

## cor of p2
corDat[which(cpg == 150)]

pdf(paste(path, "plots/accuracyOf150CGModel/model150CGComboPlots.pdf", sep = ""), height = 4, width = 11)
plot_grid(p1,
          p2,
          rel_widths = c(1,1.35),
          ncol = 2,
          labels = "AUTO")
dev.off()



## Create models using 3:n cell types #################

## For each cell type there will be a model name shorthand:
## B = B cell 
## C4 = CD4 T cell
## C8 = CD8 T cell
## G = Granulocytes
## M = Monocytes

source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))

## make design matrix of booleans for which cell type will be present
designMatrix = expand.grid(c(T,F), c(T,F), c(T,F), c(T,F), c(T,F))
designMatrix = designMatrix[apply(designMatrix, 1, sum) >= 2,]

cellTypes = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes", "Monocytes")
cellTypeShorthand = c("B", "C4", "C8", "G", "M")

# modelList = list()
# for (i in 1:nrow(designMatrix)){
#   modellingDat = CellTypeSubsetBetasAndPheno(cellTypes[unlist(designMatrix[i,])],
#                                              betas, pheno, phenoColName = "Sample.Type", justBetas = F)
#   modelList[[i]] = pickCompProbes(rawbetas = as.matrix(modellingDat[[1]]),
#                                   cellTypes = levels(as.factor(modellingDat[[2]]$Sample.Type)),
#                                   cellInd = as.factor(modellingDat[[2]]$Sample.Type),
#                                   numProbes =  150,
#                                   probeSelect = "auto")
#   names(modelList)[i] = paste("model", paste(cellTypeShorthand[unlist(designMatrix[i,])], sep = "", collapse = ""), sep = "_")
# }
# 
# save(modelList, file = paste(path, "data/VaryNCellsData.Rdata", sep = ""))
# 
# bulk = CellTypeProportionSimulator(GetModelCG(betas, modelList),
#                             pheno,
#                             phenoColName = "Sample.Type",
#                             nBulk = 15,
#                             proportionsMatrixType = "random")
# 
# stackedPlots = ModelCompareStackedBar(bulk[[1]],
#                            modelList, 
#                            nCpGPlot = F,
#                            sampleNamesOnPlots = F,
#                            trueComparison = T,
#                            trueProportions = bulk[[2]])
# 
# save(modelList, bulk, stackedPlots, file = paste(path, "data/VaryNCellsData.Rdata", sep = ""))

## load the models
load(paste(path, "data/VaryNCellsData.Rdata", sep = ""))
x = stackedPlots[[1]]$data

## add columns for which cell types in each data set and then compare specific models
modelPresent = matrix(ncol = length(cellTypes), nrow = nrow(x), data = 0)
colnames(modelPresent) = cellTypeShorthand

for(i in 1:5){
  modelPresent[grep(cellTypeShorthand[i],x$model),i] = 1
}

y = rowSums(modelPresent)
plotDat = data.frame(x,modelPresent, sums = y)

## colour by number of cells in model
plotDatBox = spread(plotDat, key = c(cellType), value = c(proportion_pred))

ggplot(plotDatBox, aes(x = as.factor(sums), y = error, fill = as.factor(sums))) +
  geom_boxplot() +
  geom_jitter() +
  theme_cowplot(18) +
  labs(x = "Number of cell types in the model", y = "DSRMSE") +
  ylim(c(0, max(plotDatBox$error))) +
  theme(legend.position = "none") 

save(plotDatBox, file = paste(path, "data/boxplotDat.Rdata", sep = ""))
### TO DO IN ISCA #####################################
load(paste(path, "data/boxplotDat.Rdata", sep = ""))
library(ggpubr)

my_comparisons <- list( c("5", "4"), c("5", "3"), c("5", "2") )
  
pdf(file = paste(path, "plots/cellTypesInModel/nCellTypeBoxplot.pdf", sep = ""), height = 6, width = 7)
ggboxplot(plotDatBox, y = "error", x = "sums", fill = "sums") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_cowplot(18) +
  geom_jitter() +
  labs(x = "Number of cell types in the model", y = "DSRMSE") +
  theme(legend.position = "none")
dev.off()

#######################################################


## general Q: 
## why do some still predict well? low proportion of that cell type? cell type less important?

## compare those with only 4 to simplify the question
plotDat4 = plotDatBox[plotDatBox$sums == 4, ]

model4Index = which(rowSums(designMatrix) == 4)
models4 = modelList[model4Index]

wantedModelNames = sapply(strsplit(names(models4), "_"), function(x){return(x[[2]])})
names(models4) = wantedModelNames

models4Compared = ModelCompareStackedBar(bulk[[1]], 
                                         modelList = models4, 
                                         trueComparison = T,
                                         noise = F,
                                         trueProportions = bulk[[2]],
                                         nCpGPlot = F,
                                         sampleNamesOnPlots = F)




pdf(file = paste(path, "plots/cellTypesInModel/4CellTypeBoxplot.pdf", sep = ""), height = 6, width = 7)
ggplot(plotDat4, aes(x = model, y = error, fill = model))+
  geom_boxplot() +
  geom_jitter() +
  theme_cowplot(18) +
  scale_x_discrete(breaks=modelNames,
                   labels=wantedModelNames) +
  labs(x = "Model", y = "DSRMSE") +
  theme(legend.position = "none") +
  ylim(c(0, max(plotDat$error)))
dev.off()


stacks = plot_grid(models4Compared[[7]] + theme(legend.position = "none"),
          models4Compared[[3]] + theme(legend.position = "none"),
          models4Compared[[4]] + theme(legend.position = "none"),
          models4Compared[[5]] + theme(legend.position = "none"),
          models4Compared[[6]] + theme(legend.position = "none"), 
          models4Compared[[2]] + theme(legend.position = "none"),
          labels = "AUTO",
          ncol = 1)

leg = get_legend(models4Compared[[7]] + theme(legend.position=c(0.15,0.8),
                                              legend.direction = "horizontal", legend.title = element_blank()))

pdf(paste(path, "plots/cellTypesInModel/stackedNCellTypes.pdf", sep = ""), height = 16, width = 12)
plot_grid(stacks, leg, rel_heights = c(0.9,0.05), ncol = 1)
dev.off()



## plot if the sum of CD4 and 8 are the same in models without them
C4C8bulk = data.frame(bulk[[2]])
C4C8bulk$CD4_CD8 = C4C8bulk$CD4.T_cells + C4C8bulk$CD8.T_cells
C4C8bulk = C4C8bulk[,-c(2,3)]
C4C8bulk = C4C8bulk[,c(1,4,2,3)]

CD4Model = plotDat4[which(plotDat4$model == "model_BC4GM"),c(1, 2, 3, 11, 12, 14, 15)]
CD4Model = CD4Model[match(rownames(C4C8bulk), CD4Model$sample),-3]
CD4Plot = TrueVSPredictedMaintainCellsPlot(CD4Model, C4C8bulk, c(T,T,F,T,T), ylabel = "Actual (CD4 = CD4 + CD8)")

CD8Model = plotDat4[which(plotDat4$model == "model_BC8GM"),c(1, 2, 3, 11, 13, 14, 15)]
CD8Model = CD8Model[match(rownames(C4C8bulk), CD8Model$sample),-3]
CD8Plot = TrueVSPredictedMaintainCellsPlot(CD8Model, C4C8bulk, c(T,F,T,T,T), ylabel = "Actual (CD8 = CD4 + CD8)")

## plot if the sum of CD4 and 8 are the same in models without them
MGbulk = data.frame(bulk[[2]])
MGbulk$M_G = MGbulk$Monocytes + MGbulk$Granulocytes
MGbulk = MGbulk[,-c(4,5)]

MModel = plotDat4[which(plotDat4$model == "model_BC4C8M"),c(1, 2, 3, 11, 12, 13, 15)]
MModel = MModel[match(rownames(MGbulk), MModel$sample),-3]
MPlot = TrueVSPredictedMaintainCellsPlot(MModel, MGbulk, c(T,T,T,F,T), ylabel = "Actual (M = M + G)")

GModel = plotDat4[which(plotDat4$model == "model_BC4C8G"),c(1, 2, 3, 11, 12, 13, 14)]
GModel = GModel[match(rownames(MGbulk), GModel$sample),-3]
GPlot = TrueVSPredictedMaintainCellsPlot(GModel, MGbulk, c(T,T,T,T,F), ylabel = "Actual (G = M + G)")


missingCellPlot = plot_grid(CD4Plot + ggtitle("CD8 missing") + theme(legend.position = "none"), 
          CD8Plot + ggtitle("CD4 missing") + theme(legend.position = "none"), 
          MPlot + ggtitle("Granulocytes missing") + theme(legend.position = "none"), 
          GPlot + ggtitle("Monocytes missing") + theme(legend.position = "none"), 
          ncol = 2, 
          labels = c("Ai", "Aii", "Bi", "Bii"))

## make legend 
plotForLeg = ggplot(x, aes(x = error, y = proportion_pred)) +
  geom_point(aes(shape = as.factor(cellType), col = as.factor(cellType))) +
  theme_cowplot(18)
leg = get_legend(plotForLeg + theme(legend.direction = "horizontal", legend.title = element_blank())) 

pdf(paste(path, "plots/cellTypesInModel/cellTypeReplacement4cell.pdf", sep = ""), height = 8, width = 9)
plot_grid(missingCellPlot, leg, ncol = 1, rel_heights = c(1,0.05))
dev.off()

## compare true proportion of B cell to error in model without B cells
BModel = plotDat4[which(plotDat4$model == "model_C4C8GM"),]
BModel = BModel[match(rownames(bulk[[2]]), BModel$sample),]
BModelCompare = cbind.data.frame(error = BModel$error, B_cells = as.data.frame(bulk[[2]])$B_cells)

pdf(paste(path, "plots/cellTypesInModel/cellTypeBModel4cell.pdf", sep = ""), height = 6, width = 6)
ggplot(BModelCompare, aes(x = B_cells, y = error)) +
  geom_point(size = 2) +
  theme_cowplot(18) +
  labs(x = "Proportion of B cells", y = "DSRMSE")
dev.off()



### Compare error when data contains noise ############
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))
load(paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))

## create model with all 5 cell types and 150 CpGs
model = pickCompProbes(rawbetas = as.matrix(betas),
                       cellTypes = levels(as.factor(pheno$Sample.Type)),
                       cellInd = as.factor(pheno$Sample.Type),
                       numProbes =  150,
                       probeSelect = "auto")

save(model, file = paste(path, "data/ERiskModel5CellTypes150CpG.Rdata", sep = ""))

## create simulated samples with increasing noise
testData = CellTypeProportionSimulator(betas = betastest, 
                                       pheno = phenotest, 
                                       phenoColName = "Sample.Type", 
                                       nBulk = 9, 
                                       proportionsMatrix = "random",
                                       noiseIn = T,
                                       proportionNoise = seq(0,0.4,0.05))

stackedWithNoise = ModelCompareStackedBar(testBetas = testData[[1]], 
                                  modelList = list(Predicted = model), 
                                  trueComparison = T,
                                  noise = T,
                                  trueProportions = testData[[2]],
                                  nCpGPlot = F,
                                  sampleNamesOnPlots = F)


plot_grid(stackedWithNoise[[1]] + theme(legend.position = "none"), 
          stackedWithNoise[[2]], 
          stackedWithNoise[[3]], 
          ncol = 1, 
          align = "v", axis = "lr", 
          rel_heights = c(0.4,1,1))


## repeate with more samples and bigger range of error
testData = CellTypeProportionSimulator(betas = betastest, 
                                       pheno = phenotest, 
                                       phenoColName = "Sample.Type", 
                                       nBulk = length(seq(0,1,0.02)), 
                                       proportionsMatrix = "random",
                                       noiseIn = T,
                                       proportionNoise = seq(0,1,0.02))

stackedWithNoise = ModelCompareStackedBar(testBetas = testData[[1]], 
                                          modelList = list(Predicted = model), 
                                          trueComparison = T,
                                          noise = T,
                                          trueProportions = testData[[2]],
                                          nCpGPlot = F,
                                          sampleNamesOnPlots = F)

## plot noise vs error
NEDat = cbind.data.frame(error = stackedWithNoise[[1]]$data$error, noise = seq(0,1,0.02))

ggplot(NEDat, aes(x = noise, y = error)) +
  geom_point(size = 2) +
  theme_cowplot(18) +
  labs(x = "Proportion of noise", y = "DSRMSE")


### Compare whole blood, simulated and buccal predictions of blood cell types ###### 
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))
load(paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))
load("/mnt/data1/Eilis/Projects/Asthma/CellTypeComparisons/Correlations_New/BetasSortedByCellType_NoY.rdat")

## create model with all 5 cell types and 150 CpGs
model = pickCompProbes(rawbetas = as.matrix(betas),
                       cellTypes = levels(as.factor(pheno$Sample.Type)),
                       cellInd = as.factor(pheno$Sample.Type),
                       numProbes =  150,
                       probeSelect = "auto")

## create simulated samples with increasing noise
testData = CellTypeProportionSimulator(betas = betastest, 
                                       pheno = phenotest, 
                                       phenoColName = "Sample.Type", 
                                       nBulk = 9, 
                                       proportionsMatrix = "random",
                                       noiseIn = F)

simStackedBar = ModelCompareStackedBar(testBetas = testData[[1]], 
                                     modelList = list(Simulated = model), 
                                     trueComparison = T,
                                     noise = F,
                                     trueProportions = testData[[2]],
                                     nCpGPlot = F,
                                     sampleNamesOnPlots = F)


nasalStackedBar = ModelCompareStackedBar(testBetas = allbetas[[8]][,which(!is.na(colnames(allbetas[[8]])))], 
                                          modelList = list(Nasal = model), 
                                          trueComparison = F,
                                          noise = F,
                                          trueProportions = NA,
                                          nCpGPlot = F,
                                          sampleNamesOnPlots = F)


wholeBloodStackedBar = ModelCompareStackedBar(testBetas = allbetas[[5]][,which(!is.na(colnames(allbetas[[5]])))], 
                                             modelList = list(Wholeblood = model), 
                                             trueComparison = F,
                                             noise = F,
                                             trueProportions = NA,
                                             nCpGPlot = F,
                                             sampleNamesOnPlots = F)

buccalStackedBar = ModelCompareStackedBar(testBetas = allbetas[[1]][,which(!is.na(colnames(allbetas[[1]])))], 
                                          modelList = list(Buccal = model), 
                                          trueComparison = F,
                                          noise = F,
                                          trueProportions = NA,
                                          nCpGPlot = F,
                                          sampleNamesOnPlots = F)

## plot error against cell type 
plotDat = rbind.data.frame(cbind.data.frame(error = simStackedBar[[1]]$data$error[1:9], cellType = "Simulated\nblood"),
                           cbind.data.frame(error = wholeBloodStackedBar[[1]]$data$error[1:28], cellType = "Whole\nblood"),
                           cbind.data.frame(error = buccalStackedBar[[1]]$data$error[1:28], cellType = "Buccal"),
                           cbind.data.frame(error = nasalStackedBar[[1]]$data$error[1:19], cellType = "Nasal"))

ggplot(plotDat, aes(x = cellType, y = error, fill = cellType)) +                           
  geom_violin() +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  labs(x = element_blank(), y = "DSRMSE")

### Check error with increased missingness of test data ####
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))
load(paste(path, "data/ERiskModel5CellTypes150CpG.Rdata", sep = ""))

## subset betas to include only those in the model
subBetas = GetModelCG(betastest, list(model))


## create function to add x NAs to betas
MakeNAsInBetas = function(proportionNA, betas){
  nToBeNA = floor(nrow(betas)*proportionNA)
  NAIndex = sample(nrow(betas), nToBeNA)
  betas[NAIndex,] = NA
  return(betas)
}

propMissing = seq(0,0.9,0.05)
x = lapply(propMissing, MakeNAsInBetas, subBetas)

plotDat = ErrorAcrossDataSets(x, model)

ggplot(plotDat, aes(x = as.factor(nCGmissing), y = error)) +
  geom_boxplot() +
  theme_cowplot(18) +
  scale_x_discrete(labels = propMissing) +
  labs(x = "Proportion of CpGs missing", y = "DSRMSE")
  


