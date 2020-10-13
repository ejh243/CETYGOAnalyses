## Script to test error metric using E-risk data set
## Dorothea Seiler Vellame
## started 09-10-2020


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

# pdf("/mnt/data1/Thea/ErrorMetric/plots/predictedAgeAcrossAllCellTypes.pdf", height = 5, width = 8)
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

# pdf("/mnt/data1/Thea/ErrorMetric/plots/trainTestSamplesAcrossPredictedAgeAndSex.pdf", height = 5, width = 5)
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

pdf("/mnt/data1/Thea/ErrorMetric/plots/PCAPerCellTypeTrainTest.pdf", height = 6, width = 15)
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
     file = "/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataFull.Rdata")

## save train and test data
betas = bloodBeta[,intersect(which(bloodPheno$traintest == "train"), which(bloodPheno$Sample.Type != "whole blood"))]
betastest = bloodBeta[,intersect(which(bloodPheno$traintest == "test"), which(bloodPheno$Sample.Type != "whole blood"))]
  
pheno = bloodPheno[intersect(which(bloodPheno$traintest == "train"), which(bloodPheno$Sample.Type != "whole blood")),]
phenotest = bloodPheno[intersect(which(bloodPheno$traintest == "train"), which(bloodPheno$Sample.Type != "whole blood")),]
  
save(betas, pheno,
     file = "/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
save(betastest, phenotest,
     file = "/mnt/data1/Thea/ErrorMetric/data/bloodEriskDatatest.Rdata")



## Find a sufficient number of CpGs to include ########
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForBrainCellProportionPrediction.r")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")
load("/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
library(genefilter)


model = pickCompProbes(rawbetas = as.matrix(betas),
                       cellTypes = levels(as.factor(pheno$Sample.Type)),
                       cellInd = as.factor(pheno$Sample.Type),
                       numProbes =  10,
                       probeSelect = "auto")













## Create models using 2:n cell types #################
