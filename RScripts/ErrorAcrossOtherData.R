## Test the blood deconvolution model created in Erisk in 
## Understanding society and EXTEND to get normal range in 
## error metric

## Dorothea Seiler Vellame 05-11-2020

### load data EXTEND and Understanding society ########
load("/mnt/data1/EPICQC/UnderstandingSociety/US_Betas_Pheno.rda")
us = dat
load("/mnt/data1/EXTEND/Methylation/QC/EXTEND_batch1_2_merged/EXTEND_batches_1_2_normalised_together.rdat")
ex = betas

rm(dat, pheno, betas)

## load model and functions
path = "/mnt/data1/Thea/ErrorMetric/"
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
load(paste(path, "data/ERiskModel5CellTypes150CpG.Rdata", sep = ""))

## get error for US
usPred = projectCellTypeWithError(GetModelCG(us, list(modelCG = model)) , model$coefEsts)

## get error for EX
exPred = projectCellTypeWithError(GetModelCG(ex, list(modelCG = model)) , model$coefEsts)

## plot!
library(ggplot2)
library(cowplot)

plotDatBC = rbind.data.frame(cbind.data.frame(usPred, dataSet = "Blood\nUnderstanding Society"),
                           cbind.data.frame(exPred, dataSet = "Blood\nEXTEND"))

# ggplot(plotDatBC aes(x = dataSet, y = error, fill = dataSet)) +
#   geom_violin() +
#   theme_cowplot(18) +
#   labs(x = "Dataset", y = "Error") +
#   theme(legend.position = "none") +
#   ylim(c(0,max(plotDat$error)))


### Ailsa's pancreas data #############################
load("/mnt/data1/EPICQC/Ailsa/QC/Ailsa_normalised.rdat") 
ap = betas.dasen
rm(betas.dasen)

apPred = projectCellTypeWithError(GetModelCG(ap, list(modelCG = model)) , model$coefEsts)

plotDatAP = rbind.data.frame(cbind.data.frame(apPred, dataSet = "Pancreas"))

# ggplot(plotDat, aes(x = dataSet, y = error, fill = dataSet)) +
#   geom_violin() +
#   theme_cowplot(18) +
#   labs(x = "Dataset", y = "Error") +
#   theme(legend.position = "none") +
#   ylim(c(0,max(plotDat$error)))

### Horvaths cancer blood data ########################
library(data.table)
cb = fread("/mnt/data1/Thea/ErrorMetric/data/externalData/GSE140038_NormalizedBetaNoob.csv")
cb = data.frame(cb)
rownames(cb) = cb[,1]
cb = cb[,-1]

## load pheno
cbPheno = read.delim("/mnt/data1/Thea/ErrorMetric/data/externalData/GSE140038_series_matrix.txt",sep = " ", skip = 14)

cbPred = projectCellTypeWithError(GetModelCG(cb, list(modelCG = model)) , model$coefEsts)

plotDatCB = rbind.data.frame(cbind.data.frame(cbPred, dataSet = "Blood with\nbreast cancer"))

# ggplot(plotDat, aes(x = dataSet, y = error, fill = dataSet)) +
  # geom_violin() +
  # theme_cowplot(18) +
  # labs(x = "Dataset", y = "Error") +
  # theme(legend.position = "none") +
  # ylim(c(0,max(plotDat$error)))


### PPMI ##############################################
load("/mnt/data1/Josh/PPMI/Methylation_data/Pipeline/finalNormBetas.Rdata")
pd = betas
rm(betas)

## load pheno 
pdPheno = read.csv("/mnt/data1/Josh/PPMI/Methylation_data/QCmetrics.csv")

## check that sample orders match
all(pdPheno$Basename == colnames(pd))

##subset for samples with pd
pd = pd[,pdPheno$ENROLL_CAT == "PD"]
pdPred = projectCellTypeWithError(GetModelCG(pd, list(modelCG = model)) , model$coefEsts)

plotDatPD = rbind.data.frame(cbind.data.frame(pdPred, dataSet = "Blood with\nParkinsons"))

# ggplot(plotDat, aes(x = dataSet, y = error, fill = dataSet)) +
#   geom_violin() +
#   theme_cowplot(18) +
#   labs(x = "Dataset", y = "Error") +
#   theme(legend.position = "none") +
#   ylim(c(0,max(plotDat$error)))



### E-risk non cell type ##############################
path = "/mnt/data1/Thea/ErrorMetric/"
source(paste(path, "RScripts/FunctionsForBrainCellProportionPrediction.r", sep = ""))
source(paste(path, "RScripts/FunctionsForErrorTesting.R", sep = ""))
# load(paste(path, "data/bloodEriskDataTrain.Rdata", sep = ""))
load(paste(path, "data/bloodEriskDatatest.Rdata", sep = ""))
load("/mnt/data1/Eilis/Projects/Asthma/CellTypeComparisons/Correlations_New/BetasSortedByCellType_NoY.rdat")
load(paste(path, "data/ERiskModel5CellTypes150CpG.Rdata", sep = ""))

ebuccPred = projectCellTypeWithError(GetModelCG(allbetas[[1]][,!is.na(colnames(allbetas[[1]]))], list(modelCG = model)) , model$coefEsts)
ebloodPred = projectCellTypeWithError(GetModelCG(allbetas[[6]][,!is.na(colnames(allbetas[[6]]))], list(modelCG = model)) , model$coefEsts)
enasalPred = projectCellTypeWithError(GetModelCG(allbetas[[8]][,!is.na(colnames(allbetas[[8]]))], list(modelCG = model)) , model$coefEsts)

simBlood = CellTypeProportionSimulator(betas = betastest, 
                                       pheno = phenotest, 
                                       phenoColName = "Sample.Type", 
                                       nBulk = 25, 
                                       proportionsMatrix = "random",
                                       noiseIn = F)

simPred = projectCellTypeWithError(GetModelCG(simBlood[[1]], list(modelCG = model)) , model$coefEsts)

plotDatErisk = rbind.data.frame(cbind.data.frame(ebuccPred, dataSet = "Buccal"),
                           cbind.data.frame(ebloodPred, dataSet = "Blood\nE-Risk"),
                           cbind.data.frame(enasalPred, dataSet = "Nasal"),
                           cbind.data.frame(simPred, dataSet = "Simulated\nblood"))

# ggplot(plotDat, aes(x = dataSet, y = error, fill = dataSet)) +
#   geom_violin() +
#   theme_cowplot(18) +
#   labs(x = "Dataset", y = "Error") +
#   theme(legend.position = "none") +
#   ylim(c(0,max(plotDat$error)))


### combine all data
plotDat = rbind.data.frame(plotDatErisk, 
                           plotDatPD,
                           plotDatCB,
                           plotDatAP,
                           plotDatBC)
plotDat$dataSet = factor(plotDat$dataSet, levels = c("Simulated\nblood", 
                                                     "Blood\nE-Risk",
                                                     "Blood\nUnderstanding Society",
                                                     "Blood\nEXTEND",
                                                     "Blood with\nParkinsons",
                                                     "Blood with\nbreast cancer",
                                                     "Buccal",
                                                     "Nasal",
                                                     "Pancreas"))

pdf("/mnt/data1/Thea/ErrorMetric/plots/testModelAcrossCellTypes/allDataError.pdf", height = 7, width = 16)
ggplot(plotDat, aes(x = dataSet, y = error, fill = dataSet)) +
  geom_violin() +
  theme_cowplot(18) +
  labs(x = "Dataset", y = "Error") +
  theme(legend.position = "none") +
  ylim(c(0,max(plotDat$error)))
dev.off()

save(plotDat, file = "/mnt/data1/Thea/ErrorMetric/data/allDataTypesPlottedError.Rdata")

