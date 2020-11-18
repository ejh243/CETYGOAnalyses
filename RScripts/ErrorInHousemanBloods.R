## Script to test error metric using Houseman data set
## Dorothea Seiler Vellame
## started 18-11-2020

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

## split into training and testing and add to phenotype cols
phenoTrain = c(rep("Test", length(pheno)))
for(i in 1:length(levels(pheno))){
  cellTypeIndex = which(pheno == levels(pheno)[i])
  set.seed(123)
  phenoTrain[sample(cellTypeIndex, ceiling(length(cellTypeIndex)/2))] = "Train"
}


## save data
betasTrain = betas[,phenoTrain == "Train"]
betasTest = betas[,phenoTrain == "Test"]
pheno = data.frame(celltype = pheno, trainTest = phenoTrain)

save(betas, betasTrain, betasTest, pheno, file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTest.Rdata")

## normalise data as train test for each normalisation method
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
save(quantileBetasTrain, quantileBetasTest, quantileBetas, 
     file = "/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTest.Rdata")


