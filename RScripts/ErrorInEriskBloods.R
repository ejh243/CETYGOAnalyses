## Script to test error metric using E-risk data set
## Dorothea Seiler Vellame
## started 09-10-2020

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

pdf("/mnt/data1/Thea/ErrorMetric/plots/predictedAgeAcrossAllCellTypes.pdf", height = 5, width = 8)
ggplot(bloodphenomerge, aes(x = PatientID, y = PredictedAge, col = Sample.Type)) +
  geom_point() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Patient ID", y = "Predicted age", col = "Cell type")
dev.off()

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

pdf("/mnt/data1/Thea/ErrorMetric/plots/trainTestSamplesAcrossPredictedAgeAndSex.pdf", height = 5, width = 5)
ggplot(wholebloodphenomerge, aes(x = sex, y = PredictedAge, col = V5)) +
  geom_jitter(width = 0.3, size = 2) +
  theme_cowplot() +
  labs(x = "Sex", y = "Predicted age", col = "")
dev.off()

## make train test PCA across each cell type to show randomness
## subset 100 sites
for (i in 1:length(bloodbeta)){
  sd = apply(bloodbeta[[i]], 1, sd)
  mostVar = 
}


## make PCA of 100 most variable sites across all cells to check no mislabeling


