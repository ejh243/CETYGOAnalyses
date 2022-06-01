#### Error across ethnicity ####

### load the datasets ####

library(wateRmelon)
# setwd("/mnt/data1/Schizophrenia/Blood/IoP")
load("/mnt/data1/Schizophrenia/Blood/IoP/FinalQCedNormalised_IOP_Filter_PC1_PC2_2SDFromMean.Rdata")
betas<-betas(mset450k.pf.dasen)
sample.annot<-sample.annot[match(colnames(betas), sample.annot$Basename),]
euro <- read.table("/mnt/data1/Schizophrenia/Blood/IoP/Genotypes/EuropeanSampleList.txt")
sample.annot$Euro <- sample.annot$V1 %in% euro$V1

source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

IoPPred = projectCellTypeWithError(betas, modelType = "ownModel", ownModelData = HousemanBlood50CpGModel)

plotDat = cbind.data.frame(Cetygo = IoPPred[,"error"],
                           Euro = sample.annot$Euro,
                           Dataset = "IoP")

# rm(list=setdiff(ls(),c("plotDat", "sample.annot", "IoPPred")))

load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Blood_WithRepeats/JustEuGEIresults/EuGEIBloodSamples_Normalised.rdat")
euro <- read.table("/mnt/data1/EuGEI/Genotypes/EUGEI_GROUP_EuroOnly_PCs_NonRefPanel_Dec17.txt", header = T)
pheno$Euro <- pheno$Geno.CHIP.Location %in% euro$FID 

source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")

EUGEIPred = projectCellTypeWithError(betas, modelType = "ownModel", ownModelData = HousemanBlood50CpGModel)

plotDat = rbind.data.frame(plotDat,
                           cbind.data.frame(Cetygo = EUGEIPred[,"error"],
                                            Euro = pheno$Euro,
                                            Dataset = "EUGEI"))
# rm(list=setdiff(ls(),c("plotDat","sample.annot", "pheno", "IoPPred", "EUGEIPred")))
library(plyr)
plotDat$Euro = revalue(as.factor(plotDat$Euro), c("TRUE" = "Yes", "FALSE" = "No"))
plotDat$Euro = factor(plotDat$Euro, levels = c("Yes", "No"))


save(plotDat, file = "/mnt/data1/Thea/ErrorMetric/data/plotDat/EthnicityData.Rdata") #?

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/EthnicityEUGEUIoPBloos.pdf", height = 5, width = 6)
ggplot(plotDat, aes(x = Euro, y = Cetygo)) +
  geom_violin() +
  facet_wrap(~Dataset) +
  theme_cowplot(18) +
  labs(x = "European?") +
  ylim(c(0,max(plotDat$Cetygo))) +
  geom_hline(yintercept = 0.1, linetype = "dashed", col = "red")
dev.off()

t.test(IoPPred[sample.annot$Euro,"error"], IoPPred[!sample.annot$Euro,"error"])$p.value
t.test(EUGEIPred[pheno$Euro,"error"], EUGEIPred[!pheno$Euro,"error"])



### REdo age sex plots with EUGEI and IoP too 
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
load("/mnt/data1/EuGEI/QC/GeorginasQC/All_Plates_Blood_WithRepeats/JustEuGEIresults/EuGEIBloodSamples_Normalised.rdat")
eu = betas
euPheno = pheno
load("/mnt/data1/Schizophrenia/Blood/IoP/FinalQCedNormalised_IOP_Filter_PC1_PC2_2SDFromMean.Rdata")
library(wateRmelon)
iop<-betas(mset450k.pf.dasen)
iopPheno<-sample.annot[match(colnames(iop), sample.annot$Basename),]

rm(dat, betas, pheno, HousemanBlood50CpGModel, sample.annot)

## get error for US
usPred = projectCellTypeWithError(us, modelType = "ownModel", ownModelData = model)

## get error for EX
exPred = projectCellTypeWithError(ex, modelType = "ownModel", ownModelData = model)

## get error for eu
euPred = projectCellTypeWithError(eu, modelType = "ownModel", ownModelData = model)

## get error for iop
iopPred = projectCellTypeWithError(iop, modelType = "ownModel", ownModelData = model)

## format sex and age

sexUS = usPheno$nsex
sexUS = ifelse(sexUS =="1", "Male", "Female")
sexIoP = iopPheno$Gender
sexIoP = ifelse(sexIoP =="M", "Male", "Female")
ageIoP = iopPheno$AgeAtBloodCollection
ageIoP[ageIoP == -99] = NA
allPheno = rbind.data.frame(cbind.data.frame(usPred, data = "Understanding\nSociety", age = usPheno$confage, sex = sexUS),
                            cbind.data.frame(exPred, data = "EXTEND", age = exPheno$Age, sex = exPheno$Sex), 
                            cbind.data.frame(euPred, data = "EUGEI", age = euPheno$Age, sex = euPheno$Sex),
                            cbind.data.frame(iopPred, data = "IoP", age = ageIoP, sex = sexIoP))

save(allPheno, file = "/mnt/data1/Thea/ErrorMetric/data/plotDat/SexAgePlotDat.Rdata")

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/sexAcrossEXandUS.pdf", height = 6, width = 7) 
ggplot(allPheno, aes(x = data, y = error, fill = sex)) +
  geom_violin() +
  geom_hline(yintercept = 0.1, col = "red", linetype = "dashed") +
  theme_cowplot(18) +
  theme(legend.position = "bottom", legend.justification = "center") +
  labs(y = "Cetygo", x = "Dataset", fill = "Sex")
dev.off()

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/ageAcrossEXandUS.pdf", height = 6, width = 7) 
ggplot(allPheno, aes(x = as.numeric(as.character(age)), y = error, col = data)) +
  geom_point() +
  geom_hline(yintercept = 0.1, col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(38-13.6, 38+13.6), linetype = "dashed") +
  geom_vline(xintercept = 38) +
  # scale_color_manual(values = colToUse) +
  theme_cowplot(18) +
  theme(legend.position = "bottom", legend.justification = "center") +
  labs(y = "Cetygo", x = "Age", col = NULL)
dev.off()

load("/mnt/data1/Thea/ErrorMetric/data/plotDat/SexAgePlotDat.Rdata")

x = t.test(allPheno[allPheno$sex == "Female","error"], allPheno[allPheno$sex == "Male","error"], alternative = "greater")
x=summary(lm(error~as.numeric(age)+sex, allPheno))

x = t.test(allPheno[which(as.numeric(allPheno$age) < 38-13.6 | as.numeric(allPheno$age) > 38+13.6),"error"], 
           allPheno[which(as.numeric(allPheno$age) > 38-13.6 | as.numeric(allPheno$age) < 38+13.6),"error"])


## X chromosome plots
load("/mnt/data1/Thea/ErrorMetric/data/cpgInModel.Rdata")

exT = ex[rownames(ex) %in% cpgInMod$IlmnID,]
usT = us[rownames(us) %in% cpgInMod$IlmnID,]
euT = eu[rownames(eu) %in% cpgInMod$IlmnID,]
iopT = iop[rownames(iop) %in% cpgInMod$IlmnID,]

exT = exT[match(cpgInMod$IlmnID, rownames(exT)),]
usT = usT[match(cpgInMod$IlmnID, rownames(usT)),]
euT = euT[match(cpgInMod$IlmnID, rownames(euT)),]
iopT = iopT[match(cpgInMod$IlmnID, rownames(iopT)),]

all(rownames(exT) == cpgInMod$IlmnID)
all(rownames(usT) == cpgInMod$IlmnID)
all(rownames(euT) == cpgInMod$IlmnID)
all(rownames(iopT) == cpgInMod$IlmnID)

## subset for those in the x chromosome
exT = exT[cpgInMod$CHR == "X",]
usT = usT[cpgInMod$CHR == "X",]
euT = euT[cpgInMod$CHR == "X",]
iopT = iopT[cpgInMod$CHR == "X",]

# sexUS = usPheno$nsex
# sexUS = ifelse(sexUS =="1", "Male", "Female")

plotDat = rbind.data.frame(data.frame(t(usT), sampleID = colnames(usT), study = "Understanding\nSociety", sex = sexUS),
                           data.frame(t(exT), sampleID = colnames(exT), study = "EXTEND", sex = exPheno$Sex),
                           data.frame(t(euT), sampleID = colnames(euT), study = "EUGEI", sex = euPheno$Sex),
                           data.frame(t(iopT), sampleID = colnames(iopT), study = "IoP", sex = sexIoP))

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
print(xtable(statDat, display = c("s","g","s","s","s","s")),
      math.style.exponents = TRUE, include.rownames=FALSE)

library(reshape2)
plotDat = melt(plotDat, id.vars = c("study", "sex", "sampleID"))

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/ErrorXchrCpGsInEXUS.pdf", height = 9, width = 7)
ggplot(plotDat, aes(x = variable, y = value, fill = sex)) +
  theme_cowplot(18) +
  geom_violin(position=position_dodge(0.5)) +
  facet_wrap(~study, ncol = 1) +
  labs(y = "Proportion of methylation", fill = "Sex", x = "X chromosome CpGs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15),
        legend.position = "bottom", legend.justification = "center")
dev.off()

#### batch plot ###
load("/mnt/data1/Thea/ErrorMetric/data/plotDat/SexAgePlotDat.Rdata")

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/modelApplicability/batchUSEXUEIOP.pdf", height = 5, width = 6)
ggplot(allPheno, aes(x = data, y = error)) +
  geom_hline(yintercept = 0.1, col = "red", linetype = "dashed") +
  geom_violin()+
  theme_cowplot(18) +
  theme(legend.position = "none") +
  labs(x = "Dataset", y = "Cetygo") 
dev.off()

# t.test(allPheno[allPheno$data == "Understanding\nSociety", "error"], 
#        allPheno[allPheno$data != "Understanding\nSociety", "error"])$p.value
# 
# t.test(allPheno[allPheno$data == "EXTEND", "error"], 
#        allPheno[allPheno$data != "EXTEND", "error"])$p.value
# 
# t.test(allPheno[allPheno$data == "EUGEI", "error"], 
#        allPheno[allPheno$data != "EUGEI", "error"])$p.value
# 
# t.test(allPheno[allPheno$data == "IoP", "error"], 
#        allPheno[allPheno$data != "IoP", "error"])$p.value

