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

t.test(IoPPred[sample.annot$Euro,"error"], IoPPred[!sample.annot$Euro,"error"])
t.test(EUGEIPred[pheno$Euro,"error"], EUGEIPred[!pheno$Euro,"error"])
