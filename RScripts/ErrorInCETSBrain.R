### Error metric tested in brain CETS data ############
## 
## noise, train test <0.1, brain compared to other tissues


### Setup ##############################################
library(minfi)
library(wateRmelon)
library(FlowSorted.DLPFC.450k)
library(gdsfmt)
library(IlluminaHumanMethylation450kmanifest)


## CETS data
referencePkg <- sprintf("FlowSorted.%s.%s", "DLPFC", "450k")
data(list = referencePkg)
referenceRGset <- get(referencePkg)
betasCETS = getBeta(preprocessRaw(referenceRGset))
phenoCETS = data.frame(Sample_ID = colnames(betasCETS),
                       Celltype = pData(referenceRGset)$CellType,
                       Sex = pData(referenceRGset)$sex,
                       Ethnicity = pData(referenceRGset)$ethnicity,
                       Age = pData(referenceRGset)$age,
                       PMI = pData(referenceRGset)$PMI,
                       SentrixID = paste(pData(referenceRGset)$Sentrix.Barcode, pData(referenceRGset)$Sample.Section, sep = "_"),
                       SentrixPlate = pData(referenceRGset)$Sentrix.Barcode,
                       Individual = as.factor(unlist(strsplit(colnames(betasCETS), "_"))[seq(1,ncol(betasCETS)*2, 2)]))
rm(list = setdiff(ls(), c("betasCETS", "phenoCETS")))

phenoCETS$Celltype = revalue(phenoCETS$Celltype, c("NeuN_pos" = "NeuN+",
                                                   "NeuN_neg" = "NeuN-"))
phenoCETS$Celltype = as.factor(as.character(phenoCETS$Celltype))


cellTypesAll = unique(c(as.character(phenoCETS$Celltype)))
colOrdeByCell = c("#4ab9b0", "#714a91")
names(colOrdeByCell) = cellTypesAll

save(colOrdeByCell, file = "/mnt/data1/Thea/humanDeconvolution/data/CETScols.Rdata")

save(betasCETS, phenoCETS, file = "/mnt/data1/Thea/humanDeconvolution/data/CETSUnnormalised.RData")

# load("/mnt/data1/Thea/humanDeconvolution/data/CETSUnnormalised.RData")
# 
# ## n needed in test
# length(levels(phenoCETS$Individual))*.5 # 14
# 
# ## subset to NeuN+ to have individuals not samples
# temp = phenoCETS[phenoCETS$Celltype == "NeuN+",]
# tempA = temp[temp$Ethnicity == "African",]
# 
# ## randomly select 2 African samples
# set.seed(1234)
# ATest = tempA$Individual[sample(nrow(tempA), 3)]
# 
# ## randomly select 7 Caucasian samples
# temp = temp[temp$Ethnicity == "Caucasian",]
# set.seed(1234)
# CTest = temp$Individual[sample(nrow(temp), 11)]
# 
# phenoCETS$TrainTestnNeeded = "Train"
# phenoCETS$TrainTestnNeeded[phenoCETS$Individual %in% c(as.numeric(as.character(ATest)), as.numeric(as.character(CTest)))] = "Test"
# 
# betasCETSTest = betasCETS[,phenoCETS$TrainTestnNeeded == "Test"]
# phenoCETSTest = phenoCETS[phenoCETS$TrainTestnNeeded == "Test",]
# 
# betasCETSTrain = betasCETS[,phenoCETS$TrainTestnNeeded == "Train"]
# phenoCETSTrain = phenoCETS[phenoCETS$TrainTestnNeeded == "Train",]
# 
# save(betasCETSTest, phenoCETSTest,betasCETSTrain,phenoCETSTrain,
#      file = "/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
# 

### PCA for CETS cell types ###########################
## the PCA of all is in HumanBrainRBDM.R
library(ggfortify)
library(ggplot2)
library(cowplot)
library(scales)

load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/humanDeconvolution/data/cellTypeColours.Rdata")

## subset data by cell type
betasNP = cbind(betasCETSTrain[,which(phenoCETSTrain$Celltype == "NeuN+")],
                betasCETSTest[,which(phenoCETSTest$Celltype == "NeuN+")])

betasNN = cbind(betasCETSTrain[,which(phenoCETSTrain$Celltype == "NeuN-")],
                betasCETSTest[,which(phenoCETSTest$Celltype == "NeuN-")])

phenoNP = rbind(phenoCETSTrain[which(phenoCETSTrain$Celltype == "NeuN+"),],
                phenoCETSTest[which(phenoCETSTest$Celltype == "NeuN+"),])

phenoNN = rbind(phenoCETSTrain[which(phenoCETSTrain$Celltype == "NeuN-"),],
                phenoCETSTest[which(phenoCETSTest$Celltype == "NeuN-"),])

## PC plots
plotDat = list() 
betaVar = apply(betasNP, 1, var, na.rm = T)
topBeta = betasNP[order(betaVar, decreasing = T)[1:1000],]
plot = prcomp(t(topBeta[complete.cases(topBeta),])) 
plotDat[[1]] = autoplot(plot, data = phenoNP, col = "Sex", shape = "TrainTest", size = 3) + 
  scale_shape_manual(values = c(21, 19)) +
  theme_cowplot(18) + 
  ggtitle("NeuN+") +
  labs(shape = "", col = "")


betaVar1 = apply(betasNN, 1, var, na.rm = T)
topBeta1 = betasNN[order(betaVar1, decreasing = T)[1:1000],]
plot1 = prcomp(t(topBeta1[complete.cases(topBeta1),])) 
plotDat[[2]] = autoplot(plot1, data = phenoNN, col = "Sex", shape = "TrainTest", size = 3) + 
  scale_shape_manual(values = c(21, 19)) +
  theme_cowplot(18) + 
  ggtitle("NeuN-") +
  theme(legend.position = "none")


leg = get_legend(plotDat[[1]]+ theme(legend.justification="center" ,legend.position = "bottom"))

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSCelltypePCA.pdf", height = 10, width = 5)
plot_grid(plotDat[[1]]+theme(legend.position = "none"),
          plotDat[[2]],
          leg, labels = c("A", "B", ""),
          ncol = 1,
          rel_heights = c(1,1,0.2))
dev.off()

# ### create model using test data ######################
# load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
# source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
# 
# CETSmodel = pickCompProbes(rawbetas = betasCETSTrain,
#                            cellTypes = levels(as.factor(as.character(phenoCETSTrain$Celltype))),
#                            cellInd = as.factor(as.character(phenoCETSTrain$Celltype)),
#                            numProbes = 50,
#                            probeSelect = "auto")
# 
# save(CETSmodel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")

### heatmap of model CpGs #############################
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")

## load data
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/humanDeconvolution/data/cellTypeColours.Rdata")

## load functions
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

modelBetas = GetModelCG(betasCETSTrain, list(CETSmodel))
library(gplots)
library(viridis)
library(scales)
library(ComplexHeatmap)

col = list(Celltype = colOrdeByCell[c(1,3)])

# Create the heatmap annotation
ha <- HeatmapAnnotation(Celltype = phenoCETSTrain$Celltype,
                        col = col)

pdf("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/heatmapForModelCpGsCETS.pdf", height = 6, width = 5)
Heatmap(modelBetas, name = "DNAm",
        top_annotation = ha, show_row_names = F, show_column_names = F, show_row_dend = F)
dev.off()

### check annotation of model CpGs to chromosomes #####
# x = read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7, header = T)
x = load("/mnt/data1/450K_reference/AllProbeIlluminaAnno.Rdata")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")

x = probeAnnot[,c("ILMNID", "CHR")]
cpgInMod = x[x$ILMNID %in% rownames(CETSmodel$coefEsts),]
table(cpgInMod$CHR)

### Compare sex application ###########################
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/cellTypeColours.Rdata")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")

pred = as.data.frame(projectCellTypeWithError(betasCETSTest, modelType = "ownModel", ownModelData = CETSmodel))
pred = cbind.data.frame(pred, phenoCETSTest)

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSSexViolin.pdf", height = 5, width = 5)
ggplot(pred, aes(x = Sex, y = error, fill = Sex)) +
  geom_violin() +
  theme_cowplot(18) +
  theme(legend.position = "none") +
  labs(y = "Cetygo") 
dev.off()

t.test(pred$error[pred$Sex == "Female"], pred$error[pred$Sex != "Female"])

### compare error with noise ##########################
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/cellTypeColours.Rdata")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

propCells = matrix(nrow = 1, byrow = T, data = c(0.5,0.5))
colnames(propCells) = levels(phenoCETSTrain$Celltype)
noise = seq(0,0.95,0.05)

## create simulated samples with increasing noise
testData = CellTypeProportionSimulator(betas = betasCETSTest, 
                                       pheno = phenoCETSTest, 
                                       phenoColName = "Celltype", 
                                       nBulk = length(noise), 
                                       proportionsMatrixType = "own",
                                       proportionsMatrix = propCells,
                                       noiseIn = T,
                                       proportionNoise = noise)

stackedWithNoise = ModelCompareStackedBar(testBetas = testData[[1]], 
                                          modelList = list(Predicted = CETSmodel), 
                                          trueComparison = T,
                                          noise = T,
                                          trueProportions = testData[[2]],
                                          nCpGPlot = F,
                                          sampleNamesOnPlots = F)


leg = get_legend(stackedWithNoise[[3]] + scale_fill_manual(values = c(colOrdeByCell[c(3,1)], "#555353"))
                 + theme(legend.justification="center" ,legend.direction = "horizontal", legend.title = element_blank())
                 + guides(fill = guide_legend(nrow = 1)))
plots = plot_grid(stackedWithNoise[[1]] + theme(legend.position = "none"),
                  stackedWithNoise[[2]] + theme(legend.position = "none") +
                    scale_y_continuous(breaks = seq(0,1,0.25))+
                    scale_fill_manual(values = colOrdeByCell[c(3,1)]),
                  stackedWithNoise[[3]] + theme(legend.position = "none") + scale_fill_manual(values = c(colOrdeByCell[c(3,1)], "#555353")), ncol = 1,
                  rel_heights = c(0.6,1,1), labels = "AUTO", axis = "rl", align = "v" )

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSsimWithNoise.pdf", height = 9, width = 7.5)     
print(plot_grid(plots, leg, ncol = 1, rel_heights = c(1,0.08)))
dev.off()


### Apply to real data #######
## load subset Essex data #######
# ## IN ESSEX SERVER
# 
# ## load gds
# library(gdsfmt)
# x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)
# 
# load("~/DSRMSE/models/CETSmodel50CpG.Rdata")
# 
# 
# betas = read.gdsn(index.gdsn(x, "rawbetas"))
# 
# pID = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
# betaIndex = pID %in% rownames(CETSmodel$coefEsts)
# 
# rownames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
# colnames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "FullBarCode"))
# 
# betaMod = betas[betaIndex,]
# betaMod = betaMod[match(rownames(CETSmodel$coefEsts), rownames(betaMod)),]
# 
# ## create column index, removing samples that have "", or unsorted
# library(bigmelon)
# pData = pData(x)
# 
# tDat = pData$Tissue
# 
# tInd = which(tDat == "" |
#                tDat == "Unsorted Tissues" |
#                tDat == "Unsorted Cell Line" |
#                tDat == "Unsorted Tumours")
# pData = pData[-tInd,]
# 
# betaMod = betaMod[,-tInd]
# 
# ## filter for samples with all 100 CpGs
# incompInd = complete.cases(t(betaMod))
# 
# betaMod = betaMod[,incompInd]
# pData = pData[incompInd,]
# 
# gfile = createfn.gds("brain.gds")
# 
# add.gdsn(gfile, "pData", pData)
# add.gdsn(gfile, "rownames", rownames(betaMod))
# add.gdsn(gfile, "colnames", colnames(betaMod))
# 
# 
# #source("~/DSRMSE/FunctionsForErrorTesting.R")
# source("~/DSRMSE/projectCellTypeWithError.R")
# model = CETSmodel
# 
# ind = c(seq(1000, 19129, 1000))
# errPred = lapply(ind, function(ind){
#   pred = projectCellTypeWithError(betaMod[,(ind - 999):ind], model = "ownModel", ownModelData = model)
#   return(pred)})
# 
# errPred[[length(errPred)+1]] = projectCellTypeWithError(betaMod[,19001:19131], model = "ownModel", ownModelData = model)
# 
# ## extract all from list and add to gds file
# pred = do.call("rbind", errPred)
# 
# add.gdsn(gfile, "Pred", pred)
# closefn.gds(gfile)
# closefn.gds(x)
# q() ######

## load
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/CETScols.Rdata")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")
library(gdsfmt)
library(bigmelon)
gfile = openfn.gds("/mnt/data1/Thea/ErrorMetric/data/EssexOutput/brain.gds")

dat = cbind.data.frame(pData(gfile), read.gdsn(index.gdsn(gfile$root, "Pred")))
colnames(dat)[11:14] = c("NeuN-", "NeuN+", "error", "nCGmissing")

## add a brain binary
dat$brain = rep(0, nrow(dat))
dat$brain[dat$Tissue == "Brain" |
            dat$Tissue == "Neuron"] = 1

dat$TissueBrain = dat$Tissue
dat$TissueBrain[dat$brain == 1] = "Brain"

closefn.gds(gfile)

## plot
library(ggplot2)
library(cowplot)
library(forcats)
library(dplyr)

dat_summary = dat %>%
  group_by(TissueBrain) %>%
  tally()

dat = merge(dat, dat_summary,  by = "TissueBrain")

dat$TissueBrain = as.factor(as.character(dat$TissueBrain))
pos = c()
for(i in 1:length(levels(dat$TissueBrain))){
  pos = c(pos, max(dat$error[dat$Tissue == levels(dat$TissueBrain)[i]]))
}

dat.pos = data.frame(TissueBrain = levels(dat$TissueBrain), pos, n = dat_summary$n)

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/ErrorCETSEssexsAllTissueBoxplot.pdf", height = 9, width = 14)
ggplot(dat, aes(x = fct_reorder(TissueBrain, brain, .fun = median, .desc =TRUE))) +
  geom_boxplot(aes(y = error, fill = as.factor(brain))) +
  theme_cowplot(18) +
  scale_fill_manual(values = c("#0A8ABA", "#BA3A0A"), name = "Brain?", labels = c("No", "Yes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = element_blank(), y = "Cetygo") +
  geom_text(data = dat.pos, aes(TissueBrain, label = n, y = pos+0.02))
dev.off()

## t test between brain and non brain samples
t.test(dat$error[dat$brain ==0], dat$error[dat$brain ==1])



### apply to CETS test and external data for benchmark ####
# library(GEOquery)
# untar(tarfile="GSE112179/GSE112179_RAW.tar", exdir="GSE112179/IDATs")
# setwd("/mnt/data1/Thea/ErrorMetric/data/externalData/GSE112179/IDATs")
# sapply(list.files(), gunzip)
# 
# GSE112179 = getGEO(filename = "../GSE112179_series_matrix.txt.gz")
# pheno = GSE112179@phenoData@data
# pheno$Basename = substring(as.character(pheno$supplementary_file), 68, 97)
# pheno = pheno[, c("age:ch1","cell.type:ch1","dist.dx:ch1",            
#                   "pmi:ch1","race:ch1","sampleID:ch1","Sex:ch1","tissue:ch1",             
#                   "tissuebank:ch1","tissuebank.id:ch1", "Basename")]
# colnames(pheno) = c("age","cell.type","dist.dx",            
#                     "pmi","race","sampleID","Sex","tissue",             
#                     "tissuebank","tissuebank.id", "Basename")
# save(pheno, file = "../phenoGSE112179.Rdata")
# 
# library(wateRmelon)
# RGset <- read.metharray.exp(base = ".", targets = pheno)
# save(RGset, file="/mnt/data1/Thea/ErrorMetric/data/externalData/GSE112179//GSE112179RGsetEPIC.Rdata")


load("/mnt/data1/Thea/ErrorMetric/data/externalData/GSE112179//GSE112179RGsetEPIC.Rdata")
library(minfi)
betas = getBeta(RGset)


## all samples are NeuN+
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

betasCETSTrain = betasCETSTrain[rownames(betasCETSTrain) %in% rownames(betas),]

CETSmodelGeo = pickCompProbes(rawbetas = betasCETSTrain,
                           cellTypes = levels(as.factor(as.character(phenoCETSTrain$Celltype))),
                           cellInd = as.factor(as.character(phenoCETSTrain$Celltype)),
                           numProbes = 50,
                           probeSelect = "auto")

pred = as.data.frame(projectCellTypeWithError(betas, modelType = "ownModel", ownModelData = CETSmodelGeo))
pred$dat = "PAI"

x = as.data.frame(projectCellTypeWithError(betasCETSTest, modelType = "ownModel", ownModelData = CETSmodelGeo))
x$dat = "CETS Test"

library(ggplot2)
library(cowplot)

plotDat = rbind.data.frame(pred, x)
# plotDat = rbind.data.frame(pred[,-which(colnames(pred) =="dist")], x)

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSERRORviolinPaiCets.pdf")
ggplot(plotDat, aes(x = dat, y = error, fill = dat)) +
  geom_violin() +
  theme_cowplot(18) +
  scale_fill_manual(values = c("#c2c2c2", "white")) +
  theme(legend.position = "none") +
  labs(x = "Dataset", y = "DSRMSE")
dev.off()

x = t.test(plotDat$error[plotDat$dat == "PAI"],plotDat$error[plotDat$dat != "PAI"])

pred$dist = abs(pred$`NeuN+`-1)
pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/ERRORCETSabsdifVSerror.pdf")
ggplot(pred, aes(x = dist, y = error))+
  geom_point() + 
  theme_cowplot(18) +
  labs(x = "Absolute difference between true and\npredicted proportion of NeuN+ ", y = "DSRMSE")
dev.off()

cor(pred$dist, pred$error)




### Create model in Caucasians and test in Caucasian and African ####
## to see if ethnic differences are picked up
load("/mnt/data1/Thea/humanDeconvolution/data/CETSUnnormalised.RData")

## n needed in test
length(levels(phenoCETS$Individual))*.5 # 14

## subset to NeuN+ to have individuals not samples
temp = phenoCETS[phenoCETS$Celltype == "NeuN+",]

## randomly select 12 Caucasian samples
temp = temp[temp$Ethnicity == "Caucasian",]
set.seed(1234)
CTest = temp$Individual[sample(nrow(temp), 12)]

phenoCETS$TrainTestnNeeded = "Test"
phenoCETS$TrainTestnNeeded[phenoCETS$Individual %in% as.numeric(as.character(CTest))] = "Train"

betasCETSTest = betasCETS[,phenoCETS$TrainTestnNeeded == "Test"]
phenoCETSTest = phenoCETS[phenoCETS$TrainTestnNeeded == "Test",]

betasCETSTrain = betasCETS[,phenoCETS$TrainTestnNeeded == "Train"]
phenoCETSTrain = phenoCETS[phenoCETS$TrainTestnNeeded == "Train",]

save(betasCETSTest, phenoCETSTest,betasCETSTrain, phenoCETSTrain,
     file = "/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTestEthnicity.RData")

## create model
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTestEthnicity.RData")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")

CETSmodel = pickCompProbes(rawbetas = betasCETSTrain,
                           cellTypes = levels(as.factor(as.character(phenoCETSTrain$Celltype))),
                           cellInd = as.factor(as.character(phenoCETSTrain$Celltype)),
                           numProbes = 50,
                           probeSelect = "auto")

save(CETSmodel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpGEthnicity.Rdata")

## apply model
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpGEthnicity.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTestEthnicity.RData")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")

pred = as.data.frame(projectCellTypeWithError(betasCETSTest, modelType = "ownModel", ownModelData = CETSmodel))

plotDat = cbind.data.frame(pred, phenoCETSTest)

library(ggplot2)
library(cowplot)

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSEthnicityViolin.pdf", height = 6, width = 6) 
ggplot(plotDat, aes(x = Ethnicity, y = error, fill = Ethnicity)) +
  geom_violin() +
  theme_cowplot(18) +
  scale_fill_manual(values = c("#c2c2c2", "white")) +
  theme(legend.position = "none") +
  labs(y = "Cetygo", x = "Ethnicity") +
  ylim(c(0,0.1))
dev.off()

t.test(plotDat$error[which(plotDat$Ethnicity == "Caucasian")], plotDat$error[which(plotDat$Ethnicity != "Caucasian")])
