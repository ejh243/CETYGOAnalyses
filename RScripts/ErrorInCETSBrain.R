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

load("/mnt/data1/Thea/humanDeconvolution/data/CETSUnnormalised.RData")
## subset to NeuN+ to have individuals not samples
temp = phenoCETS[phenoCETS$Celltype == "NeuN+",]
tempA = temp[temp$Ethnicity == "African",]

## randomly select 2 African samples
set.seed(1234)
ATest = tempA$Individual[sample(nrow(tempA), 4)]

## randomly select 7 Caucasian samples
temp = temp[temp$Ethnicity == "Caucasian",]
set.seed(1234)
CTest = temp$Individual[sample(nrow(temp), 16)]

phenoCETS$TrainTestnNeeded = "Train"
phenoCETS$TrainTestnNeeded[phenoCETS$Individual %in% c(as.numeric(as.character(ATest)), as.numeric(as.character(CTest)))] = "Test"

betasCETSTest = betasCETS[,phenoCETS$TrainTestnNeeded == "Test"]
phenoCETSTest = phenoCETS[phenoCETS$TrainTestnNeeded == "Test",]

betasCETSTrain = betasCETS[,phenoCETS$TrainTestnNeeded == "Train"]
phenoCETSTrain = phenoCETS[phenoCETS$TrainTestnNeeded == "Train",]

save(betasCETSTest, phenoCETSTest,betasCETSTrain,phenoCETSTrain,
     file = "/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")


### create model using test data ######################
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/pickCompProbes.R")

CETSmodel = pickCompProbes(rawbetas = betasCETSTrain,
                           cellTypes = levels(as.factor(as.character(phenoCETSTrain$Celltype))),
                           cellInd = as.factor(as.character(phenoCETSTrain$Celltype)),
                           numProbes = 50,
                           probeSelect = "auto")

save(CETSmodel, file = "/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")

### compare error with noise ##########################
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/CETScols.Rdata")
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

stackedWithNoise[[1]]
stackedWithNoise[[2]] + scale_fill_manual(values = c("#714a91", "#4ab9b0"))
stackedWithNoise[[3]] + scale_fill_manual(values = c("#714a91", "#4ab9b0", "#555353"))

leg = get_legend(stackedWithNoise[[3]] + scale_fill_manual(values = c("#714a91", "#4ab9b0", "#555353"))
                 + theme(legend.position=c(0.3,0.8),legend.direction = "horizontal", legend.title = element_blank())
                 + guides(fill = guide_legend(nrow = 1)))
plots = plot_grid(stackedWithNoise[[1]] + theme(legend.position = "none"),
                  stackedWithNoise[[2]] + theme(legend.position = "none") +
                    scale_y_continuous(breaks = seq(0,1,0.25))+
                    scale_fill_manual(values = c("#714a91", "#4ab9b0")),
                  stackedWithNoise[[3]] + theme(legend.position = "none") + scale_fill_manual(values = c("#714a91", "#4ab9b0", "#555353")), ncol = 1,
                  rel_heights = c(0.6,1,1), labels = "AUTO", axis = "rl", align = "v" )

pdf("/mnt/data1/Thea/ErrorMetric/plots/CETSValidation/CETSsimWithNoise.pdf", height = 9, width = 7.5)     
print(plot_grid(plots, leg, ncol = 1, rel_heights = c(1,0.08)))
dev.off()


### Apply to read data
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/CETScols.Rdata")
source("/mnt/data1/Thea/ErrorMetric/DSRMSE/projectCellTypeWithError.R")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

## load subset Essex data


