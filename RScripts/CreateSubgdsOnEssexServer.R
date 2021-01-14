## open gds file
library(gdsfmt)
x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)

betas = read.gdsn(index.gdsn(x, "rawbetas"))
#dim(betas)

load("~/DSRMSE/models/HousemanBloodModel150CpG.Rdata")

pID = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
betaIndex = pID %in% rownames(HousemanBlood150CpGModel$coefEsts)

rownames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
colnames(betas) = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "FullBarCode"))

betaMod = betas[betaIndex,]
betaMod = betaMod[match(rownames(HousemanBlood150CpGModel$coefEsts), rownames(betaMod)),]

## create column index, removing samples that have "", or unsorted
tDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Tissue"))
aDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Age"))
sDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Sex"))
dDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "DatasetOrigin"))
stDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "SubTissue"))
tInd = which(tDat == "" | 
               tDat == "Unsorted Tissues" | 
               tDat == "Unsorted Cell Line" | 
               tDat == "Unsorted Tumours")

betaMod = betaMod[,-tInd]

gfile = createfn.gds("sub.gds")

add.gdsn(gfile, "Age", aDat[-tInd])
add.gdsn(gfile, "Sex", sDat[-tInd])
add.gdsn(gfile, "Tissue", tDat[-tInd])
add.gdsn(gfile, "DatasetOrigin", dDat[-tInd])
add.gdsn(gfile, "SubTissue", stDat[-tInd])
add.gdsn(gfile, "rownames", rownames(betaMod))
add.gdsn(gfile, "colnames", colnames(betaMod))

#source("~/DSRMSE/FunctionsForErrorTesting.R")
source("~/DSRMSE/projectCellTypeWithError.R")
model = HousemanBlood150CpGModel

ind = c(seq(1000, 20960, 1000))
errPred = lapply(ind, function(ind){
  pred = projectCellTypeWithError(betaMod[,(ind - 999):ind], model = "ownModel", ownModelData = model)
  return(pred)})

errPred[[length(errPred)+1]] = projectCellTypeWithError(betaMod[,20001:20960], model = "ownModel", ownModelData = model)

## extract all from list and add to gds file
pred = errPred[[1]]
for (i in 2:length(errPred)){
  pred = rbind(pred, errPred[[i]])
}

add.gdsn(gfile, "Pred", pred)
closefn.gds(gfile)
closefn.gds(x)
q()

# scp sub.gds dSeiler@knight.ex.ac.uk:/mnt/data1/Thea/ErrorMetric/data/EssexOutput/
  

### Do a PCA of GSE89251 to see if the weird samples stand out
# GSE89251
# Low error bad prediction
# GSM2363099_9904796223_R06C01 GSM2363100_9904796223_R06C02
# 
# High error
# GSM2363086_9904796203_R05C01 GSM2363092_9904796223_R02C02 GSM2363095_9904796223_R04C01

library(gdsfmt)
x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)

betas = read.gdsn(index.gdsn(x, "rawbetas"))

dDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "DatasetOrigin"))

dat = betas[, which(dDat == "GSE89251.gds")]

## double check that they're all T cells
tDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Tissue"))
unique(tDat[which(dDat == "GSE89251.gds")]) # Yep!


library(ggfortify)
library(ggplot2)
library(cowplot)
library(scales)

## make weirdness a phenotype 
pheno = rep(0, ncol(dat))
pDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "FullBarCode"))[which(dDat == "GSE89251.gds")]
pheno[pDat == "GSM2363099_9904796223_R06C01" | 
        pDat == "GSM2363100_9904796223_R06C02"] = "Low true, low error"
pheno[pDat == "GSM2363086_9904796203_R05C01" |
        pDat == "GSM2363092_9904796223_R02C02" |
        pDat == "GSM2363095_9904796223_R04C01"] = "High error"
pheno = as.factor(pheno)

dat = dat[complete.cases(dat),]
betaVar = apply(dat, 1, var, na.rm = T)
topBeta = dat[order(betaVar, decreasing = T)[1:1000],]
plot = cbind.data.frame(prcomp(t(topBeta))$x[,1:2], pheno)

ggplot(plot,  aes(x = PC1, y = PC2, col = pheno)) +
  geom_point() +
  theme_cowplot(18)



### Compare my prediction with minfi for blood cell types
## load the data
library(gdsfmt)
library(minfi)

x = openfn.gds("/storage/st05d/deepmelon/GEOClod.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)
pheno = cbind.data.frame(dDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "DatasetOrigin")),
                         tDat = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "Tissue")))

cellIndex = pheno$tDat == "B Cells" |
  pheno$tDat == "Granulocyes" |
  pheno$tDat == "NK" |
  pheno$tDat == "T Cells"

pheno = pheno[cellIndex,]
pheno$tDat = as.factor(as.character(pheno$tDat))
pheno$dDat = as.factor(as.character(pheno$dDat))

betasIN = read.gdsn(index.gdsn(x, "rawbetas"))[,cellIndex]
rownames(betasIN) = read.gdsn(index.gdsn(index.gdsn(x$root, "fData"), "Probe_ID"))
colnames(betasIN) = read.gdsn(index.gdsn(index.gdsn(x$root, "pData"), "FullBarCode"))[cellIndex]

## load the reference data
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
betasREF = referenceRGset[,index]
betasREF = getBeta(preprocessRaw(betasREF))

## match rows in the data
betasIN = betasIN[which(rownames(betasIN) %in% rownames(betasREF)),]
betasIN = betasIN[match(rownames(betasREF), rownames(betasIN)),]
all(rownames(betasIN) == rownames(betasREF))


datasetIndex = 1

source("~/DSRMSE/pickCompProbes.R")
source("~/DSRMSE/projectCellTypeWithError.R")



normNPred = function(datasetIndex, betasREF, betasIN, pheno, phenoDat){
  ## per dataset, combine with minfi data, normalise using quantile and predict proportions
  
  ## subset for study
  betasStudy = betasIN[, pheno$dDat == levels(pheno$dDat)[datasetIndex]]
  
  ## normalise together
  betasNorm = normalizeQuantiles(cbind(betasREF, betasStudy))
  
  ## remove all rows with NAs
  betasNorm = betasNorm[which(apply(betasNorm, 1, function(x){sum(is.na(x))==0})),]
  
  ## create model
  model = pickCompProbes(rawbetas = betasNorm[,1:ncol(betasREF)],
                         cellTypes = levels(phenoDat),
                         cellInd = phenoDat,
                         numProbes =  150,
                         probeSelect = "auto")
 
  ## project cell proportions 
  return(projectCellTypeWithError(betasNorm[,37:ncol(betasNorm)], 
                                  model = "ownModel", 
                                  ownModelData = model)) 
}

pred = lapply(1:length(levels(pheno$dDat)), normNPred, betasREF, betasIN, pheno, phenoDat)

pheno$dDat = as.factor(unlist(strsplit(as.character(pheno$dDat), ".g"))[seq(1,nrow(pheno)*2, 2)])

## unlist pred
pdat = data.frame(pred[[1]][,1:6], DatasetOrigin = rep(levels(pheno$dDat)[1],nrow(pred[[1]])))
for (i in 2:length(pred)){
  pdat = rbind.data.frame(pdat, 
                          data.frame(pred[[i]][,1:6], 
                                     DatasetOrigin = rep(levels(pheno$dDat)[i],nrow(pred[[i]]))))
}

save(pdat, file = "minfiPred.Rdata")

## make sample a column
pdat$Sample = rownames(pdat)

library(reshape2)

minfiDat = melt(pdat, id.vars = c("Sample", "DatasetOrigin"), 
                meause.vars = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))
colnames(minfiDat)[3:4] = c("CellType", "minfi.Value")


## load my predictions
dsPred = openfn.gds("~/sub.gds", readonly=TRUE, allow.duplicate=FALSE, allow.fork=FALSE)

dsDat = data.frame(read.gdsn(index.gdsn(dsPred$root, "Pred"))[,1:7], 
                   Sample = read.gdsn(index.gdsn(dsPred$root, "colnames")))
colnames(dsDat) = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "error", "Sample")

dsDat = dsDat[dsDat$Sample %in% pdat$Sample,]
dsDat = dsDat[match(pdat$Sample, dsDat$Sample),]

errorDat = melt(dsDat, id.vars = c("Sample", "error"), 
                meause.vars = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"))

colnames(errorDat)[3:4] = c("CellType", "DSRMSE.Value")

## merge together
plotDat = merge(x = minfiDat, y = errorDat, by = c("Sample", "CellType"), all = T)

library(ggplot2)
library(cowplot)
library(viridis)

pdf("MinfiVSMyPredForValidation.pdf", height = 7, width = 8)
ggplot(plotDat, aes(x = DSRMSE.Value, y = minfi.Value, col = error, shape = CellType)) +
  geom_point() + 
  scale_color_viridis() +
  theme_cowplot(18) +
  labs(x = "Unnormalised model predictions", y = "Minfi predictions", 
       col = "DSRMSE", shape = "Cell type") +
  annotate("text", x = 0.1, y = 1, size = 5,
           label = paste("Cor =", signif(cor(plotDat$DSRMSE.Value, plotDat$minfi.Value), 3)))
dev.off()



