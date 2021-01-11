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

