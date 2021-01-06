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

gfile = createfn.gds("sub.gds")

add.gdsn(gfile, "Age", aDat[-tInd])
add.gdsn(gfile, "Sex", sDat[-tInd])
add.gdsn(gfile, "Tissue", tDat[-tInd])
add.gdsn(gfile, "DatasetOrigin", dDat[-tInd])
add.gdsn(gfile, "SubTissue", stDat[-tInd])

betaMod = betaMod[,-tInd]
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

#scp sub.gds dSeiler@knight.ex.ac.uk:/mnt/data1/Thea/ErrorMetric/data/EssexOutput/
  