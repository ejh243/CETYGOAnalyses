## train model for whole blood for use with CETYGO functions


cellTypes.450k = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")

library(minfi)
library(FlowSorted.Blood.450k)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(CETYGO)

### select probes for decovolution of whole blood using Reinus 6 purified cell types

### Load 450K data #########
referencePkg <- "FlowSorted.Blood.450k"
data(list = referencePkg)
referenceRGset <- get(referencePkg)
phenoDat.450k <- pData(referenceRGset)

## only keep 6 commonly used cell types
index <- phenoDat.450k$CellType %in% cellTypes.450k
phenoDat.450k <- phenoDat.450k[index,] 
betas.450k <- referenceRGset[,index]
probeAnno.450k <- getAnnotation(betas.450k)
autoProbes.450k <- rownames(probeAnno.450k)[probeAnno.450k$chr %in% paste0("chr", 1:22)]

## filter to autosomal sites
betas.450k  = subsetByLoci(betas.450k , includeLoci = autoProbes.450k)

## filter to overlapping epic
library(IlluminaHumanMethylationEPICmanifest)
data(IlluminaHumanMethylationEPICmanifest)
epicProbes<-c(getProbeInfo(IlluminaHumanMethylationEPICmanifest, type = c("I"))$Name, getProbeInfo(IlluminaHumanMethylationEPICmanifest, type = c("II"))$Name)
betas.450k  = subsetByLoci(betas.450k , includeLoci = epicProbes)

# convert to matrix
betas.450k = getBeta(preprocessRaw(betas.450k))

model.blood <- pickCompProbesMatrix(rawbetas = betas.450k,
									   cellTypes = cellTypes.450k,
									   cellInd = phenoDat.450k$CellType,
									   numProbes = 100,
									   probeSelect = "auto")
									   
save(model.blood, file = "Models/Blood.Reinus.6CT.rdata")									   