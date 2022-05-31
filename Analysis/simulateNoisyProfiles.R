## run simulations to test reconstructed whole blood profiles with added noise
## to avoid leakage don't normalise


set.seed(9278)
cellTypes.450k = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")
cellTypes.epic = c("Bcell", "CD4T", "CD8T", "Neu", "Mono", "NK")

## mean proportions of each cell type in whole blood (from Reinius2012)
meanBloodProp.450k = c(3.01,13.4, 6.13, 68.8, 5.4, 2.43)
names(meanBloodProp.450k)<-cellTypes.450k

meanBloodProp.epic = c(3.01,13.4, 6.13, 64.9, 5.4, 2.43)
names(meanBloodProp.epic)<-cellTypes.epic

nSims<-10 ## for epic where samples not grouped by individual
noise = c(seq(0,0.1,0.01), seq(0.12,0.5,0.02))

library(minfi)
library(FlowSorted.Blood.450k)
library(FlowSorted.Blood.EPIC)
library(ExperimentHub)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(ggplot2)
library(ggpubr)
library(CETYGO)

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

# convert to matrix
betas.450k = getBeta(preprocessRaw(betas.450k))


### Load EPIC data #########
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]

## only keep the 6 commonly used cell types NOTE only has Neu rather than Gran
index <- FlowSorted.Blood.EPIC$CellType %in% cellTypes.epic
phenoDat.epic <- pData(FlowSorted.Blood.EPIC)[index,] 
betas.epic <- FlowSorted.Blood.EPIC[,index]
probeAnno.epic <- getAnnotation(betas.epic)
autoProbes.epic <- rownames(probeAnno.epic)[probeAnno.epic$chr %in% paste0("chr", 1:22)]

## filter to autosomal sites
betas.epic  = subsetByLoci(betas.epic , includeLoci = autoProbes.epic)

# convert to matrix
betas.epic = getBeta(preprocessRaw(betas.epic))


## create matrix of proportions of cell types and noise for simulations
matrixBloodProp.450k<-generateCellProportions(meanBloodProp.450k, noise)
matrixBloodProp.epic<-generateCellProportions(meanBloodProp.epic, noise)

predPropAll.450k<-NULL
for(ind in unique(phenoDat.450k$SampleID)){
	## select individual to keep as test data
	testIndex <- which(phenoDat.450k$SampleID == ind)
	## ensure test beta columns in same order as cell types
	testIndex<-testIndex[order(phenoDat.450k$CellType[testIndex])] 

	predProp<-runNoiseSimulations(betas.450k[,-testIndex],cellTypes.450k,
		phenoDat.450k$CellType[-testIndex], betas.450k[,testIndex], matrixBloodProp.450k)
		
	## save results
	predPropAll.450k<-rbind(predPropAll.450k, cbind(noise, predProp))
}

# add some summary stats
predPropAll.450k<-as.data.frame(predPropAll.450k)

## calculate RMSE for prediction
predPropAll.450k<-cbind(predPropAll.450k, apply(predPropAll.450k[, cellTypes.450k],1, RMSE, meanBloodProp.450k/sum(meanBloodProp.450k)))

## calculate total prediction
predPropAll.450k<-cbind(predPropAll.450k, rowSums(predPropAll.450k[, cellTypes.450k]))
predPropAll.450k<-cbind(predPropAll.450k, "450k") 
colnames(predPropAll.450k)[10:12]<-c("RMSE", "TotalComposition", "Array")


### run for EPIC data
predPropAll.epic<-NULL
for(i in 1:nSims){
	## select samples to keep as test data
	testIndex<-NULL
	for(each in cellTypes.epic){
		testIndex<-c(testIndex, sample(which(phenoDat.epic$CellType == each), 1))
	}
	predProp<-runNoiseSimulations(betas.epic[,-testIndex],cellTypes.epic,
		phenoDat.epic$CellType[-testIndex], betas.epic[,testIndex], matrixBloodProp.epic)
			
	predPropAll.epic<-rbind(predPropAll.epic, cbind(noise, predProp))
}


# add some summary stats

## calculate RMSE for prediction
predPropAll.epic<-cbind(predPropAll.epic, apply(predPropAll.epic[, cellTypes.epic],1, RMSE, meanBloodProp.epic/sum(meanBloodProp.epic)))

## calculate total prediction
predPropAll.epic<-cbind(predPropAll.epic, rowSums(predPropAll.epic[, cellTypes.epic]))
predPropAll.epic<-as.data.frame(predPropAll.epic)
predPropAll.epic<-cbind(predPropAll.epic, "EPIC") 

colnames(predPropAll.epic)[10:12]<-c("RMSE", "TotalComposition", "Array")

## merge results into single dataframe for plotting
colnames(predPropAll.epic)<-colnames(predPropAll.450k)
predPropAll<-rbind(predPropAll.450k, predPropAll.epic)

save(predPropAll, file = "RData/ResultsNoisyProfilesWithinBatch")


## train in 450K test in EPIC with noise
## limit to shared cell types
cellTypes.share<-intersect(cellTypes.450k, cellTypes.epic)
matrixBloodProp.share<-generateCellProportions(meanBloodProp.epic[cellTypes.share], noise)

## limit to shared sites
sites.share<-intersect(rownames(betas.450k), rownames(betas.epic))
betas.450k<-betas.450k[sites.share,]
betas.epic<-betas.epic[sites.share,]

predPropAll.4te<-NULL
for(i in 1:nSims){
	## select samples to keep as test data
	testIndex<-NULL
	for(each in cellTypes.share){
		testIndex<-c(testIndex, sample(which(phenoDat.epic$CellType == each), 1))
	}
	predProp<-runNoiseSimulations(betas.450k,cellTypes.share,
		phenoDat.450k$CellType, betas.epic[,testIndex], matrixBloodProp.share)
			
	predPropAll.4te<-rbind(predPropAll.4te, cbind(noise, predProp))
}

## add some summary stats

## calculate RMSE for prediction
predPropAll.4te<-cbind(predPropAll.4te, apply(predPropAll.4te[, cellTypes.share],1, RMSE, meanBloodProp.epic[cellTypes.share]/sum(meanBloodProp.epic[cellTypes.share])))

## calculate total prediction
predPropAll.4te<-cbind(predPropAll.4te, rowSums(predPropAll.4te[, cellTypes.share]))
predPropAll.4te<-as.data.frame(predPropAll.4te)
predPropAll.4te<-cbind(predPropAll.4te, "450K : EPIC", "Unnormalised") 

colnames(predPropAll.4te)[9:12]<-c("RMSE", "TotalComposition", "Array", "Normalised")




predPropAll.et4<-NULL
## train in epic test in 450k$CellType
for(ind in unique(phenoDat.450k$SampleID)){
	## select individual to keep as test data
	testIndex <- which(phenoDat.450k$SampleID == ind)
	## ensure test beta columns in same order as cell types
	testIndex<-testIndex[match(cellTypes.share,phenoDat.450k$CellType[testIndex])] 
	
	predProp<-runNoiseSimulations(betas.epic,cellTypes.share,
		phenoDat.epic$CellType, betas.450k[,testIndex], matrixBloodProp.share)
			
	predPropAll.et4<-rbind(predPropAll.et4, cbind(noise, predProp))
}

## calculate RMSE for prediction
predPropAll.et4<-cbind(predPropAll.et4, apply(predPropAll.et4[, cellTypes.share],1, RMSE, meanBloodProp.epic[cellTypes.share]/sum(meanBloodProp.epic[cellTypes.share])))

## calculate total prediction
predPropAll.et4<-cbind(predPropAll.et4, rowSums(predPropAll.et4[, cellTypes.share]))
predPropAll.et4<-as.data.frame(predPropAll.et4)
predPropAll.et4<-cbind(predPropAll.et4, "EPIC : 450K", "Unnormalised") 

colnames(predPropAll.et4)[9:12]<-c("RMSE", "TotalComposition", "Array", "Normalised")


## normalise together and repeat
### Load 450K data #########
referencePkg <- "FlowSorted.Blood.450k"
data(list = referencePkg)
referenceRGset <- get(referencePkg)
phenoDat.450k <- pData(referenceRGset)

## only keep shared cell types
index <- phenoDat.450k$CellType %in% cellTypes.share
phenoDat.450k <- phenoDat.450k[index,] 
betas.450k <- referenceRGset[,index]


### Load EPIC data #########
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]

## only keep the 6 commonly used cell types NOTE only has Neu rather tha Gran
index <- FlowSorted.Blood.EPIC$CellType %in% cellTypes.share
phenoDat.epic <- pData(FlowSorted.Blood.EPIC)[index,] 
betas.epic <- FlowSorted.Blood.EPIC[,index]

betas.merge<-combineArrays(betas.450k,betas.epic)

quantileBetas = getBeta(preprocessQuantile(betas.merge, fixOutliers = TRUE,
                                            removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                            quantileNormalize = TRUE, stratified = TRUE, 
                                            mergeManifest = FALSE, sex = NULL))



## filter to autosomal sites
betas.merge <- quantileBetas[sites.share,]


## separate
betas.450k.norm<-betas.merge[,colnames(betas.merge) %in% colnames(betas.450k)]

betas.epic.norm<-betas.merge[,colnames(betas.merge) %in% colnames(betas.epic)]

predPropAll.4te.norm<-NULL
for(i in 1:nSims){
	## select samples to keep as test data
	testIndex<-NULL
	for(each in cellTypes.share){
		testIndex<-c(testIndex, sample(which(phenoDat.epic$CellType == each), 1))
	}
	predProp<-runNoiseSimulations(betas.450k.norm,cellTypes.share,
		phenoDat.450k$CellType, betas.epic.norm[,testIndex], matrixBloodProp.share)
			
	predPropAll.4te.norm<-rbind(predPropAll.4te.norm, cbind(noise, predProp))
}

## add some summary stats

## calculate RMSE for prediction
predPropAll.4te.norm<-cbind(predPropAll.4te.norm, apply(predPropAll.4te.norm[, cellTypes.share],1, RMSE, meanBloodProp.epic[cellTypes.share]/sum(meanBloodProp.epic[cellTypes.share])))

## calculate total prediction
predPropAll.4te.norm<-cbind(predPropAll.4te.norm, rowSums(predPropAll.4te.norm[, cellTypes.share]))
predPropAll.4te.norm<-as.data.frame(predPropAll.4te.norm)
predPropAll.4te.norm<-cbind(predPropAll.4te.norm, "450K : EPIC", "Normalised") 

colnames(predPropAll.4te.norm)[9:12]<-c("RMSE", "TotalComposition", "Array", "Normalised")




predPropAll.et4.norm<-NULL
## train in epic test in 450k$CellType
for(ind in unique(phenoDat.450k$SampleID)){
	## select individual to keep as test data
	testIndex <- which(phenoDat.450k$SampleID == ind)
	## ensure test beta columns in same order as cell types
	testIndex<-testIndex[match(cellTypes.share,phenoDat.450k$CellType[testIndex])] 
	
	predProp<-runNoiseSimulations(betas.epic.norm,cellTypes.share,
		phenoDat.epic$CellType, betas.450k.norm[,testIndex], matrixBloodProp.share)
			
	predPropAll.et4.norm<-rbind(predPropAll.et4.norm, cbind(noise, predProp))
}

## calculate RMSE for prediction
predPropAll.et4.norm<-cbind(predPropAll.et4.norm, apply(predPropAll.et4.norm[, cellTypes.share],1, RMSE, meanBloodProp.epic[cellTypes.share]/sum(meanBloodProp.epic[cellTypes.share])))

## calculate total prediction
predPropAll.et4.norm<-cbind(predPropAll.et4.norm, rowSums(predPropAll.et4.norm[, cellTypes.share]))
predPropAll.et4.norm<-as.data.frame(predPropAll.et4.norm)
predPropAll.et4.norm<-cbind(predPropAll.et4.norm, "EPIC : 450K", "Normalised") 

colnames(predPropAll.et4.norm)[9:12]<-c("RMSE", "TotalComposition", "Array", "Normalised")

predPropArrays<-rbind(predPropAll.et4, predPropAll.4te, predPropAll.et4.norm, predPropAll.4te.norm) 

save(predPropArrays, file = "RData/ResultsNoisyProfilesAcrossBatch.rdata")