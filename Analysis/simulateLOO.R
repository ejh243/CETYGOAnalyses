## test a model trained on less cell types than included in bulk tissue generation

set.seed(9278)
cellTypes.450k = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")


## mean proportions of each cell type in whole blood (from Reinius2012)
meanBloodProp.450k = c(3.01,13.4, 6.13, 68.8, 5.4, 2.43)
names(meanBloodProp.450k)<-cellTypes.450k

## proportion of missing cell type in bulk profile
propMissing<-seq(0.1,0.9,0.1)

library(minfi)
library(FlowSorted.Blood.450k)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
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

predPropLOO<-NULL
for(ind in unique(phenoDat.450k$SampleID)){
	## select individual to keep as test data
	testIndex <- which(phenoDat.450k$SampleID == ind)
	## ensure test beta columns in same order as cell types
	testIndex<-testIndex[order(phenoDat.450k$CellType[testIndex])] 

	## create 5 cell type model dropping a each cell type in turn
	for(i in 1:length(cellTypes.450k)){
		
		## generate profiles of cellular proportions
		matrixSimProp<-generateCellProportions(meanBloodProp.450k[-i], propMissing)
		## relabel column
		colnames(matrixSimProp)[length(cellTypes.450k)]<-cellTypes.450k[i]
		## reorganise
		matrixSimProp<-matrixSimProp[,cellTypes.450k]
	
		model <- pickCompProbesMatrix(rawbetas = betas.450k[,-testIndex],
									   cellTypes = cellTypes.450k[-i],
									   cellInd = phenoDat.450k$CellType[-testIndex],
									   numProbes = 50,
									   probeSelect = "auto")

		## create test data
		testBulkBetas <-createBulkProfiles(betas.450k[rownames(model$coef),testIndex], matrixSimProp)
		## do cellular prediction with error
		predProp<-projectCellTypeWithError(testBulkBetas, model$coef)
		
		predProp<-cbind(predProp[,"CETYGO"], propMissing, cellTypes.450k[i] ) 	
		## save results
		predPropLOO<-rbind(predPropLOO, predProp)
	}
}

save(predPropLOO, file = "RData/ResultsLOOModels.rdata")


