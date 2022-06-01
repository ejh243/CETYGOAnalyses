## test a model trained on less cell types than included in bulk tissue generation

set.seed(9278)
cellTypes.450k = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK")


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

## simulate full spectrun of blood profiles where each cell type is present at at least 10%
intervals<-seq(0.1,0.5,0.1) ## max proportion is 0.5 if all other cell types must be present at 0.1
matrixSimProp<-expand.grid(rep(list(intervals), length(cellTypes.450k)))
matrixSimProp<-matrixSimProp[which(rowSums(matrixSimProp) == 1),]

## make design matrix of booleans for which cell type will be present
designMatrix = expand.grid(c(T,F), c(T,F), c(T,F), c(T,F), c(T,F), c(T,F))
designMatrix = designMatrix[apply(designMatrix, 1, sum) >= 3 ,]
designMatrix<-as.matrix(designMatrix[-1,])


predPropCombined<-NULL
for(ind in unique(phenoDat.450k$SampleID)){
	## select individual to keep as test data
	testIndex <- which(phenoDat.450k$SampleID == ind)
	## ensure test beta columns in same order as cell types
	testIndex<-testIndex[order(phenoDat.450k$CellType[testIndex])] 

	## test all models with between 3-5 cell types
	for(i in 1:nrow(designMatrix)){
		model <- pickCompProbesMatrix(rawbetas = betas.450k[,-testIndex],
									   cellTypes = cellTypes.450k[designMatrix[i,]],
									   cellInd = phenoDat.450k$CellType[-testIndex],
									   numProbes = 50,
									   probeSelect = "auto")

		## create test data
		testBulkBetas <-createBulkProfiles(betas.450k[rownames(model$coef),testIndex], matrixSimProp)
		## do cellular prediction with error
		predProp<-projectCellTypeWithError(testBulkBetas, model$coef)
		
		nCT<-sum(designMatrix[i,])
		propRepresented<-rowSums(matrixSimProp[,designMatrix[i,]])
		predProp<-cbind(predProp[,"CETYGO"], nCT,propRepresented) 	
		## save results
		predPropCombined<-rbind(predPropCombined, predProp)
	}
}

save(predPropCombined, file = "RData/ResultsIncompleteModels.rdata")


