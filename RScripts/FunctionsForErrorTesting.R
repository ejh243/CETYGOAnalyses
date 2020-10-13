### Functions for error testing

### CellTypeSubsetBetasAndPheno #######################
### subsetting betas and phenotype using cell type information given the cell type names wanted to keep 

##  INPUT: array of wanted cell type names
##         betas
##         pheno
##         cell type column name
##         bothOrJustBetas - can return just betas and ignore pheno
## The order of betas and pheno must be matching

## OUTPUT: betas and pheno in a list 

# ## test data
# load("/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
# betas = betas[1:100,]
# cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes")
# cellTypeNames = c("B_cells")
# # cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes", "testFake")
# phenoColName = "Sample.Type"
# 
# x = CellTypeSubsetBetasAndPheno(cellTypeNames, betas, pheno, phenoColName, bothOrJustBetas = T)

CellTypeSubsetBetasAndPheno = function(cellTypeNames, betas, pheno, phenoColName, bothOrJustBetas = F){
  
  ## ensure that phenocol is a factor
  phenoCol = as.factor(pheno[,phenoColName])
  
  ## check that cell types wanted are present
  if(all(cellTypeNames %in% phenoCol) == F){
    message(paste("At least one cell type given is not in the cell type column:", cellTypeNames[cellTypeNames %in% phenoCol == F], sep = " "))
  }
  
  betasOut = betas[, phenoCol %in% cellTypeNames]
  phenoOut = pheno[phenoCol %in% cellTypeNames, ]
  
  if (bothOrJustBetas == T){
    return(betasOut)
  }
  
  return(list(betasOut, phenoOut))
}



### randomNSumToProp ##################################
## calculates a row of random proportions that sum to totalProportion/100
## if random is included totalProportion should be 100-proportionNoise 
randomNSumToProp = function(totalProportion = 100, n = nCellTypes){
  x = sample(totalProportion+n-1, n-1)
  x = sort(x)
  prop = c(x[1] - 1)
  for (i in 2:(n-1)){
    prop[i] = x[i]-x[i-1]-1
  }
  prop[n] = totalProportion + n - 1 - x[n-1]
  prop = prop/100
  
  return(prop)
}



### CellTypeProportionSimulator #######################
### Simulate proportions of each cell type 

##  INPUT: betas       
##         pheno
##         phenoColName
##         proportionsMatrixType = "random" - can also be a matrix (nCellType * nBulk) 
##                                        colnames MUST match cell types
##         proportionsMatrix = NA - must have input if not "random"
##         nBulk - number of bulk samples wanted
##         noiseIn = F
##         proportionNoise = vector of noise same length as nBulk MUST HAVE VALUE IF noiseIn != F


## OUTPUT: matrix of bulk
##         matrix of proportions

## test data
load("/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
betas = betas[1:100,]
cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes")
cellTypeNames = c("B_cells")
# cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes", "testFake")
phenoColName = "Sample.Type"
nBulk = 12
noiseIn = F
proportionNoise = c(0.1,0.2,0.3,0.2,0.1,0.3)
proportionsMatrixType = "random"

CellTypeProportionSimulator(betas, pheno, phenoColName, nBulk, proportionsMatrix = "random")

CellTypeProportionSimulator = function(betas, pheno, phenoColName, nBulk, 
                                       proportionsMatrixType = "random", proportionsMatrix = NA, 
                                       noiseIn = F, proportionNoise = NA){
  
  ## ensure that phenocol is a factor
  phenoCol = as.factor(pheno[,phenoColName])
  
  ## get cell types
  cellTypes = levels(phenoCol)
  nCellTypes = length(cellTypes)
  
  ## separate each cell type out
  cellSamplesAlone = list()
  for(i in 1:nCellTypes){
    cellSamplesAlone[[i]] = CellTypeSubsetBetasAndPheno(cellTypes[i], betas, pheno, phenoColName, bothOrJustBetas = T)
  }
  
  ## If not input, create proportionsMatrix
  if (proportionsMatrixType == "random"){
    if (noiseIn == F){
      proportionsMatrix = t(replicate(nBulk, randomNSumToProp(totalProportion = 100)))
      colnames(proportionsMatrix) = cellTypes
    } else{
      proportionsMatrix = cbind(t(sapply(100 - proportionNoise*100, randomNSumToProp)), proportionNoise)
      colnames(proportionsMatrix) = c(cellTypes, "Noise")
    }
  }
  
  ## randomly pick one from each and sum with proportionsMatrix 
  CalculateBulk = function(proportionRow,
                           noiseInCB = noiseIn, 
                           nCellTypesCB = nCellTypes, 
                           cellSamplesAloneCB = cellSamplesAlone,
                           proportionsMatrixCB = proportionsMatrix){
    if (noiseInCB == F){
      simBulk = matrix(ncol = nCellTypesCB, nrow = nrow(cellSamplesAloneCB[[1]]))
      for (i in 1:nCellTypesCB){
        simBulk[,i] = cellSamplesAloneCB[[i]][,sample(ncol(cellSamplesAloneCB[[i]]), 1)] * proportionsMatrixCB[proportionRow,i]
      }
      simBulk = rowSums(simBulk)
      return(simBulk)
    }else{
      simBulk = matrix(ncol = nCellTypesCB+1, nrow = nrow(cellSamplesAloneCB[[i]]))
      for (i in 1:nCellTypesCB){
        simBulk[,i] = cellSamplesAloneCB[[i]][,sample(ncol(cellSamplesAloneCB[[i]]), 1)] * proportionsMatrixCB[proportionRow,i]
      }
      noiseSample = sample(nCellTypesCB, 1)
      simBulk[,nCellTypesCB+1] = cellSamplesAloneCB[[noiseSample]][,sample(ncol(cellSamplesAloneCB[[noiseSample]]), 1)] 
      
      simBulk[,nCellTypesCB+1] = simBulk[sample(nrow(simBulk), nrow(simBulk)), nCellTypesCB+1] * proportionsMatrixCB[proportionRow,nCellTypesCB+1]
      simBulk = rowSums(simBulk)
      return(simBulk)
    }
  }
  
  bulk = sapply(1:nBulk, CalculateBulk)
  colnames(bulk) = rownames(proportionsMatrix) = paste("S", 1:ncol(bulk), sep = "")
  
  return(list(bulk, proportionsMatrix))
}




### GetModelCG ########################################
### Subset the CpGs in betas for those used in a model
GetModelCG = function(betas, model){
  return(betas[rownames(betas) %in% rownames(model$coefEsts),])
}



### PredictionErrorAndResiduals #######################
### From model, predict proportions with error and output error vs residuals of true vs actual in testing

##  INPUT: model
##         testBetas
##         testPheno
##         phenoColName
##         

## OUTPUT: plot x = abs(true- pred), y = error, col = cell type

## test data
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForBrainCellProportionPrediction.r")
library(quadprog)
library(tidyr)
library(ggplot2)
library(cowplot)

# model = pickCompProbes(rawbetas = as.matrix(betas),
#                        cellTypes = levels(as.factor(pheno$Sample.Type)),
#                        cellInd = as.factor(pheno$Sample.Type),
#                        numProbes =  10,
#                        probeSelect = "auto")
# 
# test = CellTypeProportionSimulator(GetModelCG(betas, model), 
#                                         pheno, 
#                                         phenoColName = "Sample.Type", 
#                                         nBulk = 15, 
#                                         proportionsMatrixType = "random")
# testBetas = test[[1]]
# testPheno = test[[2]]

PredictionErrorAndResiduals = function(model, testBetas, testPheno){
  
  predictedPheno = projectCellTypeWithError(YIN = testBetas, coefCellTypeIN = model$coefEsts, sampleDup = 0)
  
  x = gather(cbind.data.frame(predictedPheno, sample = rownames(predictedPheno)), key = "cellType", value = "proportion_pred", 
             -one_of(c("error","nCGmissing", "sample")))
  y = gather(cbind.data.frame(testPheno, rownames(testPheno)), key = "cellType", value = "proportion_true", 
             -one_of(c("rownames(testPheno)")))

  plotDat = cbind.data.frame(x[,c("cellType", "proportion_pred","error","sample")] ,y["proportion_true"])  
  plotDat$residuals = abs(plotDat$proportion_pred - plotDat$proportion_true)
  
  print(ggplot(plotDat, aes(x = residuals, y = error, col = cellType, shape = cellType)) +
    geom_point(size = 2) +
    theme_cowplot(18) +
    ylim(c(0,max(plotDat$error))))
}










































