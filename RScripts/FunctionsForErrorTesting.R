### Functions for error testing

library(quadprog)
library(genefilter)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)

### CellTypeSubsetBetasAndPheno #######################
### subsetting betas and phenotype using cell type information given the cell type names wanted to keep 

##  INPUT: array of wanted cell type names
##         betas
##         pheno
##         cell type column name
##         justBetas - can return just betas and ignore pheno
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
# x = CellTypeSubsetBetasAndPheno(cellTypeNames, betas, pheno, phenoColName, justBetas = T)

CellTypeSubsetBetasAndPheno = function(cellTypeNames, betas, pheno, phenoColName, justBetas = F){
  
  ## ensure that phenocol is a factor
  phenoCol = as.factor(pheno[,phenoColName])
  
  ## check that cell types wanted are present
  if(all(cellTypeNames %in% phenoCol) == F){
    message(paste("At least one cell type given is not in the cell type column:", cellTypeNames[cellTypeNames %in% phenoCol == F], sep = " "))
  }
  
  betasOut = betas[, phenoCol %in% cellTypeNames]
  phenoOut = pheno[phenoCol %in% cellTypeNames, ]
  
  if (justBetas == T){
    return(betasOut)
  }
  
  return(list(betasOut, phenoOut))
}



### randomNSumToProp ##################################
## calculates a row of random proportions that sum to totalProportion/100
## if random is included totalProportion should be 100-proportionNoise 
randomNSumToProp = function(totalProportion = 100, n ){
  
  if(n==2){
    x = sample(100,1)
    prop = c(x/100, (100-x)/100)
    return(prop)
  }
  
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

# ## test data
# load("/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
# betas = betas[1:100,]
# cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes")
# cellTypeNames = c("B_cells")
# # cellTypeNames = c("B_cells", "CD4.T_cells", "CD8.T_cells", "Granulocytes", "testFake")
# phenoColName = "Sample.Type"
# nBulk = 8
# noiseIn = F
# proportionNoise = c(0.1,0.2,0.3,0.2,0.1,0.3)
# proportionsMatrixType = "random"
# 
# CellTypeProportionSimulator(betas, pheno, phenoColName, nBulk, proportionsMatrix = "random")

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
    cellSamplesAlone[[i]] = CellTypeSubsetBetasAndPheno(cellTypes[i], betas, pheno, phenoColName, justBetas = T)
  }
  
  ## If not input, create proportionsMatrix
  if (proportionsMatrixType == "random"){
    if (noiseIn == F){
      proportionsMatrix = t(replicate(nBulk, randomNSumToProp(totalProportion = 100, nCellTypes)))
      colnames(proportionsMatrix) = cellTypes
    } else{
      proportionsMatrix = cbind(t(sapply(100 - proportionNoise*100, randomNSumToProp, nCellTypes)), proportionNoise)
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
  rownames(bulk) = rownames(betas)
  
  return(list(bulk, proportionsMatrix))
}




### GetModelCG ########################################
### Subset the CpGs in betas for those used in a model
## model must be a list
GetModelCG = function(betas, modelList){
  ## get unique CpGs of models
  rownamesAll = c()
  for (i in 1:length(modelList)){
    model = modelList[[i]]
    rownamesAll = c(rownamesAll, rownames(model[["coefEsts"]]))
  }
  
  rownamesAll = unique(rownamesAll)
  
  betasMismatched = betas[rownames(betas) %in% rownamesAll,]
  
  betas = betasMismatched[match(rownamesAll, rownames(betasMismatched)),]
  return(betas)
}



### PredictionErrorAndResiduals #######################
### From model, predict proportions with error and output error vs residuals of true vs actual in testing

##  INPUT: model
##         testBetas - only those in model
##         testPheno

## OUTPUT: plot x = abs(true- pred), y = error, col = cell type

# ## test data
# source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForBrainCellProportionPrediction.r")
# library(quadprog)
# library(tidyr)
# library(ggplot2)
# library(cowplot)
# 
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
          ylim(c(0,max(plotDat$error)))+
          labs(x = "|Actual - Predicted|", y = "Error", shape = "Cell type", col = "Cell type"))
}




### ModelCompareStackedBar ############################

##  INPUT: testBetas - can containing all betas present in any model
##         modelList - must be list form even if only 1 model
##                     Should be a named list do that names can be in plots
##         trueComparison = F - also plot true proportions?
##                              if T, the proportions must be in trueProportions
##         noise = F - If there's noise in trueProportions change to T
##                     will be coloured in grey not a cell type
##         trueProportions = NA Matrix of true proportions from CellTypeProportionSimulator
##         nCpGPlot = F - plot number of CpGs?
##         sampleNamesOnPlots = F

## OUTPUT: stackedBarCharts per model for nBulk samples

# ## test data
# source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForBrainCellProportionPrediction.r")
# load("/mnt/data1/Thea/ErrorMetric/data/bloodEriskDataTrain.Rdata")
# library(quadprog)
# library(genefilter)
# library(tidyr)
# library(ggplot2)
# library(cowplot)
# library(scales)
# 
# model1 = pickCompProbes(rawbetas = as.matrix(betas),
#                         cellTypes = levels(as.factor(pheno$Sample.Type)),
#                         cellInd = as.factor(pheno$Sample.Type),
#                         numProbes =  10,
#                         probeSelect = "auto")
# 
# model2 = pickCompProbes(rawbetas = as.matrix(betas),
#                         cellTypes = levels(as.factor(pheno$Sample.Type)),
#                         cellInd = as.factor(pheno$Sample.Type),
#                         numProbes =  20,
#                         probeSelect = "auto")
# 
# phenoNewCellType = pheno
# phenoNewCellType[pheno$Sample.Type == "CD4.T_cells", "Sample.Type"] = "Other_cell"
# 
# model3 = pickCompProbes(rawbetas = as.matrix(betas),
#                         cellTypes = levels(as.factor(phenoNewCellType$Sample.Type)),
#                         cellInd = as.factor(phenoNewCellType$Sample.Type),
#                         numProbes =  10,
#                         probeSelect = "auto")
# 
# test = CellTypeProportionSimulator(GetModelCG(betas, list(model1, model2, model3)),
#                                    pheno,
#                                    phenoColName = "Sample.Type",
#                                    nBulk = 8,
#                                    proportionsMatrixType = "random")
# nBulk = 8
# testNoise = CellTypeProportionSimulator(GetModelCG(betas, list(model1, model2, model3)),
#                                         pheno,
#                                         phenoColName = "Sample.Type",
#                                         nBulk = 8,
#                                         proportionsMatrixType = "random",
#                                         noiseIn = T,
#                                         proportionNoise = seq(0.01, 0.1, 0.01)[1:nBulk])
# 
# testBetas = testNoise[[1]]
# trueProportions = testNoise[[2]]
# modelList = list(NameTest1 = model1, SecondModel = model2, thirds = model3)
# i = 1
# modelList = list(NameTest1 = model1)
# 
# trueComparison = T
# noise = T
# nCpGPlot = F
# sampleNamesOnPlots = F
# 
# x = ModelCompareStackedBar(testBetas,
#                            list(model1), nCpGPlot = F, 
#                            sampleNamesOnPlots = F)
# x = ModelCompareStackedBar(testBetas,
#                            modelList,
#                            trueComparison = T,
#                            noise = T,
#                            trueProportions = trueProportions)

ModelCompareStackedBar = function(testBetas, 
                                  modelList, 
                                  trueComparison = F,
                                  noise = F,
                                  trueProportions = NA,
                                  nCpGPlot = F,
                                  sampleNamesOnPlots = F){
  
  plotList = list()
  
  if(is.null(names(modelList))){
    message("Name your models!!")
  }
  
  ## predict cell proportions with each model and save in a data frame
  predictions = data.frame(matrix(ncol = 6, nrow = 0))
  for(i in 1:length(modelList)){
    model = modelList[[i]]
    newPred = projectCellTypeWithError(YIN = GetModelCG(testBetas, list(modelList[[i]])), 
                                       coefCellTypeIN = model[["coefEsts"]])
    newPred = cbind.data.frame(newPred, sample = rownames(newPred), model = names(modelList)[i])
    newPred = gather(newPred, key = "cellType", value = "proportion_pred", 
                     -one_of(c("error","nCGmissing", "sample", "model")))
    
    predictions = rbind.data.frame(predictions, newPred)
  }
  
  ## plot error per sample per model
  if(sampleNamesOnPlots){
    plotList[[1]] = ggplot(predictions, aes(x = sample, y = error, col = model)) +
      geom_point() +
      theme_cowplot(18) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(x = "Sample", y = "Error") +
      ylim(c(0, max(predictions$error)))
  }else{plotList[[1]] = ggplot(predictions, aes(x = sample, y = error, col = model)) +
    geom_point() +
    theme_cowplot(18) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(x = "Sample", y = "Error") +
    ylim(c(0, max(predictions$error)))
  }
  
  if(nCpGPlot == T){
    ## plot nCpG per sample per model
    if(sampleNamesOnPlots){
      plotList[[2]] = ggplot(predictions, aes(x = sample, y = nCGmissing, col = model)) +
        geom_point() +
        theme_cowplot(18) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(x = "Sample", y = "CGs missing") +
        ylim(c(0, max(predictions$error)))
    }else{
      plotList[[2]] = ggplot(predictions, aes(x = sample, y = nCGmissing, col = model)) +
        geom_point() +
        theme_cowplot(18) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        labs(x = "Sample", y = "CGs missing") +
        ylim(c(0, max(predictions$error)))
    }
  }
  
  if(trueComparison == T){
    predictions = predictions[,-which(colnames(predictions) %in% c("error", "nCGmissing"))]
    
    trueProportions = cbind.data.frame(trueProportions, sample = rownames(trueProportions), model = "True proportions")
    trueProportions = gather(trueProportions, key = "cellType", value = "proportion_pred", 
                             -one_of(c("sample", "model")))
    predictions = rbind.data.frame(predictions, trueProportions)
    
  }
  
  if(noise == F){
    predictions$cellType = as.factor(predictions$cellType)
    allCellTypes = levels(predictions$cellType)
    
    ## create one colour for each cell type
    colPerCell = hue_pal()(length(allCellTypes))
    
    
    ## plot stacked  bar chart
    for(i in 1:length(levels(as.factor(predictions$model)))){
      plotDat = predictions[predictions$model == levels(as.factor(predictions$model))[i],]
      
      colNeeded = which(allCellTypes %in% levels(as.factor(as.character(plotDat$cellType))))
      
      if(sampleNamesOnPlots){
        plotList[[length(plotList)+1]] = ggplot(plotDat, aes(x = sample, y = proportion_pred, fill = cellType)) +
          geom_bar(position="stack", stat="identity") +
          theme_cowplot(18) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(x = "Sample", y = "Proportion", fill = "Cell type") +
          ggtitle(levels(as.factor(predictions$model))[i]) +
          scale_fill_manual(values = c(colPerCell[colNeeded]))
      }else{
        plotList[[length(plotList)+1]] = ggplot(plotDat, aes(x = sample, y = proportion_pred, fill = cellType)) +
          geom_bar(position="stack", stat="identity") +
          theme_cowplot(18) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          labs(x = "Sample", y = "Proportion", fill = "Cell type") +
          ggtitle(levels(as.factor(predictions$model))[i]) +
          scale_fill_manual(values = c(colPerCell[colNeeded]))
      }
    }
  }else{
    predictions$cellType = as.factor(predictions$cellType)
    allCellTypes = levels(predictions$cellType)
    
    ## create one colour for each cell type with grey for noise
    hues = hue_pal()(length(allCellTypes)-1)
    colPerCell = c(matrix(nrow = 1, ncol = length(allCellTypes)))
    colPerCell[which(allCellTypes == "Noise")] = "#555353"
    colPerCell[-which(allCellTypes == "Noise")] = hues
    
    for(i in 1:length(levels(as.factor(predictions$model)))){
      plotDat = predictions[predictions$model == levels(as.factor(predictions$model))[i],]
      
      colNeeded = which(allCellTypes %in% levels(as.factor(as.character(plotDat$cellType))))
      
      if(sampleNamesOnPlots){
        plotList[[length(plotList)+1]] = ggplot(plotDat, aes(x = sample, y = proportion_pred, fill = cellType)) +
          geom_bar(position="stack", stat="identity") +
          theme_cowplot(18) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          labs(x = "Sample", y = "Proportion", fill = "Cell type") +
          ggtitle(levels(as.factor(predictions$model))[i]) +
          scale_fill_manual(values = c(colPerCell[colNeeded]))
      }else{
        plotList[[length(plotList)+1]] = ggplot(plotDat, aes(x = sample, y = proportion_pred, fill = cellType)) +
          geom_bar(position="stack", stat="identity") +
          theme_cowplot(18) +
          theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          labs(x = "Sample", y = "Proportion", fill = "Cell type") +
          ggtitle(levels(as.factor(predictions$model))[i]) +
          scale_fill_manual(values = c(colPerCell[colNeeded]))
      }
    }
  }
  return(plotList)
}



### TrueVSPredictedPlot ###############################


##  INPUT: pred - predicted cell types
##         true - true proportions

## OUTPUT: plot true vs pred

# ## test data
# pred 
# true = bulkPheno


TrueVSPredictedPlot = function(pred, true){
  
  plotDat = cbind.data.frame(gather(data.frame(pred), key = "cellType", value = "proportion_pred", 
                                    -one_of(c("error","nCGmissing"))),
                             gather(data.frame(true), key = "cellType_true", value = "proportion_true"))
  
  
  print(ggplot(plotDat, aes(x = proportion_pred, y = proportion_true, col = cellType, shape = cellType)) +
    geom_point() +
    theme_cowplot(18) +
    labs(x= "Predicted", y = "Actual", shape = "Cell type", col = "Cell type"))
  
}


### TrueVSPredictedMaintainCellsPlot ###############################


##  INPUT: pred - predicted cell types
##         true - true proportions
##         cells = c(T,T,T,T,T) vector of 5 booleans to show which cell type in factor are present

## OUTPUT: plot true vs pred

# ## test data
# pred = CD4Model
# true = C4C8bulk
# cells = c(T,T,F,T,T)


TrueVSPredictedMaintainCellsPlot = function(pred, true, cells = c(T,T,T,T,T), ylabel = "Actual"){
  
  shapeOptions = c(16, 17, 15, 3, 7)
  shapes = shapeOptions[cells]
  
  colourOptions = hue_pal()(5)
  colours = colourOptions[cells]
  
  plotDat = cbind.data.frame(gather(data.frame(pred), key = "cellType", value = "proportion_pred", 
                                    -one_of(c("error","nCGmissing"))),
                             gather(data.frame(true), key = "cellType_true", value = "proportion_true"))
  
  
  print(ggplot(plotDat, aes(x = proportion_pred, y = proportion_true, col = cellType, shape = cellType)) +
          geom_point() +
          theme_cowplot(18) +
          labs(x= "Predicted", y = ylabel, shape = "Cell type", col = "Cell type") +
          scale_shape_manual(values=shapes) +
          scale_color_manual(values=colours))
  
}



### ErrorAcrossDataSets ###############################

##  INPUT: list of data sets with nCG from model
##         model to be applied

## OUTPUT: error per dataset in one data frame

# ## test data
# dataList = x
# i = 2

ErrorAcrossDataSets = function(dataList, model){
  
  outDat = matrix(ncol = 7, nrow = 0)
  for (i in 1:length(dataList)){
    outDat = rbind(outDat, projectCellTypeWithError(YIN = GetModelCG(dataList[[i]], list(model)), 
                             coefCellTypeIN = model[["coefEsts"]]))
  }
  return(as.data.frame(outDat))
}