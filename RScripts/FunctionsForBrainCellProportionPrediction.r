## Main script by Eilis Hannon, adapted by Dorothea Seiler Vellame
## functions for estimate cell counts edited to take matrices rather than mset
# library(genefilter)
# library(quadprog)
# 
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]

  if(is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }

  ## Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()

  if(verbose)
    cat("[validationCellType] ")
  for(j in 1:M) { # For each CpG
    ## Remove missing methylation values
    ii <- !is.na(Y[j,])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]

    if(j%%round(M/10)==0 && verbose)
      cat(".") # Report progress

    try({ # Try to fit a mixed model to adjust for plate
      if(!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
        OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
      } else
        OLS <- TRUE

      if(OLS) {
        fit <- lm(modelFix, data=pheno[ii,])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j,] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)

      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }
  if(verbose)
    cat(" done\n")
  ## Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1

  ## Get P values corresponding to F statistics
  Pval <- 1-pf(Fstat, sizeModel, degFree)

  out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
              sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
              orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
              degFree=degFree)

  out
}


pickCompProbes <- function(rawbetas, cellInd, cellTypes = NULL, numProbes = 50, probeSelect = probeSelectIN) {
  ## p is matrix of beta values
  ## cellInd is vector denoting cell type
  splitit <- function(x) {
    split(seq(along=x), x)
  }

  #   p <- getBeta(mSet)
  #   pd <- as.data.frame(colData(mSet))
  if(!is.null(cellTypes)) {
    if(!all(cellTypes %in% as.character(cellInd)))
      stop("elements of argument 'cellTypes' is not part of 'cellInd'")
    keep <- which(as.character(cellInd) %in% cellTypes)
    rawbetas <- rawbetas[,keep]
    cellInd<-cellInd[keep]
  }
  ## make cell type a factor
  cellInd <- factor(cellInd)
  ffComp <- rowFtests(rawbetas, cellInd)
  prof <- sapply(splitit(cellInd), function(i) rowMeans(rawbetas[,i]))
  r <- matrixStats::rowRanges(rawbetas)
  compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range")
  tIndexes <- splitit(cellInd)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(rawbetas))
    x[i] <- 1
    return(rowttests(rawbetas, factor(x)))
  })

  if (probeSelect == "any"){
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[,"p.value"] < 1e-8,]
      yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]
      c(rownames(yAny)[1:(numProbes*2)])
    })
  } else {
    probeList <- lapply(tstatList, function(x) {
      y <- x[x[,"p.value"] < 1e-8,]
      yUp <- y[order(y[,"dm"], decreasing=TRUE),]
      yDown <- y[order(y[,"dm"], decreasing=FALSE),]
      c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
  }

  trainingProbes <- na.omit(unique(unlist(probeList)))
  rawbetas <- rawbetas[trainingProbes,]

  pMeans <- colMeans(rawbetas)
  names(pMeans) <- cellInd

  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(cellInd), collapse="+")))
  phenoDF <- as.data.frame(model.matrix(~cellInd-1))
  colnames(phenoDF) <- sub("cellInd", "", colnames(phenoDF))
  if(ncol(phenoDF) == 2) { # two group solution
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(rawbetas))
  } else { # > 2 group solution
    tmp <- validationCellType(Y = rawbetas, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }

  out <- list(coefEsts = coefEsts, compTable = compTable,
              sampleMeans = pMeans)
  return(out)
}
# 
# 
# ##  INPUT: Y - raw data subset by CpGs in coefCellType
# ##         coefCellType - output of pickCompProbes containing the predictive CpGs and column per cell type
# 
# ## OUTPUT: predicted value per sample
# 
# projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
#   if(is.null(contrastCellType))
#     Xmat <- coefCellType
#   else
#     Xmat <- tcrossprod(coefCellType, contrastCellType) 
#   
#   nCol <- dim(Xmat)[2]
#   if(nCol == 2) {
#     Dmat <- crossprod(Xmat)
#     mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
#     colnames(mixCoef) <- colnames(Xmat)
#     return(mixCoef)
#   } else {
#     nSubj <- dim(Y)[2]
#     
#     mixCoef <- matrix(0, nSubj, nCol)
#     rownames(mixCoef) <- colnames(Y)
#     colnames(mixCoef) <- colnames(Xmat)
#     
#     if(nonnegative){
#       if(lessThanOne) {
#         Amat <- cbind(rep(-1, nCol), diag(nCol))
#         b0vec <- c(-1, rep(0, nCol))
#       } else {
#         Amat <- diag(nCol)
#         b0vec <- rep(0, nCol)
#       }
#       for(i in 1:nSubj) {
#         obs <- which(!is.na(Y[,i])) 
#         Dmat <- crossprod(Xmat[obs,])
#         mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
#       }
#     } else {
#       for(i in 1:nSubj) {
#         obs <- which(!is.na(Y[,i])) 
#         Dmat <- crossprod(Xmat[obs,])
#         mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
#       }
#     }
#     return(mixCoef)
#   }
# }
# 
# 
# #### Dorothea's edits
# 
# ## load packages
# library(ggplot2)
# # install.packages("ggfortify")
# library(ggfortify)
# # install.packages("cowplot")
# library(cowplot)
# # install.packages("ggrepel")
# library(ggrepel)
# # BiocManager::install("minfi")
# library(minfi)
# # BiocManager::install("genefilter")
# library(genefilter)
# 
# 
# 
# ## create function that randomly combines +ve and -ve samples at set proportions
# 
# ##  INPUT: posCellProp - vecor of the NeuN+ proportion to be simulated. Within [0,1]
# ##         betaList - List of 2 objects, +ve and -ve NeuN
# ##         nSamplesOut - number of simulated samples to be returned per proportion
# 
# ## OUTPUT: simSamples - matrix of betas whos cols match the entries in posCellProp
# 
# simulateProportionsF1 = function(posCellProp, betaList, nSamplesOut){
#   simSamples = matrix(ncol = nSamplesOut, nrow = nrow(betaList[[1]]))
#   
#   for (nSample in  1:nSamplesOut){
#     posSample = betaList[[1]][,sample(1:ncol(betaList[[1]]), 1)]
#     negSample = betaList[[2]][,sample(1:ncol(betaList[[2]]), 1)]
#     
#     simSamples[,nSample] = posSample*posCellProp + negSample*(1 - posCellProp)
#     
#   }
#   return(simSamples)
# }
# 
# 
# 
# 
# ### wrapper function to format the output of simulateProportions
# 
# ##  INPUT: posCellProp - vecor of the NeuN+ proportion to be simulated. Within [0,1]
# ##         betaList - List of 2 objects, +ve and -ve NeuN
# ##         nSamplesOut - number of simulated samples to be returned per proportion
# 
# ## OUTPUT: simBetasMatrix - nrow = nCpG, ncol = nSamplesOut*length(posCellProp)
# ##         propVec - vector of cell proportions per output column 
# 
# 
# simulateProportions = function(posCellProp, betaList, nSamplesOut){
#   
#   simSamples = sapply(posCellProp, simulateProportionsF1, betaList = betaList, nSamplesOut = nSamplesOut)
#   
#   propVec = rep(posCellProp, each = nSamplesOut)
#   
#   
#   
#   ## reshape the data so that each simulated sample has it's own column
#   simBetasMatrix = matrix(nrow = nrow(betaList[[1]]), ncol = length(posCellProp)*nSamplesOut)
#   for (i in 1:length(posCellProp)){
#     for (j in 1:nSamplesOut){
#       simBetasMatrix[ ,(i - 1)*nSamplesOut + j] = simSamples[(nrow(betaList[[1]]) * (j - 1) + 1):(nrow(betaList[[1]]) * (j)),i]
#     }
#   }
#   rownames(simBetasMatrix) = rownames(betaList[[1]])
#   
#   return(list(simBetasMatrix, propVec))
# }
# 
# 
# ## make function where the number of probes can be varied
# 
# ##  INPUT: Number of probes to be used in the model per cell type
# ##         
# 
# ## OUTPUT: R squared of the model
# 
# checkNCpGsForModel = function(numProbes){
#   rawbetas=as.matrix(betas[,trainInd])
#   rownames(rawbetas) = rownames(betas)
#   
#   compData <- pickCompProbes(rawbetas = rawbetas,
#                              as.factor(as.character(phenoNuc[trainInd, 12])), 
#                              cellTypes = levels(as.factor(as.character(phenoNuc[trainInd, 12]))), 
#                              numProbes = numProbes, 
#                              probeSelect = "auto")
#   
#   coefs <- compData$coefEsts
#   
#   testBetas=as.matrix(betas[,-trainInd])
#   rownames(testBetas) = rownames(betas)
#   testCounts <- projectCellType(testBetas[rownames(coefs), ], coefs)
#   
#   # calculate R^2
#   return(cor(testCounts[,2], as.numeric(phenoNuc[-trainInd,'Tissue'])))
# }
# 
# 
# 
# 
# ### Function to build models in train data with varying numbers of CpGs
# 
# ##  INPUT: posNCpG - vector of nCpGs to use in the model per cell type
# ##         betas - full betas containing both train and test data
# ##         pheno - containing columns $TrainTest and $Tissue
# 
# ## OUTPUT: neuNPosPRED - predicted NeuN positive from each model 
# ##         nueNPosTRUE - actual proportion of NeuN+
# ##         rSquared - ncol = 2, rSquared and nCpG
# ##         plot(rSquared[,1], rSquared[,2])
# 
# ## function building data
# # posNCpG = seq(1, 10)
# # betas = betas
# # pheno = phenoNuc
# 
# 
# 
# buildModelWithVaryingNumbersOfCpGs = function(posNCpG, betas, pheno, nSamplesOut = 500){
#   
#   ## subset train data
#   train = as.matrix(betas[,which(pheno$TrainTest == "Train")])
#   rownames(train) = rownames(betas)
#   phenoTrain = pheno[which(pheno$TrainTest == "Train"), ]
#   
#   ## create validation data and merge with actual purified test data
#   test = as.matrix(betas[ ,which(pheno$TrainTest == "Test")])
#   rownames(test) = rownames(betas)
#   phenoTest = pheno[which(pheno$TrainTest == "Test"), ]
#   
#   ## create inputs for simulateProportions
#   betaList = list(test[, which(phenoTest$Tissue == "nuclei_neuN_pos")],
#                   test[, which(phenoTest$Tissue == "nuclei_neuN_neg")])
#   posCellProp = seq(0.1, 0.9, 0.1)
#   
#   ## simulate data
#   simDat = simulateProportions(posCellProp, betaList, nSamplesOut)
#   simProp = simDat[[2]]
#   simDat = simDat[[1]]
#   
#   ## merge simulated data matrix with real thing and add real Neun+ values to [[2]]
#   simDatTest = cbind(simDat, test)
#   x = matrix(nrow = nrow(phenoTest), ncol = 1)
#   x[which(phenoTest$Tissue == "nuclei_neuN_pos"),] = 1
#   x[which(phenoTest$Tissue == "nuclei_neuN_neg"),] = 0
#   nueNPosTRUE = c(simProp, x)
#   
#   coefData = list()
#   nueNPosPRED = matrix(nrow = ncol(simDatTest), ncol = length(posNCpG))
#   for(nCpG in posNCpG){
#     
#     ## create models for each number of CpGs
#     coefData[[nCpG]] <- pickCompProbes(rawbetas = train,
#                                        as.factor(as.character(phenoTrain$Tissue)), 
#                                        cellTypes = levels(as.factor(as.character(phenoTrain$Tissue))), 
#                                        numProbes = nCpG, 
#                                        probeSelect = "auto")$coefEsts
#     
#     ## apply models to the testing data
#     nueNPosPRED[,nCpG] = projectCellType(simDatTest[rownames(coefData[[nCpG]]), ], 
#                                          coefData[[nCpG]])[, 2] # 2 is NeuN+
#   }
#   
#   ## calculate R Squared for each model 
#   rSquaredPerCol = function(nueNPosPREDcol, nueNPosTRUE){
#     return(cor(nueNPosPREDcol, nueNPosTRUE)^2)
#   }
#   
#   rSquaredNCpG = cbind(apply(nueNPosPRED, 2, rSquaredPerCol, nueNPosTRUE), 
#                        posNCpG)
#   colnames(rSquaredNCpG) = c("rSquared", "posNCpG")
#   
#   
#   print(ggplot(data = rSquaredNCpG, aes(x = posNCpG, y = rSquared)) +
#           geom_point() +
#           theme_cowplot(18) +
#           xlab("Number of CpGs per cell type") +
#           ylab("R squared"))
#   
#   return(list(nueNPosPRED, nueNPosTRUE, rSquaredNCpG))
#   
# }


library(quadprog)




projectCellTypeWithError = function(YIN, coefCellTypeIN, sampleDup = 0){
  
  ## if there's only 1 sample, pretend there are 2 so it doesn't break
  if (ncol(as.matrix(YIN)) == 1){ 
    sampleDup = 1
    YIN = cbind(YIN, YIN)
  }
  
  
  projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
      Xmat <- coefCellType
    else
      Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
      obs <- which(apply(Y, 1, function(x){return(sum(is.na(x)) == 0)}))
      Dmat <- crossprod(Xmat[obs,])
      mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat[obs,], x[obs])) }))
      colnames(mixCoef) <- colnames(Xmat)
      return(mixCoef)
    } else {
      nSubj <- dim(Y)[2]
      
      mixCoef <- matrix(0, nSubj, nCol)
      rownames(mixCoef) <- colnames(Y)
      colnames(mixCoef) <- colnames(Xmat)
      
      if(nonnegative){
        if(lessThanOne) {
          Amat <- cbind(rep(-1, nCol), diag(nCol))
          b0vec <- c(-1, rep(0, nCol))
        } else {
          Amat <- diag(nCol)
          b0vec <- rep(0, nCol)
        }
        for(i in 1:nSubj) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
        }
      } else {
        for(i in 1:nSubj) {
          obs <- which(!is.na(Y[,i])) 
          Dmat <- crossprod(Xmat[obs,])
          mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
      }
      return(mixCoef)
    }
  }
  
  ## error function
  getErrorPerSample = function(applyIndex,
                               predictedIN = mixCoef,
                               coefDataIN = coefCellTypeIN,
                               betasBulkIN = YIN){
    
    trueBulk = matrix(ncol = 1, nrow = nrow(coefDataIN), data = 0)
    
    RMSE = function(m, o){
      sqrt(mean((m - o)^2))
    }
    
    for (i in 1:ncol(coefDataIN)){
      
      trueBulk[,1] = trueBulk[,1] + coefDataIN[,i]*predictedIN[applyIndex,i]
    }
    
    betasBulkIN = t(apply(betasBulkIN, 1, function(x){x[is.na(x)] = 0; return(x)}))
    
    error = RMSE(trueBulk, betasBulkIN[,applyIndex])
    return(error)
  }
  
  mixCoef = projectCellType(YIN, coefCellTypeIN)
  error = sapply(1:nrow(mixCoef), getErrorPerSample)
  nCGmissing = apply(YIN, 2, function(x){sum(is.na(x))})
  mixCoef = cbind(mixCoef, error, nCGmissing)
  
  if (sampleDup == 1){
    mixCoef = mixCoef[1,]
  }
  
  return(mixCoef)
  
}
