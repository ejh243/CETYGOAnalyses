## compare unnorm model cpgs for RMSe between norm with all in train data

load("/mnt/data1/Thea/ErrorMetric/data/Houseman/unnormalisedBetasTrainTestMatrix.Rdata")
load("/mnt/data1/Thea/ErrorMetric/data/Houseman/quantileNormalisedBetasTrainTestMatrix.Rdata")

load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel150CpG.Rdata")
model = HousemanBlood150CpGModel

rm(list = setdiff(ls(), c("betasTrain", "quantileBetasTrain", "model")))

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")

modelBeta = GetModelCG(betasTrain, list(model))
qmodelBeta = GetModelCG(quantileBetasTrain, list(model))

dat = cbind(c(betasTrain), c(quantileBetasTrain))
dat = dat[complete.cases(dat),]
allrmse = RMSE(dat[,1], dat[,2])

modelrmse = RMSE(c(modelBeta), c(qmodelBeta))
