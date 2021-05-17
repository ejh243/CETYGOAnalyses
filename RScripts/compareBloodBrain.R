### Brain blood comparison

## load brain
load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
load("/mnt/data1/Thea/humanDeconvolution/data/cellTypeColours.Rdata")
source("/mnt/data1/Thea/ErrorMetric/RScripts/FunctionsForErrorTesting.R")
library(gdsfmt)
library(bigmelon)
gfile = openfn.gds("/mnt/data1/Thea/ErrorMetric/data/EssexOutput/brain.gds")

brainDat = cbind.data.frame(pData(gfile), read.gdsn(index.gdsn(gfile$root, "Pred")))
colnames(brainDat)[11:14] = c("NeuN-", "NeuN+", "error", "nCGmissing")

## add a brain binary
brainDat$brain = rep(0, nrow(brainDat))
brainDat$brain[brainDat$Tissue == "Brain" |
            brainDat$Tissue == "Neuron"] = 1

closefn.gds(gfile)

library(gdsfmt)
gfile = openfn.gds("/mnt/data1/Thea/ErrorMetric/data/EssexOutput/sub.gds")

bloodDat = cbind.data.frame(read.gdsn(index.gdsn(gfile$root, "Pred")),
                       Age = read.gdsn(index.gdsn(gfile$root, "Age")),
                       Sex = read.gdsn(index.gdsn(gfile$root, "Sex")),
                       Tissue = read.gdsn(index.gdsn(gfile$root, "Tissue")),
                       SubTissue = read.gdsn(index.gdsn(gfile$root, "SubTissue")),
                       Sample = read.gdsn(index.gdsn(gfile$root, "colnames")),
                       DatasetOrigin = read.gdsn(index.gdsn(gfile$root, "DatasetOrigin")))
load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/HousemanBloodModel50CpG.Rdata")
colnames(bloodDat)[1:8] = c(colnames(HousemanBlood50CpGModel$coefEsts), "error", "nCGmissing")

bloodDat$blood = rep(0, nrow(bloodDat))
bloodDat$blood[bloodDat$Tissue == "Blood" |
            bloodDat$Tissue == "B Cells" |
            bloodDat$Tissue == "Granulocyes" |
            bloodDat$Tissue == "Neutrophils" |
            bloodDat$Tissue == "NK" |
            bloodDat$Tissue == "Lymph Node" |
            bloodDat$Tissue == "T Cells"] = 1

closefn.gds(gfile)


### subset blood to blood and brain to brain
blood = bloodDat[which(bloodDat$blood == 1 | bloodDat$Tissue == "Brain"),]
brain = brainDat[which(brainDat$brain == 1 | brainDat$Tissue == "Blood"),]

library(dplyr)
dat = inner_join(blood, brain, by = c("Sample" = "FullBarCode"))


ggplot(dat, aes(x = error.x, y = error.y, col = Tissue.x)) +
  geom_point(size = 2, shape = 21) +
  geom_hline(yintercept = 0.1, col = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.1, col = "red", linetype = "dashed") +
  theme_cowplot(18) +
  labs(x = "Cetygo in blood RBDM", y = "Cetygo in brain RBDM", col = "Tissue")


unique(dat$DatasetOrigin.x[which(dat$error.y <0.1 & dat$Tissue.x == "Blood")])
badblood = dat[dat$Tissue.x == "Blood" & dat$error.y <0.1,]
badbrain = dat[dat$Tissue.x == "Brain" & dat$error.x <0.1,]


badblood$sums = badblood$"NeuN-" + badblood$"NeuN+"
badbrain$sums = apply(badbrain[,c("Bcell","CD4T","CD8T","Gran","Mono","NK")],1,sum)


plotDat = rbind.data.frame(badblood[,c("Tissue.x", "sums")],
                           badbrain[,c("Tissue.x", "sums")])

ggplot(plotDat, aes(x = Tissue.x, y = sums)) +
  geom_boxplot() +
  theme_cowplot(18) +
  scale_x_discrete(labels=c("Blood" = "Blood samples\npredicted as brain",
                            "Brain" = "Brain samples\npredicted as blood")) +
  labs(x = element_blank(), y = "Sum of cell type\nproportions in RBDM")



only01Brain = dat[which(dat$error.y<0.1 & dat$Tissue.x == "Blood"),]
apply(only01Brain[,c("Bcell","CD4T","CD8T","Gran","Mono","NK","NeuN-","NeuN+")],2,mean, na.rm = T)
apply(dat[which(dat$Tissue.x == "Blood"),c("Bcell","CD4T","CD8T","Gran","Mono","NK","NeuN-","NeuN+")],2,mean, na.rm = T)       
