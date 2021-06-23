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

png("/mnt/data1/Thea/ErrorMetric/plots/ValidateInitialModel/brainBloodPredicted.png",
    height = 550, width = 600)
ggplot(dat, aes(x = error.x, y = error.y, col = Tissue.x)) +
  geom_point(size = 2.5, shape = 21) +
  geom_hline(yintercept = 0.1, col = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.1, col = "red", linetype = "dashed") +
  theme_cowplot(18) +
  labs(x = "Cetygo in blood RBDM", y = "Cetygo in brain RBDM", col = "Tissue")
dev.off()

## get sample sizes of each tissue in each quadrant
head(dat)
sum(dat$Tissue.x[dat$error.x <0.1] == "Blood")
sum(dat$Tissue.x[dat$error.x <0.1] == "Brain")

sum(dat$Tissue.x[dat$error.y <0.1] == "Blood")
sum(dat$Tissue.x[dat$error.y <0.1] == "Brain")

sum(dat$Tissue.x[dat$error.y >0.1 & dat$error.x >0.1] == "Blood")
sum(dat$Tissue.x[dat$error.y >0.1 & dat$error.x >0.1] == "Brain")

# error.x = blood error
# error.y = brain error

## proportion of the dataset that is 'brain'
blOnly = dat[dat$Tissue.x == "Blood",]
blOnly$DatasetOrigin.x = as.factor(as.character(blOnly$DatasetOrigin.x))

## blood predicted as brain
blAsBr = blOnly[which(blOnly$error.y <0.1 & blOnly$Tissue.x == "Blood"),]

propWrongBrainTissue = 1-(table(blOnly$DatasetOrigin.x)-table(blAsBr$DatasetOrigin.x))/table(blOnly$DatasetOrigin.x)
propWrongBrainTissue = propWrongBrainTissue[propWrongBrainTissue!=0]

propWrongBrainTissue = cbind.data.frame(propWrongBrainTissue,
                                     totWrong = table(blAsBr$DatasetOrigin.x)[table(blAsBr$DatasetOrigin.x)!=0])[,c(2,1,4)]
colnames(propWrongBrainTissue) = c("prop", "dat", "n")
propWrongBrainTissue$int = 0
propWrongBrainTissue$prop = as.numeric(as.character(propWrongBrainTissue$prop))

library(viridis)
ggplot(propWrongBrainTissue, aes(y = prop)) +
  geom_boxplot() +
  geom_jitter(aes(x = int,col = n)) +
  scale_color_viridis() +
  theme_cowplot(18) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "Proprotion of blood samples\npredicted as brain", col = "Number\nof samples")



## look at Erisk as a whole: How many samples? Are they removed in QC?
eriskwb = (dat[dat$DatasetOrigin.x == "GSE105018.gds",])


## load erisk pheno file 
eriskwbPheno = read.csv("/mnt/data1/E-risk/FinalQC_Liv/All_chips94_95_Samplesheet.csv", skip = 3, header = T)
eriskwbPheno$Barcode = paste(eriskwbPheno$Chip, eriskwbPheno$Position, sep = "_")

sum(eriskwb$barcode %in% eriskwbPheno$Barcode)

sum(blAsBr$barcode %in% eriskwbPheno$Barcode)



















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




### Save the data ####
## load brain
# load("/mnt/data1/Thea/humanDeconvolution/data/CETSTrainTest.RData")
# load("/mnt/data1/Thea/ErrorMetric/DSRMSE/models/CETSmodel50CpG.Rdata")
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


## save datasets real quick
save(brainDat, bloodDat, file = "/mnt/data1/Thea/ErrorMetric/data/EssexBrainBloodPredicted.Rdata")


# #### Check cell type summing for well predicted non brain/blood tissues ####
load(file = "/mnt/data1/Thea/ErrorMetric/data/EssexBrainBloodPredicted.Rdata")

# ## add sum column to each cell model
# brainDat$sum = abs(brainDat$`NeuN+` + brainDat$`NeuN-` -1)
# bloodDat$sum = abs(bloodDat$Bcell + bloodDat$CD4T + bloodDat$CD8T + bloodDat$Gran + bloodDat$Mono + bloodDat$NK -1)
# 
# ## subset both for only those with good error
# brainDat = brainDat[brainDat$error <0.1,]
# bloodDat = bloodDat[bloodDat$error <0.1,]
# 
# ## boxplot the difference in sub between actual and other tissue
# library(ggplot2)
# library(cowplot)
# library(forcats)
# library(dplyr)
# 
# ## blood
# 
# bloodDat$Tissue = as.character(bloodDat$Tissue)
# bloodDat[bloodDat$blood ==1, "Tissue"] = "Blood"
# bloodDat$Tissue = as.factor(bloodDat$Tissue)
# 
# dat_summary = bloodDat %>%
#   group_by(Tissue) %>%
#   tally()
# 
# dat = merge(bloodDat, dat_summary,  by = "Tissue")
# 
# dat$Tissue = as.factor(as.character(dat$Tissue))
# pos = c()
# for(i in 1:length(levels(dat$Tissue))){
#   pos = c(pos, max(dat$sum[dat$Tissue == levels(dat$Tissue)[i]]))
# }
# 
# dat.pos = data.frame(Tissue = levels(dat$Tissue), pos, n = dat_summary$n)
# 
# 
# ggplot(bloodDat, aes(x = fct_reorder(Tissue, blood, .fun = median, .desc =TRUE))) +
#   geom_boxplot(aes(y = sum, fill = as.factor(blood))) +
#   theme_cowplot(18) +
#   scale_fill_manual(values = c("#0A8ABA", "#BA3A0A"), name = "Blood?", labels = c("No", "Yes")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = element_blank(), y = "Cell type sum - 1") +
#   geom_text(data = dat.pos, aes(Tissue, label = n, y = pos+0.02))
# 
# 
# 
# ## brain
# brainDat$Tissue = as.character(brainDat$Tissue)
# brainDat[brainDat$brain ==1, "Tissue"] = "Brain"
# brainDat$Tissue = as.factor(brainDat$Tissue)
# 
# dat_summary = brainDat %>%
#   group_by(Tissue) %>%
#   tally()
# 
# dat = merge(brainDat, dat_summary,  by = "Tissue")
# 
# dat$Tissue = as.factor(as.character(dat$Tissue))
# pos = c()
# for(i in 1:length(levels(dat$Tissue))){
#   pos = c(pos, max(dat$sum[dat$Tissue == levels(dat$Tissue)[i]]))
# }
# 
# dat.pos = data.frame(Tissue = levels(dat$Tissue), pos, n = dat_summary$n)
# 
# 
# ggplot(brainDat, aes(x = fct_reorder(Tissue, brain, .fun = median, .desc =TRUE))) +
#   geom_boxplot(aes(y = sum, fill = as.factor(brain))) +
#   theme_cowplot(18) +
#   scale_fill_manual(values = c("#0A8ABA", "#BA3A0A"), name = "Brain?", labels = c("No", "Yes")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = element_blank(), y = "Cell type sum - 1") +
#   geom_text(data = dat.pos, aes(Tissue, label = n, y = pos+0.02))































































































































































 















