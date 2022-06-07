## apply CETYGO to range of different datasets including blood and non-blood
## assumes executed from the CETYGOAnalyses folder

source("R/SmokingScoreFunction.r")

processDataset<-function(dataName, dataPath,  dataTissue, dataArray){
	## assumes dataPath contains two R objects betas and pheno
	## assumes in matched order
	## asummes pheno contains "Age" and "Sex"
	## checks for "Plate" and "Chip" columns
	## adds tissue column unless the name of a valid column is provided.
	load(dataPath)

	## check if Tissue column is present
	if(dataTissue %in% colnames(pheno)){
		pheno$Tissue<-pheno[,dataTissue]
	} else{
		pheno$Tissue<-dataTissue
	}

	## check Plate column is present
	if("Plate" %in% colnames(pheno)){
		paste0(dataName, pheno$Plate)
	}else{
		pheno$Plate<-NA
	}
	## check Chip column is present	
	if(sum(c("Chip", "ChipID", "MethArray_Chip") %in% colnames(pheno)) == 0){
		pheno$Chip<-NA
	}else{	
		if(sum(colnames(pheno) %in% c("ChipID", "MethArray_Chip")) > 0){
			colnames(pheno)[colnames(pheno) %in% c("ChipID", "MethArray_Chip")]<-"Chip"
		}
	}
	## calculate cell composition and cetygo
	rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
	predProp<-projectCellTypeWithError(betas, modelBloodCoef[rowIndex,])
	
	## calculate DNAmAge
	dnamage<-agep(betas)
	
	## calculate smoking score	
	smokeScore <- smokingscore(betas)
	
	return(data.frame("Dataset"=dataName, "Tissue"=pheno$Tissue, "Array" = dataArray,"CETYGO" = predProp[,"CETYGO"], "Age" = pheno$Age,  "DNAmAge" = dnamage,"Sex" = pheno$Sex, "Smoking Score" = smokeScore, "Plate" = pheno$Plate, "Chip" = pheno$Chip))
	

}


library(CETYGO)
library(quadprog)
library(wateRmelon)
library(ggplot2)
library(data.table)

realData<-NULL
datasets<-read.csv("RealDatasetsMetaData.csv", stringsAsFactors = FALSE)
for(i in 1:nrow(datasets)){
	realData<-rbind(realData, processDataset(dataName = datasets[i,"Dataset"], dataPath = datasets[i,"Path"], dataTissue = datasets[i,"Tissue"], dataArray = datasets[i,"Array"]))
}


## reformat strings
group<-rep("Blood derivative", nrow(realData))
group[realData$Tissue %in% c("Buccal","Nasal","Neocortex","Cerebellum", "Hippocampus", "Superior temporal gyrus", "Pancreas", "Liver")]<-"Non-blood"
realData$group<-factor(group, levels = c("Blood derivative", "Non-blood"))

realData$Tissue[which(realData$Tissue == "whole blood")]<-"WholeBlood"

realData$Tissue<-factor(realData$Tissue, levels= unique(realData[order(realData$group), "Tissue"]))

realData$Age<-as.numeric(realData$Age)

realData$Sex[which(realData$Sex == "")]<-NA
realData$Sex[tolower(realData$Sex) %in% c("male")]<-"Male"
realData$Sex[tolower(realData$Sex) %in% c("female")]<-"Female"
realData$Sex<-factor(realData$Sex)

realData$Smoking.Score<-as.numeric(realData$Smoking.Score)

save(realData, file = "RData/realData.rdata")
## plot as function of tissue type


figa1 <- ggplot(realData, aes(x=reorder(Tissue, group), y=CETYGO, fill = group)) + 
    geom_violin(scale = "width") + geom_boxplot(width = 0.1) + 
	labs(y = "CETYGO", x = "", fill = "") +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key.size = unit(0.5, 'cm'))
figa1
ggsave("plots/ViolinplotsCetygoByTissue.pdf", width = 20, height = 15, units = "cm")
	
realData.arisk<-subset(realData, Dataset == "A-Risk" &  group == "Blood derivative")
figa2 <- ggplot(realData.arisk, aes(x=Tissue, y=CETYGO, fill = group)) + 
    geom_violin(scale = "width") + geom_boxplot(width = 0.1) + 
	labs(y = "CETYGO", x = "", fill = "") +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")
figa2

ggsave("plots/ViolinplotsAriskCetygoByBloodTissue.pdf", width = 10, height = 10, units = "cm")
		
	
realData.blood<-subset(realData, Tissue == "WholeBlood")
realData.blood$Age2<-realData.blood$Age^2
	
figb1 <- ggplot(realData.blood, aes(x=Dataset, y=CETYGO, fill = Array)) + 
    geom_violin() + geom_boxplot(width = 0.1) +
	labs(y = "CETYGO", x = "Dataset")	+
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 45, margin = margin(t = 20)), legend.key.size = unit(0.5, 'cm')) + scale_x_discrete(
        labels=c("1", "2", "3", "4", "5", "6", "7", "8"))
figb1
ggsave("plots/ViolinplotsCetygoWBByDataset.pdf", width = 16, height = 8, units = "cm")

## sex differences?
model<-lm(CETYGO ~ Age + Age2 + Sex + Smoking.Score + Array + Dataset, data = realData.blood)
dataset.null<-lm(CETYGO ~ Age + Age2 + Sex + Smoking.Score + Array, data = realData.blood)
anova(model, dataset.null)

figb2 <- ggplot(subset(realData.blood, !is.na(Sex)), aes(x=Sex, y=CETYGO, fill = Array)) + 
    geom_violin() + 
	labs(y = "CETYGO", x = "")
	
figb3 <- ggplot(realData.blood, aes(x=Age, y=CETYGO, color = Array, alpha=I(0.25))) + 
    geom_point(size = 0.5) +
	labs(y = "CETYGO", x = "Age") + 
  geom_smooth(method='lm', formula= y~x + I(x^2)) +
  theme(legend.position="none")
	
figb4 <- ggplot(realData.blood, aes(x=Smoking.Score, y=CETYGO, color = Array, alpha=I(0.25))) + 
    geom_point(size = 0.5) +
	labs(y = "CETYGO", x = "Smoking score") + 
  geom_smooth(method='lm', formula= y~x) +
  theme(legend.position="none")

ggarrange(figb2, figb3, figb4, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
		  
ggsave("plots/PlotsCetygoWBAgainstBiologicalFactors.pdf", width = 25, height = 10, units = "cm")		  
		  
		  
## is cetygo indicative of greater "errors" on other epigenetic predictors?

realData.blood$deltaAge<- realData.blood$DNAmAge - realData.blood$Age

modelDA<-lm(deltaAge ~ CETYGO + Age + Age2 + Sex + Smoking.Score + Array + Dataset, data = realData.blood)

modelDAAbs<-lm(abs(deltaAge) ~ CETYGO + Age + Age2 + Sex + Smoking.Score + Array + Dataset, data = realData.blood)

figc <- ggplot(realData.blood, aes(x=abs(deltaAge), y=CETYGO)) + 
	labs(y = "CETYGO", x = "Delta Age")+
  geom_hex() + 
  geom_smooth(method='lm', formula= y~x) +
  theme_bw()
figc
ggsave("plots/HeatscatterplotsCetygoAgainstDeltaAge.pdf", width = 10, height = 8, units = "cm")	