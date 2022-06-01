## create plots for manuscript from saved results files
## assumes executed from CETYGOAnalyses folder

library(Rmisc)
library(ggplot2)
library(ggpubr)

#### Plot results from simulations with noise

load("RData/ResultsNoisyProfilesWithinBatch.rdata")
sumpredProp <- summarySE(predPropAll, measurevar="RMSE", groupvars=c("noise", "Array"))
figa1<-ggplot(sumpredProp, aes(x=noise, y=RMSE, colour=Array)) + 
    geom_line() +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Proportion noise") +
    ylab("RMSE of prediction") +
    expand_limits(y=0) +                        # Expand y range
    theme_bw() +
    theme(legend.position="none")  +             # Position legend in bottom right
    geom_ribbon(aes(ymin=RMSE-ci, ymax=RMSE+ci, colour=Array), linetype=2, alpha=0.1)

sumpredProp <- summarySE(predPropAll, measurevar="CETYGO", groupvars=c("noise", "Array"))
figa2<-ggplot(sumpredProp, aes(x=noise, y=CETYGO, colour=Array)) + 
    geom_line() +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Proportion noise") +
    ylab("CETYGO") +
    expand_limits(y=0) +                        # Expand y range
    theme_bw() +
	theme(legend.position = "none") +
    geom_ribbon(aes(ymin=CETYGO-ci, ymax=CETYGO+ci, colour=Array), linetype=2, alpha=0.1)

sumpredProp <- summarySE(predPropAll, measurevar="TotalComposition", groupvars=c("noise", "Array"))
figa3<-ggplot(sumpredProp, aes(x=noise, y=TotalComposition, colour=Array)) + 
    geom_line() +
    geom_point(size=1.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Proportion noise") +
    ylab("Sum of cellular composition") +
	theme_bw() + theme(legend.position = "right") +
    geom_ribbon(aes(ymin=TotalComposition-ci, ymax=TotalComposition+ci, colour=Array), linetype=2, alpha=0.1)

figa4<-ggplot(predPropAll, aes(x=TotalComposition, y=CETYGO, colour=Array)) + 
    geom_point(size=1.5, shape=21) + 
	geom_smooth(method = "loess", se = TRUE) +
    ylab("CETYGO") +
    xlab("Sum of cellular composition") +
    expand_limits(y=0) +
	theme_bw() +
	theme(legend.position = "none")
	
ggarrange(figa1, figa2, figa3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1, widths = c(4,4,5.5))
ggsave("plots/LineGraphCETYGONoisyModelsWithinBatch.pdf", width = 24, height = 8, units = "cm")
		  
### compare results across batches
load("RData/ResultsNoisyProfilesAcrossBatch.rdata")

	
sumpredProp <- summarySE(predPropArrays, measurevar="RMSE", groupvars=c("noise", "Array", "Normalised"))
supfiga1<-ggplot(sumpredProp, aes(x=noise, y=RMSE, colour=Array, shape = Normalised)) + 
    geom_line() +
    geom_point() + 
    xlab("Proportion noise") +
    ylab("RMSE of prediction") +
    expand_limits(y=0) +                        # Expand y range
    theme_bw() +
	theme(legend.position = "none")

sumpredProp <- summarySE(predPropArrays, measurevar="CETYGO", groupvars=c("noise", "Array", "Normalised"))
supfiga2<-ggplot(sumpredProp, aes(x=noise, y=CETYGO, colour=Array, shape = Normalised)) + 
    geom_line() +
    geom_point() + # 21 is filled circle
    xlab("Proportion noise") +
    ylab("CETYGO") +
    expand_limits(y=0) +                        # Expand y range
    theme_bw() +
	theme(legend.position = "none")
	
sumpredProp <- summarySE(predPropArrays, measurevar="TotalComposition", groupvars=c("noise", "Array", "Normalised"))
supfiga3<-ggplot(sumpredProp, aes(x=noise, y=TotalComposition, colour=Array, shape = Normalised)) + 
    geom_line() +
    geom_point() + # 21 is filled circle
    xlab("Proportion noise") +
    ylab("Sum of cellular composition") +
	theme_bw() +
	theme(legend.position = "right") 
	
supfiga4<-ggplot(predPropArrays, aes(x=TotalComposition, y=CETYGO, colour=Array, shape = Normalised)) + 
    geom_point() + 
	geom_smooth(method = "loess", se = TRUE) +
    ylab("CETYGO") +
    xlab("Sum of cellular composition") +
    expand_limits(y=0) +
	theme_bw() +
	theme(legend.position = "right")
	
ggarrange(supfiga1, supfiga2, supfiga3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1, widths = c(4,4,5.5))
ggsave("plots/LineGraphCETYGONoisyModelsAcrossBatch.pdf", width = 24, height = 8, units = "cm")
		  		  

		  
#### Plot results from simulations with incomplete models


load("RData/ResultsIncompleteModels.rdata")
colnames(predPropCombined)[1]<-"CETYGO"

predPropCombined<-as.data.frame(predPropCombined)

figb1 <- ggplot(predPropCombined, aes(x=factor(nCT), y=CETYGO)) + 
    geom_violin() + 
	geom_boxplot(width=0.1) +
	labs(y = "CETYGO", x = "Number of cell types")


figb2 <- ggplot(predPropCombined, aes(x=factor(propRepresented), y=CETYGO)) + 
    geom_violin() + 
	geom_boxplot(width=0.1) +
	labs(y = "CETYGO", x = "Proportion represented in model")



ggarrange(figb1, figb2,  
          labels = c("A", "B"),
          ncol = 2, nrow = 1, widths = c(1,2))
		  
ggsave("plots/ViolinPlotCETYGOIncompleteModels.pdf", width = 18, height = 8, units = "cm")

load("RData/ResultsLOOModels.rdata")

predPropLOO<-as.data.frame(predPropLOO)
predPropLOO$V1<-as.numeric(as.character(predPropLOO$V1))
colnames(predPropLOO)<-c("CETYGO", "PropMissing", "MissingCellType") 

sumpredProp <- summarySE(predPropLOO, measurevar="CETYGO", groupvars=c("PropMissing", "MissingCellType"))
sumpredProp$PropMissing<-as.numeric(as.character(sumpredProp$PropMissing))
figc1<-ggplot(sumpredProp, aes(x=PropMissing, y=CETYGO, colour=MissingCellType)) + 
    geom_line() +
    geom_point() + # 21 is filled circle
    xlab("Proportion of missing cell type") +
    ylab("CETYGO") +
    expand_limits(y=0) +
	theme_bw() +
	theme(legend.position="right", legend.title = element_text(size = 10), 
               legend.text = element_text(size = 10)) +
    geom_ribbon(aes(ymin=CETYGO-ci, ymax=CETYGO+ci, colour=MissingCellType), linetype=2, alpha=0.1)


ggsave("plots/LineGraphCETYGOLOOModels.pdf", width = 12, height = 10, units = "cm")

