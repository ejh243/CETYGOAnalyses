## for same set of samples compare CETYGO calculated
# 1.  after normalising ref data and test data together (using estimateCellCountsWithError())
# 2.  after normalisation separately (using projectCellTypeWithError())

args<-commandArgs(trailingOnly=TRUE)

idatPath<-args[1]
normPath<-args[2]


library(ggplot2)
library(ggpubr)
library(gridExtra)
library(CETYGO)
#library(minfi)
#library(quadprog)

load(normPath)
RGset <- read.metharray.exp(base = idatPath, targets = pheno, force = TRUE)


## compare cell type estimates across both methods
cellCompOrig<-estimateCellCounts(RGset)
cellCompAdpt<-estimateCellCountsWithError(RGset)


cellCompMatrix<-projectCellTypeWithError(betas, modelBloodCoef[rownames(modelBloodCoef) %in% rownames(betas),])

## compare predictions
cellTypes<-colnames(cellCompOrig)

p <- list()
for(each in cellTypes){
	dat<-data.frame("Original" = cellCompOrig[,each], "Adapted" = cellCompMatrix[, each]) 
	p[[each]] <- ggplot(dat, aes(x=Original, y=Adapted)) + 
    geom_point() +
	labs(y = "Normalised separately", x = "Normalised together", title = each)
}

do.call(grid.arrange,p)
ggsave("plots/PlotsCompareNormalisationStrategies.pdf", width = 20, height = 20, units = "cm")		  

dat<-data.frame("Original" = cellCompAdpt[,"error"], "Adapted" = cellCompMatrix[, "error"], "Sex" = pheno.pass$Sex, "Age" = pheno.pass$Age) 
figb <- ggplot(dat, aes(x=Original, y=Adapted, color = Sex)) + 
    geom_point() +
	labs(y = "Normalised separately", x = "Normalised together", title = "CETYGO")
figb 
ggsave("plots/PlotsCompareNormalisationStrategiesCETYGO.pdf", width = 20, height = 20, units = "cm")		  

	
figc1<-ggplot(dat, aes(x=Sex, y=Original)) + 
    geom_violin() + 
	geom_boxplot(width=0.1)+ 
	labs(y = "CETYGO", x = "", title = "Normalised together")
figc2<-ggplot(dat, aes(x=Sex, y=Adapted)) + 
    geom_violin() + 
	geom_boxplot(width=0.1)+ 
	labs(y = "CETYGO", x = "", title = "Normalised separately")
	

  
ggarrange(figc1, figc2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 2)
		  
ggsave("plots/ViolinplotsCompareNormalisationStrategiesBySex.pdf", width = 20, height = 20, units = "cm")		  
		  