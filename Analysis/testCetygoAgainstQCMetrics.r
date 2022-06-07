## compare cetygo to QC metrics
## reload dataset with samples that failed QC
## assumes executed from CETYGOAnalyses folder

args<-commandArgs(trailingOnly=TRUE)

idatPath<-args[1]
idatList<-args[2]

library(wateRmelon)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(CETYGO)
library(quadprog)


sampToLoad<-read.csv(idatList)
msetEPIC <- readEPIC(idatPath, barcodes=sampToLoad[,1])

m_intensities<-methylated(msetEPIC)
u_intensities<-unmethylated(msetEPIC)


M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

bs<-bscon(msetEPIC)

msetEPIC.pf <- pfilter(msetEPIC)
pFilterPass<-colnames(betas(msetEPIC)) %in% colnames(betas(msetEPIC.pf))

rm(msetEPIC.pf)
msetEPIC.dasen <- dasen(msetEPIC)
betas.norm<-betas(msetEPIC.dasen)

rm(msetEPIC.dasen)

## define good quality samples:
qcPass<-(pFilterPass & (M.median > 500 | U.median > 500) & bs > 80)

## calculate CETYGO with cell composition prediction

predProp<-as.data.frame(projectCellTypeWithError(betas.norm, modelBloodCoef))
predProp$M.median<-M.median
predProp$U.median<-U.median
predProp$bs<-bs
predProp$pFilterPass<-pFilterPass
predProp$qcPass<-qcPass


figa1 <- ggplot(predProp, aes(x=qcPass, y=CETYGO, colour=qcPass)) + 
    geom_violin() + 
	geom_boxplot(width=0.02)+
	labs(y = "CETYGO", x = "Pass QC", color = "Pass QC")
	
figa2<-ggplot(predProp, aes(x=M.median, y=CETYGO, colour=qcPass)) + 
    geom_point() +
    ylab("CETYGO") +
    xlab("Median M intensity") +
    expand_limits(y=0)+
    expand_limits(x=0) +
	theme_bw()  +
  labs(color = "Pass QC")

	
figa3<-ggplot(predProp, aes(x=U.median, y=CETYGO, colour=qcPass)) + 
    geom_point() +
    ylab("CETYGO") +
    xlab("Median U intensity") +
    expand_limits(y=0)+
    expand_limits(x=0) +
	theme_bw() +
  labs(color = "Pass QC")

	
		
figa4<-ggplot(predProp, aes(x=bs, y=CETYGO, colour=qcPass)) + 
    geom_point() +
    ylab("CETYGO") +
    xlab("Bisulfite conversion (%)") +
    expand_limits(y=0)+
    expand_limits(x=100) +
	theme_bw() +
  labs(color = "Pass QC")


ggarrange(figa1, figa2, figa3, figa4, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
		  
ggsave("plots/ScatterplotsCETYGOAgainstQCMetrics.pdf", width = 16, height = 16, units = "cm")

save(predProp, file = "RData/CETYGOPoorDataQuality.rdata")
