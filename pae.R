##this script is to calculate and plot the distance variation matrix for all alpha carbons

#library(bio3d)
library(gplots)

pae <- read.csv("pae.csv", sep="", header=F)

dist.variant <- as.matrix(pae)

labs <- c(1,rep(NA,98),100,rep(NA,99),200,rep(NA,99),300,rep(NA,99),400,
          rep(NA,99),500, rep(NA,99),600, rep(NA,99),700, rep(NA,22))

max(dist.variant)

cols <- colorRampPalette(c("darkgreen","lightgreen", "white"))(n=256)

png(paste("E", ".pae.", "heatmap2.png", sep=""), 
    width=6, height=7.5, res=600, units="in")
heatmap.2(dist.variant, col=cols, trace='none',
          key.xlab=NA, key.title=NA, key.ylab=NA,
	  xlab="Residue No.", ylab="Residue No.",
	  labRow=labs, labCol=labs,
	  srtCol=45,
	  cexRow = 1,
	  cexCol = 1,
          key.xtickfun = NULL, key.ytickfun = NULL, Rowv = F, Colv = F,
          dendrogram = "none",
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lwid = c(1.5,4), lhei=c(1.5,4,1))
dev.off()

