##this script is to calculate and plot the distance variation matrix for all alpha carbons

library(bio3d)
library(gplots)

# load the trajectory (dcd) and protein (pdb)
dcd <- read.dcd("protein.rst1.dcd")
pdb <- read.pdb("protein.pdb") 

#total number of atoms (both pdb and dcd)
total.an <- dim(dcd)[2]
#frame number
frame.number <- dim(dcd)[1]

# all alpha carbon indices
ca.inds <- atom.select(pdb, "calpha")

#number of alpha carbons
an <- length(ca.inds$atom)

#define the matrix for distance variant
#create a matrix for distance variations (for those among alpha carbons)
#derived from all 100 frames in the trajectory (condensed!)
dist.variant <- matrix(nrow=an, ncol=an)

#the array for alpha-alpha distances
dist.array <- array(dim=c(frame.number/100, an, an))

#calculate distance every 100 frame (10k/100), set a random number
#the sampling should be repeatable
#the random number varies between 0 to 100
#as every 100 frame would be chosen from a 10k traj
#the statistics is for 100 sequential frames
rdm <- sample(c(0:99))[1]

# loop, note the random number may vary in each run
for (i in 1:(frame.number/100)) {
  for (k in 1:an) {
    kserial <- ca.inds$atom[k]
    kxcor <- dcd[i*100-rdm, kserial*3-2]
    kycor <- dcd[i*100-rdm, kserial*3-1]
    kzcor <- dcd[i*100-rdm, kserial*3]
    kcoor <- c(kxcor, kycor, kzcor)
    for (j in 1:an) {
      jserial <- ca.inds$atom[j]
      jxcor <- dcd[i*100-rdm, jserial*3-2]
      jycor <- dcd[i*100-rdm, jserial*3-1]
      jzcor <- dcd[i*100-rdm, jserial*3]
      jcoor <- c(jxcor, jycor, jzcor)
      dist.array[i, k, j] <- dist.xyz(kcoor, jcoor)
    }
  }
}

# uncertainty (error) for the alpha-alpha distances
# median * IQR
# because the 1D simplyficatin ignore orientational differences (angles)
# this may be compensated by the weighting by median radius
# Another option is the maximal - minimal (but there might be outlier biases)
for (s in 1:an){
  for (t in 1:an){
#      dist.variant[s,t] = round(sqrt(median(dist.array[,s,t])*IQR(dist.array[,s,t])),2)
    dist.variant[s,t] = round(IQR(dist.array[,s,t]),2)
  }
}

#normalization -- it does not contribute to the heatmap appearance
#dist.normalize <- (dist.variant-min(dist.variant))/(max(dist.variant)-min(dist.variant))

resid <- pdbseq(pdb)
resname <- c()
for (i in 1:length(resid)) {
	resname <- rbind(resname, paste(resid[i], i , sep=""))
}
#the residue names can be used as the row/col names
write.table(dist.variant, file="DV.matrix.csv", sep=",", row.names=resname, col.names=resname, quote=F)

cols <- colorRampPalette(c("darkgreen", "white"))(n=256)

labs <- c(1,rep(NA,98),100,rep(NA,99),200,rep(NA,99),300,rep(NA,99),400,
          rep(NA,99),500, rep(NA,99),600, rep(NA,99),700, rep(NA,22))

png(paste("E", ".dv.iqr.", rdm, ".heatmap2.png", sep=""), 
    width=6, height=7.5, res=600, units="in")
heatmap.2((dist.variant), col=cols, trace='none',
	  xlab="Residue No.", ylab="Residue No.",
	  srtCol=45,
	  cexCol=1, cexRow=1,
	  labRow=labs, labCol=labs,
          key.xlab=NA, key.title=NA, key.ylab=NA,
          key.xtickfun = NULL, key.ytickfun = NULL, Rowv = F, Colv = F,
          dendrogram = "none",
          lmat=rbind(c(0,3),c(2,1),c(0,4)), lwid = c(1.5,4), lhei=c(1.5,4,1))
dev.off()
#normalization does not affect the appearance!
max(dist.variant) # 16.85 Angs.

pae <- read.csv("pae.csv", sep="", header=F)

pae.mat <- as.matrix(pae)

max(pae.mat)

cor.test(c(dist.variant), c(pae.mat))
