library(bio3d)

dcd <- read.dcd("protein.rst1.dcd")
pdb <- read.pdb("protein.pdb") #dimer

#plot a Ramachaindrion plot
tor <- torsion.pdb(pdb)
pdf("rama.pdf", height=6, width=6, paper='special')
plot(tor$phi, tor$psi, xlab="Phi", ylab="Psi", col="blue")
dev.off()

#set the alpha carbons and backbone atoms
ca.inds <- atom.select(pdb, "calpha")
bb.inds <- atom.select(pdb, elety=c("N", "CA", "C", "O", "OT1", "OT2"))
xyz <- fit.xyz(fixed = pdb$xyz, mobile = dcd,
               fixed.inds = bb.inds$xyz,
               mobile.inds = bb.inds$xyz)

#plot the RMSD
pdf("rmsd.pdf", width=8, height=4, paper='special')
rd <- rmsd(xyz[1,bb.inds$xyz], xyz[, bb.inds$xyz])
plot(rd, type="l", ylab="RMSD (Angs.)", xlab="Time (ns)")
points(lowess(rd), type="l", col="red", lty=2, lwd=2)
dev.off()

#plot RMSF
rf <- rmsf(xyz[,ca.inds$xyz])
pdf("rmsf.pdf", width=8, height=4, paper='special')
#plot(rf, ylab="RMSF (Angs.)", xlab="Residue Position", type="l")
plot.bio3d(rf, resno=pdb, type="l", lwd=1.5,
	  col="red", ylab="RMSF (Angs.)")
#points(rf, typ="l", col="red")
dev.off()
write.table(rf, file="protein.rmsf.dat", row.names=F, col.names=F, quote=F)

#PCA over backbone atoms
pc.traj <- pca.xyz(xyz[, bb.inds$xyz], mass=atom2mass(pdb,bb.inds))
pdbseq(pdb)

seq <- c()
for (i in 1:length(pdbseq(pdb))) {
	seq <- c(seq, rep(i,4))
}
seq <- c(seq, length(pdbseq(pdb)))

resid <- c()
resname <- aa123(pdbseq(pdb))
for (j in 1:length(pdbseq(pdb))) {
	resid <- c(resid, rep(resname[j], 4))
}
resid <- c(resid, resname[length(pdbseq(pdb))])

resid <- c(resid, resid)

an <- c(rep(c("N", "CA", "C", "O"), length(pdbseq(pdb))-1), 
        c("N", "CA", "C", "OT1", "OT2"))

#top three trajectories
mktrj.pca(pc.traj, pc=1, file="pc1.pdb", elety = an, resno = seq, resid = resid)
mktrj.pca(pc.traj, pc=2, file="pc2.pdb", elety = an, resno = seq, resid = resid)
mktrj.pca(pc.traj, pc=3, file="pc3.pdb", elety = an, resno = seq, resid = resid)

#print top 6
print(pc.traj, nmodes=6)

#residue cross-correlation
cij <- dccm(xyz[,ca.inds$xyz])

png("residue.correlation.png", width=6, height=5, res=600, units="in")
plot(cij)
dev.off()

#visualize cross-correlation using pymol
pymol.dccm(cij, pdb, type="launch")

