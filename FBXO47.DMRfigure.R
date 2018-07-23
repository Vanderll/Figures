rm(list=ls())
library(genemodel)
library(RColorBrewer)


######## DMR Results #############
#individual probe info
load(file="Y:/LaurenV_random/Norris/data/T1Dfigures/DMR.probeResults.v2.Rdata")
#BED file with comb-p output
dmrEPIC = read.table(file="Y:/LaurenV_random/Norris/data/T1Dfigures/mmt1d_filt_EPIC.anno.hg19.bed", sep="\t",header=TRUE)
dmr450 = read.table(file="Y:/LaurenV_random/Norris/data/T1Dfigures/mmt1d_filt.anno.hg19.bed", sep="\t",header=TRUE)

dmrEPIC.want = dmrEPIC[which(dmrEPIC$z_sidak_p<0.05 & dmrEPIC$n_probes>1),]
dmr450.want = dmr450[which(dmr450$z_sidak_p<0.05 & dmr450$n_probes>1),]

dmrEPIC.want = dmrEPIC.want[which(dmrEPIC.want$refGene_name %in% dmr450.want$refGene_name),]
dmr450.want = dmr450.want[which(dmr450.want$refGene_name %in% dmrEPIC.want$refGene_name),]

#### FBXO47 ###
#get the DMR range;
range450 = dmr450.want[which(dmr450.want$refGene_name=="FBXO47"),c("start", "end")]
rangeEPIC = dmrEPIC.want[which(dmrEPIC.want$refGene_name=="FBXO47"),c("start", "end")]

bothDMRrange = c(min(c(range450[,"start"], rangeEPIC[,"start"])), max(c(range450[,"end"], rangeEPIC[,"end"])))

##methylation (probe) results
FBXO47.meth = combined[grep("FBXO47", combined$UCSC_RefGene_Name),]

##gene structure info for plotting
FBXO47 = read.csv(file="Y:/LaurenV_random/Norris/data/T1Dfigures/FBXO47_structureHg19.csv")

##find what to adjust by (note: the "End" position is the lowest because FBXO47 is on the reverse strand)
adj = min(c(FBXO47.meth$Start, FBXO47$End))
hg19.start = adj

adj.zoom = bothDMRrange[1]
DMRrange = bothDMRrange[2] - bothDMRrange[1]

## Adjust the 2 datasets;
#first the gene structure info full view
FBXO47.adj = FBXO47
FBXO47.adj$Start = FBXO47.adj$Start-adj
FBXO47.adj$End = FBXO47.adj$End-adj

toScale = DMRrange/(max(FBXO47.adj$Start))
FBXO47.adj.scaled = FBXO47.adj
FBXO47.adj.scaled$Start = FBXO47.adj.scaled$Start*toScale
FBXO47.adj.scaled$End = FBXO47.adj.scaled$End*toScale


#gene structure zoomed view
FBXO47.adj.zoom = FBXO47
FBXO47.adj.zoom$Start = FBXO47.adj.zoom$Start-adj.zoom
FBXO47.adj.zoom$End = FBXO47.adj.zoom$End-adj.zoom


#then the methylation (zoomed)
FBXO47.meth.2 = data.frame(FBXO47.meth)
FBXO47.meth.2$xLoc = FBXO47.meth.2$pos - adj.zoom
FBXO47.meth.2 = FBXO47.meth.2[order(FBXO47.meth.2$xLoc),]

### adjust start site to include all CpG sites;
FBXO47.toPlot = data.frame(FBXO47.adj[,1], paste((FBXO47.adj[,2]), "-", (FBXO47.adj[,3]), sep=""))
colnames(FBXO47.toPlot) = c("type", "coordinates")

FBXO47.toPlot.zoom = data.frame(FBXO47.adj.zoom[,1], paste((FBXO47.adj.zoom[,2]), "-", (FBXO47.adj.zoom[,3]), sep=""))
colnames(FBXO47.toPlot.zoom) = c("type", "coordinates")
FBXO47.toPlot.zoom = FBXO47.toPlot.zoom[-which(FBXO47.adj.zoom$Start<0 | FBXO47.adj.zoom$End>DMRrange),]

FBXO47.toPlot.scaled = data.frame(FBXO47.adj.scaled[,1], paste((FBXO47.adj.scaled[,2]), "-", (FBXO47.adj.scaled[,3]), sep=""))
colnames(FBXO47.toPlot.scaled) = c("type", "coordinates")

### adjust DMR tracks;

range450.adj.zoom = range450-adj.zoom
rangeEPIC.adj.zoom = rangeEPIC-adj.zoom


### get the SNPs for a 3rd track;
FBXO47_snps = combined[grep("FBXO47", combined$UCSC_RefGene_Name), ]
FBXO47_snps_SBE = FBXO47_snps[!is.na(FBXO47_snps$SBE_rs),]
FBXO47_snps_Probe = FBXO47_snps[!is.na(FBXO47_snps$Probe_rs),]

FBXO47_snps.2 = rbind(FBXO47_snps_SBE, FBXO47_snps_Probe)
FBXO47_snps.2$toPlot = FBXO47_snps.2$pos-adj.zoom

#### Plot it;
#colors = brewer.pal("Set1", n=4)
colors = brewer.pal("Paired", n=12)
source("Y:/LaurenV_random/Norris/programs/T1Dpaper/genemodel.plot.LV.v3.R")
#par(mfrow=c(2,1))

png(file="Y:/LaurenV_random/Norris/data/T1Dfigures/FBXO47.v6.png", res=350, height=480*4, width=480*4)

plot.new()
par(mfrow = c(6,1)
    #par(mfrow = c(1,1))
    ,
    oma = c(0.05,1,0.05,0.05) + 0.1,
    mar = c(0.05,4,0.05,0.05) + 0.1)
layout(matrix(c(1,2,3,4,5,6), byrow=TRUE), widths = c(1,1,1,1,1,1), heights=c(1.5,0.4, 0.4, 0.4,0.4,1.25))

#Scatter Plot
plot(x=as.numeric(FBXO47.meth.2$xLoc), y=FBXO47.meth.2$Bcontrol.850, lwd=1.5, type="b", col=colors[1], pch=19, ylim=c(0,1), ylab="Beta",axes = F, main="", xlab='', xlim=c(0,max(FBXO47.meth.2$xLoc)))

lines(x=FBXO47.meth.2$xLoc, y=FBXO47.meth.2$Bcase.850, lwd=1.5, type="b", col=colors[7], pch=19)
lines(x=FBXO47.meth.2$xLoc, y=FBXO47.meth.2$Bcontrol, lwd=1.5, type="b", col=colors[2], pch=19)
lines(x=FBXO47.meth.2$xLoc, y=FBXO47.meth.2$Bcase, lwd=1.5, type="b", col=colors[8], pch=19)

legend("topleft", fill=colors[c(8,2,7,1)], legend=c("Cases 450K", "Controls 450K", "Cases EPIC", "Controls EPIC"), cex=1)

axis(side=2, labels=TRUE)
#mtext(c("*", "???"), side = 3, line=0, at=c(157697058-adj, 157454062-adj), cex=1.2)

#SNPs
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='SNPs',xlab='', xlim=c(0,max(FBXO47.meth.2$xLoc)), ylim=c(0,0.1)) 
for(i in 1:nrow(FBXO47_snps.2)){abline(v=FBXO47_snps.2[i,"toPlot"], col="black")}

#450K DMR
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim=c(0,max(FBXO47.meth.2$xLoc)), ylim=c(0,0.5))
rect(range450.adj.zoom$start,0,range450.adj.zoom$end,1, col= "gray", border = FALSE)
text(range450.adj.zoom$start, 0.25, labels = '450K DMR', pos=4)

#EPIC DMR
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim=c(0,max(FBXO47.meth.2$xLoc)), ylim=c(0,0.1))
rect(rangeEPIC.adj.zoom$start,0,rangeEPIC.adj.zoom$end,1, col= "gray",border = FALSE)
text(rangeEPIC.adj.zoom$start, 0.05, labels = 'EPIC DMR', pos=4)

#Zoomed in legend;
plot(1, type="n", xlab="", ylab="", xlim=c(0, max(FBXO47.meth.2$xLoc)), ylim=c(0, 1), frame.plot=FALSE,axes=FALSE)
axis(1, at=seq(0,300,by=50), labels=seq(0,300,by=50)+adj.zoom)
#text(20, 0.5, labels = 'Chr17', pos=2)
mtext("Chr 17",2, adj=0, cex=0.75, las=1)

##Gene Structure
#par(mfrow = c(1,1), mar=c(5,5,5,5))
genemodel.plot.LV(model=FBXO47.toPlot, start=min(FBXO47[,3]), bpstop=max(FBXO47[,2]), orientation="reverse", xaxis=TRUE, gene="FBXO47")
#rect(bothDMRrange[1],-1,bothDMRrange[2],1, col= "#0000FF80", border=FALSE)
rect(bothDMRrange[1],-1,bothDMRrange[2],1, col= "gray", border=FALSE)
#mtext("Chr 17",3, adj=0, cex=0.75)

dev.off()

