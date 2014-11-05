setwd("/Users/beissinger/Documents/PolygenicAdaptation/AllPhenos")
source ( "Scripts/CreateTraitFile.R")
source ( "Scripts/functions_modified.R")

######################### DO ANALYSIS FOR GWAS P-VALUES < 1e-4 #########################
Out <- list()
for(i in 1:22){
traitFile = paste("Trait_Data/",dir("Trait_Data")[i+22],sep="")
freqsFile = paste("Trait_Data/",dir("Trait_Data")[i],sep="")
pathDir = strsplit(traitFile,split="[.]")[[1]][3]

Out[[i]] <- PolygenicAdaptationFunction (	gwas.data.file = traitFile,
				freqs.file = freqsFile ,
				env.var.data.files = list ( "EnvVar/Environment.txt") ,
							match.pop.file = "Genome_Data/match.pop.file.txt" ,
							full.dataset.file = "Genome_Data/full.dataset.file.txt" ,
							path = paste("OUTPUT/",pathDir,sep="") ,
							match.categories = c ( "FRQ") ,
							match.bins = list ( seq ( 0 , 1 , 0.02 ) ) ,
							cov.SNPs.per.cycle = 5000 ,
							cov.cycles = 1 ,
							null.phenos.per.cycle = 1000 ,
							null.cycles = 10 ,
							load.cov.mat = T ,
							sim.null = F ,
							check.allele.orientation = F
							)
names(Out)[i] <- pathDir
print("------------------------------------------------------------------")
print(Out[[i]])
}


### Make table of asymptotic Qx p-values
QxPvalues <- matrix(NA,ncol=6,nrow=22)
QxPvalues <- data.frame(QxPvalues)
names(QxPvalues) <- c("Phenotype","Overdispersion-P","MH-P","ML-P","SAH-P","SAL-P")
for(i in 1:22){
QxPvalues[i,1] <- names(Out)[i]
QxPvalues[i,2] <- Out[[i]]$asymptotic.p.vals$Qx
QxPvalues[i,3:6] <- Out[[i]]$asymptotic.p.vals$ind.Z
}


### Plot results
pdf("Summary_Plot_10-4.pdf",height=8.5,width=11)
par(mar=c(7,4,4,2))
plot(-log10(QxPvalues[,2]),xaxt="n",xlab="",ylab=expression("-log"[10]*" p-value"),pch=19,ylim=c(0,max(-log10(QxPvalues[,2:6]))),cex=2,main="Test for polygenic adaptation")
axis(1,at=1:22,labels=F)
text(x=1:22,y=-0.3,labels=QxPvalues[,1],srt=45,adj=1,xpd=T,cex=0.75)
abline(h=-log10(0.05),lty=2)
points(-log10(QxPvalues[,3]),col="red",pch=16)
points(-log10(QxPvalues[,4]),col="orange",pch=16)
points(-log10(QxPvalues[,5]),col="blue",pch=16)
points(-log10(QxPvalues[,6]),col="purple",pch=16)

legend("topright","(x,y)",c("Overdispersion", "MH","ML","SAH","SAL"),pch=19,pt.cex=c(2,1,1,1,1),col=c("black","red","orange","blue","purple"))

dev.off()













######################### DO ANALYSIS FOR LOWEST 100 GWAS P-VALUES #########################
Out100 <- list()
for(i in 1:22){
traitFile = paste("Trait_Data100/",dir("Trait_Data")[i+22],sep="")
freqsFile = paste("Trait_Data100/",dir("Trait_Data")[i],sep="")
pathDir = strsplit(traitFile,split="[.]")[[1]][3]

Out100[[i]] <- PolygenicAdaptationFunction (	gwas.data.file = traitFile,
				freqs.file = freqsFile ,
				env.var.data.files = list ( "EnvVar/Environment.txt") ,
							match.pop.file = "Genome_Data/match.pop.file.txt" ,
							full.dataset.file = "Genome_Data/full.dataset.file.txt" ,
							path = paste("OUTPUT100/",pathDir,sep="") ,
							match.categories = c ( "FRQ") ,
							match.bins = list ( seq ( 0 , 1 , 0.02 ) ) ,
							cov.SNPs.per.cycle = 5000 ,
							cov.cycles = 1 ,
							null.phenos.per.cycle = 1000 ,
							null.cycles = 10 ,
							load.cov.mat = T ,
							sim.null = F ,
							check.allele.orientation = F
							)
names(Out100)[i] <- pathDir
print("------------------------------------------------------------------")
print(Out100[[i]])
}


### Make table of asymptotic Qx p-values
QxPvalues100 <- matrix(NA,ncol=6,nrow=22)
QxPvalues100 <- data.frame(QxPvalues100)
names(QxPvalues100) <- c("Phenotype","Overdispersion-P","MH-P","ML-P","SAH-P","SAL-P")
for(i in 1:22){
QxPvalues100[i,1] <- names(Out100)[i]
QxPvalues100[i,2] <- Out100[[i]]$asymptotic.p.vals$Qx
QxPvalues100[i,3:6] <- Out100[[i]]$asymptotic.p.vals$ind.Z
}


### Plot results
pdf("Summary_Plot_100.pdf",height=8.5,width=11)
par(mar=c(7,4,4,2))
plot(-log10(QxPvalues100[,2]),xaxt="n",xlab="",ylab=expression("-log"[10]*" p-value"),pch=19,ylim=c(0,max(-log10(QxPvalues[,2:6]))),cex=2,main="Test for polygenic adaptation")
axis(1,at=1:22,labels=F)
text(x=1:22,y=-0.3,labels=QxPvalues[,1],srt=45,adj=1,xpd=T,cex=0.75)
abline(h=-log10(0.05),lty=2)
points(-log10(QxPvalues100[,3]),col="red",pch=16)
points(-log10(QxPvalues100[,4]),col="orange",pch=16)
points(-log10(QxPvalues100[,5]),col="blue",pch=16)
points(-log10(QxPvalues100[,6]),col="purple",pch=16)

legend("topright","(x,y)",c("Overdispersion", "MH","ML","SAH","SAL"),pch=19,pt.cex=c(2,1,1,1,1),col=c("black","red","orange","blue","purple"))

dev.off()



Warning messages:
1: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
2: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
3: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
4: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
5: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
6: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
7: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
8: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
9: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
10: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
11: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
12: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
13: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
14: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
15: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
16: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
17: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
18: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
19: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
20: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
21: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
22: In mapply(rep, x = split.sampled.cov.data, times = sampled.SNPs.count) :
  longer argument not a multiple of length of shorter
