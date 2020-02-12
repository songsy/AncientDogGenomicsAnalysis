A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/Jackel_qp3pop_matrix.txt",header=TRUE)
library(fields)
par(mar=c(4,4.5,3,7))
image(B,col="white",xlab=expression(paste("Split time/Years(x",10^4,")")),ylab=expression(paste("Migration rate(x",10^-5,")")),main="total variation distance between simulated to real TMRCA distribution",axes=FALSE)
axis(1, at=seq(0,1,length=10), labels=c("6","7","8", "9", "10", "11", "12", "13", "14", "15"),lwd = 2)
axis(2, at=seq(0,1,length=11), labels=c("0","2","4", "6", "8", "10", "12", "14","16","18","20"),lwd = 2)
image.plot(B,axes=FALSE,add=TRUE,legend.args=list( text=expression(-log10(variation)),col="black", cex=1, side=3, line=0))

rownames(A)=colnames(A)
B=as.matrix(A)
heatmap(B)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  B[i,i]=20
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=30-B[i,j]
    distance[j,i]=30-B[i,j]
  }
}
distance=as.dist(distance)
a=hclust(distance)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_dendropgram.pdf",12,7)
plot(a)
dev.off()


distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  B[i,i]=20
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=25-B[i,j]
    distance[j,i]=25-B[i,j]
  }
  distance[i,i]=0
}
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_heatmap.pdf",8,7)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()



A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/Jackel_qp3pop_matrix_raw.txt",header=TRUE)

rownames(A)=colnames(A)
B=as.matrix(A)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=2.15-B[i,j]
    distance[j,i]=2.15-B[i,j]
  }
  distance[i,i]=0
}
distance2=as.dist(distance)
distance2=as.dist(distance)
a=hclust(distance2)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_aDNA_raw_f3_dendropgram.pdf",12,7)
plot(a,cex=0.8)
dev.off()

pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_aDNA_raw_f3_heatmap.pdf",8,7)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()

for(i in 1:dim(B)[1]){
  B[i,i]=max(B)
}
heatmap.2(B,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))

## aDNA caller
A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/GoldenJackal_aDNA_qp3pop_dstats_matrix.txt",header=TRUE)

rownames(A)=colnames(A)
B=as.matrix(A)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=2.15-B[i,j]
    distance[j,i]=2.15-B[i,j]
  }
  distance[i,i]=0
}
distance2=as.dist(distance)
a=hclust(distance2)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_aDNA_raw_f3_dendropgram.pdf",12,7)
plot(a)
dev.off()

pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_aDNA_raw_f3_heatmap.pdf",8,7)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()


A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/wolfRed_aDNA_qp3pop_dstats_matrix.txt",header=TRUE)

rownames(A)=colnames(A)
B=as.matrix(A)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=0.9-B[i,j]
    distance[j,i]=0.9-B[i,j]
  }
  distance[i,i]=0
}
distance2=as.dist(distance)
a=hclust(distance2)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_wolfred_aDNA_raw_f3_dendropgram.pdf",12,7)
plot(a)
dev.off()

pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_woflred_aDNA_raw_f3_heatmap.pdf",8,7)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()


## Just plot wolf
A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/GoldenJackal_aDNA_qp3pop_dstats_matrix.txt",header=TRUE)

rownames(A)=colnames(A)
C=A[c(37,38,40:48),1:36]
C=as.matrix(C)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/qp3pop_aDNA_raw_f3_wolf_vs_dog_heatmap.pdf",8,7)
heatmap.2(C,scale="none",Colv=FALSE,key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()

C=A[1:2,3:36]
C=as.matrix(C)
heatmap.2(C,scale="none",Colv=FALSE,key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))

## Dstats
dog_list=c('HXHdog', 'CTCdog', 'ToyPoodle', 'Chihuahua', 'PembrokeWe', 'Saluki', 'Beagle', 'AfghanHound', 'GreatDane', 'ChineseCrest', 'Bulldog', 'BosniaCaucasianOvcharka', 'Pekingnese', 'Boxer', 'BosniaTornjak', 'ScottishTerrier', 'EnglishCockerSpaniel', 'SiberianHusky', 'Dingo', 'Mastiff', 'ChowChow', 'LabradorRetriever', 'Basenji', 'BosniaIstrianShorthairedHound', 'Xolo', 'IndiaTibetanMastiff', 'StandardPoodle', 'FlatCoatedRetriever', 'KerryBlueTerry', 'SharPei','Europe', 'Mideast_LB', 'Oceania_PG', 'India', 'SubSaharan', 'China_N', 'Mideast_EG', 'China_S', 'China_Kazakhstan', 'Vietnam', 'Oceania_Borneo', 'China_MongoliaSheperd', 'Mideast_QA', 'Oceania_TW', 'China_Kunming')
for(i in 1:36){
  A<-read.table(paste("Documents/GoogleDrive/LabProject/dog/Dstats/Jackel_",dog_list[i],"_wolfs_dstats_Zscore_matrix.txt",sep=""),header=TRUE)
  rownames(A)=colnames(A)
  B=as.matrix(A)
  sort_order=order(rowSums(B),decreasing=TRUE)
  B_sorted=B[sort_order,sort_order]
  colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  pdf(paste("Documents/GoogleDrive/LabProject/dog/Dstats/Dstats_",dog_list[i],"_2wolfs.pdf",sep=""))
  heatmap.2(B_sorted,scale="none",Rowv=FALSE,Colv=FALSE,key=TRUE,trace="none",col=my_palette,density.info="none",symbreaks=T, cexRow=0.6,cexCol=0.6,margins=c(6,5),main=dog_list[i])
  dev.off()
}

library(gplots)
wolf_list=c('wolfYellowstone', 'wolfSpain', 'wolfgreatLakes', 'IsraeliWolf', 'wolfChina', 'wolfIndia', 'wolfItaly', 'wolfPortugal', 'wolfIberian', 'wolfIran', 'wolfMexico', 'CroatianWolf')
for(i in 1:12){
  A<-read.table(paste("Documents/GoogleDrive/LabProject/dog/Dstats/Jackel_",wolf_list[i],"_dogs_dstats_Zscore_matrix.txt",sep=""),header=TRUE)
  rownames(A)=colnames(A)
  B=as.matrix(A)
  sort_order=order(rowSums(B),decreasing=TRUE)
  B_sorted=B[sort_order,sort_order]
  colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  #pdf(paste("Documents/GoogleDrive/LabProject/dog/Dstats/Dstats_",dog_list[i],"_2wolfs.pdf",sep=""))
  pdf(paste("Documents/GoogleDrive/LabProject/dog/Dstats/Dstats_",wolf_list[i],"_2dogs.pdf",sep=""))
  heatmap.2(B_sorted,scale="none",Rowv=FALSE,Colv=FALSE,key=TRUE,trace="none",col=my_palette,density.info="none",symbreaks=T, cexRow=0.6,cexCol=0.6,margins=c(6,5),main=wolf_list[i])
  dev.off()
  
}



## qp3pop stats on Shannon dataset, aDNAcaller, results5
A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_GoldenJackal_qp3pop_matrix_raw.txt",header=TRUE)

rownames(A)=colnames(A)
B=as.matrix(A)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=2.8-B[i,j]
    distance[j,i]=2.8-B[i,j]
  }
  distance[i,i]=0
}
distance2=as.dist(distance)
a=hclust(distance2)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_raw_f3_dendropgram.pdf",12,7)
par(cex=0.5, mar=c(5, 8, 4, 1))
plot(a)
dev.off()

pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_raw_f3_heatmap.pdf",12,10)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.5,cexCol=0.4,margins=c(6,5))
dev.off()


A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_wolfRed_qp3pop_matrix_raw.txt",header=TRUE)

rownames(A)=colnames(A)
B=as.matrix(A)
distance=as.matrix(A)
for(i in 1:dim(B)[1]){
  distance[i,i]=0
  for(j in 1:dim(B)[1]){
    distance[i,j]=1.3-B[i,j]
    distance[j,i]=1.3-B[i,j]
  }
  distance[i,i]=0
}
distance2=as.dist(distance)
a=hclust(distance2)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_wolfred_aDNA_raw_f3_dendropgram.pdf",12,7)
plot(a)
dev.off()

pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_woflred_aDNA_raw_f3_heatmap.pdf",8,7)
#par(mar=c(4,5,5,4))
#heatmap.3(distance, scale="column",key=TRUE)
heatmap.2(distance,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()


## Just plot wolf
A<-read.table("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_GoldenJackal_qp3pop_matrix_raw.txt",header=TRUE)

rownames(A)=colnames(A)
C=A[c(179,181:182,185:187,190),1:177]
C=as.matrix(C)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_aDNA_raw_f3_wolf_subset_vs_dog_heatmap.pdf",12,8)
par(mar=c(5, 8, 4, 1))
heatmap.2(C,scale="none",Rowv=FALSE,Colv=FALSE,key=TRUE,trace="none",cexRow=0.6,cexCol=0.5,margins=c(6,5))
dev.off()


C=A[1:2,3:177]
C=as.matrix(C)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_aDNA_raw_f3_dog_vs_ancientdog_heatmap.pdf",10,12)
par(mar=c(5, 8, 4, 1))
heatmap.2(C,scale="none",key=TRUE,trace="none",cexRow=0.6,cexCol=0.4,margins=c(6,5))
dev.off()


C=A[1:2,165:177]
C=as.matrix(C)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_aDNA_raw_f3_village_vs_ancientdog_heatmap.pdf",12,8)
par(mar=c(5, 8, 4, 1))
heatmap.2(C,scale="none",key=TRUE,Rowv=FALSE,Colv=FALSE,trace="none",cexRow=0.6,cexCol=0.6,margins=c(6,5))
dev.off()

C=A[c(1:2,165:177),3:164]
C=as.matrix(C)
pdf("Documents/GoogleDrive/LabProject/dog/Dstats/Shannon_qp3pop_aDNA_raw_f3_dog_vs_village_dog_heatmap.pdf",10,12)
par(mar=c(5, 8, 4, 1))
heatmap.2(C,scale="none",key=TRUE,Rowv=FALSE,Colv=FALSE,trace="none",cexRow=0.6,cexCol=0.4,margins=c(6,5))
dev.off()
