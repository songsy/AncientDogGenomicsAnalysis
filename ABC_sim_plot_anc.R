####R plotting

Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/aDNA_Model1_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

mydata.pca <- prcomp(Tbl2, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl2)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/Model1_combo_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,1][10102:20101], scores[,2][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][10002:20001], scores[,2][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,1][2:10001], scores[,2][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,3][10102:20101], scores[,4][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,3][102:10101], scores[,4][102:10101],col="black",pch=20,cex=0.2)
points(scores[,3][2:10001], scores[,4][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][2:10001], scores[,6][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,7][10102:20101], scores[,8][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,7][102:10101], scores[,8][102:10101],col="black",pch=20,cex=0.2)
points(scores[,7][2:10001], scores[,8][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,9][10102:20101], scores[,10][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,9][102:10101], scores[,10][102:10101],col="black",pch=20,cex=0.2)
points(scores[,9][2:10001], scores[,10][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/aDNA_Model1_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][1:20001],main=param[j])
#  abline(v = max(Tbl2[,j][10002:20001]), col = "blue")
#  abline(v = min(Tbl2[,j][10002:20001]), col = "blue")
  abline(v = max(Tbl2[,j][2:20001]), col = "green")
  abline(v = min(Tbl2[,j][2:20001]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  if (Tbl2[,j][1]-min(Tbl2[,j][2:20001])<0| max(Tbl2[,j][2:20001])-Tbl2[,j][1]<0){
    print(c(param[j],Tbl2[,j][1],min(Tbl2[,j][2:20001]),max(Tbl2[,j][2:20001])))
  }
}
dev.off()


A=matrix(c(0,1,1,0,0,1,0,1,1,1,1,1,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0),5,5)
j=80
j=84
j=19

####R plotting

Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

mydata.pca <- prcomp(Tbl2, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl2)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,1][30002:40001], scores[,2][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,1][20002:30001], scores[,2][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,1][10002:20001], scores[,2][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,1][2:10001], scores[,2][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,3][30002:40001], scores[,4][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,3][20002:30001], scores[,4][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,3][10002:20001], scores[,4][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,3][2:10001], scores[,4][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][30002:40001], scores[,6][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,5][20002:30001], scores[,6][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,5][10002:20001], scores[,6][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,5][2:10001], scores[,6][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,7][30002:40001], scores[,8][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,7][20002:30001], scores[,8][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,7][10002:20001], scores[,8][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,7][2:10001], scores[,8][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,9][30002:40001], scores[,10][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,9][20002:30001], scores[,10][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,9][10002:20001], scores[,10][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,9][2:10001], scores[,10][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][1:40001],main=param[j])
  abline(v = max(Tbl2[,j][2:10001]), col = "blue")
  abline(v = min(Tbl2[,j][2:10001]), col = "blue")
  abline(v = max(Tbl2[,j][10002:20001]), col = "black")
  abline(v = min(Tbl2[,j][10002:20001]), col = "black")
  abline(v = max(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = min(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = max(Tbl2[,j][30002:40001]), col = "green")
  abline(v = min(Tbl2[,j][30002:40001]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  if (Tbl2[,j][1]-min(Tbl2[,j][2:40001])<0| max(Tbl2[,j][2:40001])-Tbl2[,j][1]<0){
    print(c(param[j],Tbl2[,j][1],min(Tbl2[,j][2:1018]),max(Tbl2[,j][2:1018])))
  }
}
dev.off()


## HXH
Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

mydata.pca <- prcomp(Tbl2, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl2)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,1][30002:40001], scores[,2][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,1][20002:30001], scores[,2][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,1][10002:20001], scores[,2][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,1][2:10001], scores[,2][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,3][30002:40001], scores[,4][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,3][20002:30001], scores[,4][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,3][10002:20001], scores[,4][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,3][2:10001], scores[,4][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][30002:40001], scores[,6][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,5][20002:30001], scores[,6][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,5][10002:20001], scores[,6][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,5][2:10001], scores[,6][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,7][30002:40001], scores[,8][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,7][20002:30001], scores[,8][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,7][10002:20001], scores[,8][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,7][2:10001], scores[,8][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,9][30002:40001], scores[,10][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,9][20002:30001], scores[,10][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,9][10002:20001], scores[,10][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,9][2:10001], scores[,10][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][1:40001],main=param[j])
  abline(v = max(Tbl2[,j][2:10001]), col = "blue")
  abline(v = min(Tbl2[,j][2:10001]), col = "blue")
  abline(v = max(Tbl2[,j][10002:20001]), col = "black")
  abline(v = min(Tbl2[,j][10002:20001]), col = "black")
  abline(v = max(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = min(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = max(Tbl2[,j][30002:40001]), col = "green")
  abline(v = min(Tbl2[,j][30002:40001]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  #  if (Tbl2[,j][1]-min(Tbl2[,j][2:1018])<0| max(Tbl2[,j][2:1018])-Tbl2[,j][1]<0){
  #    print(c(param[j],Tbl2[,j][1],min(Tbl2[,j][2:1018]),max(Tbl2[,j][2:1018])))
  #  }
}
dev.off()

##CTC with mig
Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_mig_v2_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

Tbl3=Tbl2[,-c(22,46,47,66,67,81,82,91,92,96,97,121,145,146,165,166,180,181,190,191,195,196)]
mydata.pca <- prcomp(Tbl3, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl3)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_mig_v2_rm_pri7_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,1][30002:40001], scores[,2][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,1][20002:30001], scores[,2][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,1][10002:20001], scores[,2][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,1][2:10001], scores[,2][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,3][30002:40001], scores[,4][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,3][20002:30001], scores[,4][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,3][10002:20001], scores[,4][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,3][2:10001], scores[,4][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][30002:40001], scores[,6][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,5][20002:30001], scores[,6][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,5][10002:20001], scores[,6][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,5][2:10001], scores[,6][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,7][30002:40001], scores[,8][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,7][20002:30001], scores[,8][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,7][10002:20001], scores[,8][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,7][2:10001], scores[,8][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,9][30002:40001], scores[,10][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,9][20002:30001], scores[,10][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,9][10002:20001], scores[,10][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,9][2:10001], scores[,10][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/CTC_model1234_mig_v2_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][1:40001],main=param[j])
  abline(v = max(Tbl2[,j][2:10001]), col = "blue")
  abline(v = min(Tbl2[,j][2:10001]), col = "blue")
  abline(v = max(Tbl2[,j][10002:20001]), col = "black")
  abline(v = min(Tbl2[,j][10002:20001]), col = "black")
  abline(v = max(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = min(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = max(Tbl2[,j][30002:40001]), col = "green")
  abline(v = min(Tbl2[,j][30002:40001]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  if (Tbl2[,j][1]-min(Tbl2[,j][2:40001])<0| max(Tbl2[,j][2:40001])-Tbl2[,j][1]<0){
      print(c(param[j],Tbl2[,j][1],min(Tbl2[,j][2:1018]),max(Tbl2[,j][2:1018])))
    }
}
dev.off()

## HXH with mig
Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_mig_v2_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

mydata.pca <- prcomp(Tbl2, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl2)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_mig_v2_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,1][30002:40001], scores[,2][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,1][20002:30001], scores[,2][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,1][10002:20001], scores[,2][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,1][2:10001], scores[,2][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,3][30002:40001], scores[,4][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,3][20002:30001], scores[,4][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,3][10002:20001], scores[,4][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,3][2:10001], scores[,4][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
#points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
#points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][30002:40001], scores[,6][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,5][20002:30001], scores[,6][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,5][10002:20001], scores[,6][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,5][2:10001], scores[,6][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,7][30002:40001], scores[,8][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,7][20002:30001], scores[,8][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,7][10002:20001], scores[,8][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,7][2:10001], scores[,8][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,9][30002:40001], scores[,10][30002:40001],col="green",pch=20,cex=0.2)
points(scores[,9][20002:30001], scores[,10][20002:30001],col="orange",pch=20,cex=0.2)
points(scores[,9][10002:20001], scores[,10][10002:20001],col="black",pch=20,cex=0.2)
points(scores[,9][2:10001], scores[,10][2:10001],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/HXH_model1234_mig_v2_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][1:40001],main=param[j])
  abline(v = max(Tbl2[,j][2:10001]), col = "blue")
  abline(v = min(Tbl2[,j][2:10001]), col = "blue")
  abline(v = max(Tbl2[,j][10002:20001]), col = "black")
  abline(v = min(Tbl2[,j][10002:20001]), col = "black")
  abline(v = max(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = min(Tbl2[,j][20002:30001]), col = "orange")
  abline(v = max(Tbl2[,j][30002:40001]), col = "green")
  abline(v = min(Tbl2[,j][30002:40001]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  if (Tbl2[,j][1]-min(Tbl2[,j][2:40001])<0| max(Tbl2[,j][2:40001])-Tbl2[,j][1]<0){
    print(c(param[j],Tbl2[,j][1],min(Tbl2[,j][2:1018]),max(Tbl2[,j][2:1018])))
  }
}
dev.off()

CTC with mig
"fix14_m"       "0.46916853627" "0.482086"      "0.581413"     
[1] "fix16_m"        "0.439322863553" "0.463254"       "0.571487"      
[1] "fix74_m"        "0.259433000204" "0.27541"        "0.487049"      
[1] "pri75_m"        "0.455163505337" "0.0958597"      "0.431572"      
[1] "pri31_v"       "1.30220602607" "0.55575"       "1.23061"      
[1] "pri41_v"       "1.96417106448" "1.2329"        "1.73362"      
[1] "pri51_v"       "2.08832323641" "1.38419"       "1.99402"      
[1] "fix16_v"        "0.749955583635" "0.77326"        "1.09036"       
[1] "fix74_v"        "0.459310891716" "0.513933"       "0.964089"      
[1] "fix75_v"        "0.384026812878" "0.386777"       "0.767624"      

# after use mig_v2
[1] "fix16_m"        "0.439322863553" "0.462234"       "0.568971"      
[1] "fix16_v"        "0.749955583635" "0.763769"       "1.10774" 

HXH with mig
"fix14_m"        "0.499669443343" "0.514941"       "0.604786"      
[1] "fix16_m"        "0.470911014148" "0.494314"       "0.606307"      
[1] "pri73_m"        "0.571929128653" "0.11748"        "0.526841"      
[1] "pri74_m"        "0.577680814492" "0.112654"       "0.504826"      
[1] "pri75_m"        "0.550178500595" "0.0953325"      "0.46106"       
[1] "pri76_m"        "0.572458019305" "0.106109"       "0.486976"      
[1] "fix14_v"        "0.873297439529" "0.893928"       "1.23291"       
[1] "fix16_v"        "0.811232392693" "0.847898"       "1.18647"      

# after use mig_v2
[1] "fix16_m"        "0.470911014148" "0.500661"       "0.605249" 


####R plotting

Tbl <-read.table("Documents/GoogleDrive/LabProject/dog/simulation/Model1234_combo_PCA", header=TRUE)

stdev=apply(Tbl,2,sd,na.rm=T) #normalize columns (SNPs)
stdev.check=sum(stdev==0) #want this to be 0
#remove monomorphic sites because stdev == 0 for these
if(stdev.check!=0){
  Tbl2=Tbl[,which(stdev!=0)]  # remove SNPs with 0 MAF
}else{
  Tbl2=Tbl
}

mydata.pca <- prcomp(Tbl2, retx=TRUE, center=TRUE,scale.=TRUE)
sd <- mydata.pca$sdev
loadings <- mydata.pca$rotation
rownames(loadings) <- colnames(Tbl2)
scores <- mydata.pca$x

correlations <- t(loadings)*sd

eigenvalues <- sd^2

#plot(eigenvalues/sum(eigenvalues),ylab="% variance", xlab="principle component",type="b", pch=16)


pdf("Documents/GoogleDrive/LabProject/dog/simulation/Model1234_combo_PCA.pdf")
xlab=paste("PCA 1 (",as.character(round(eigenvalues[1]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 2 (",as.character(round(eigenvalues[2]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,1], scores[,2], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,1][20102:30101], scores[,2][20102:30101],col="pink",pch=20,cex=0.2)
points(scores[,1][10102:20101], scores[,2][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][102:10101], scores[,2][102:10101],col="black",pch=20,cex=0.2)
points(scores[,1][2:101], scores[,2][2:101],col="blue",pch=20,cex=0.2)
points(scores[,1][1],scores[,2][1],col="red",pch=20)


xlab=paste("PCA 3 (",as.character(round(eigenvalues[3]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 4 (",as.character(round(eigenvalues[4]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,3], scores[,4], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,3][10102:20101], scores[,4][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][20102:30101], scores[,2][20102:30101],col="pink",pch=20,cex=0.2)
points(scores[,3][102:10101], scores[,4][102:10101],col="black",pch=20,cex=0.2)
points(scores[,3][2:101], scores[,4][2:101],col="blue",pch=20,cex=0.2)
points(scores[,3][1],scores[,4][1],col="red",pch=20)


xlab=paste("PCA 5 (",as.character(round(eigenvalues[5]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 6 (",as.character(round(eigenvalues[6]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,5], scores[,6], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,5][10102:20101], scores[,6][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][20102:30101], scores[,2][20102:30101],col="pink",pch=20,cex=0.2)
points(scores[,5][102:10101], scores[,6][102:10101],col="black",pch=20,cex=0.2)
points(scores[,5][2:101], scores[,6][2:101],col="blue",pch=20,cex=0.2)
points(scores[,5][1],scores[,6][1],col="red",pch=20)


xlab=paste("PCA 7 (",as.character(round(eigenvalues[7]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 8 (",as.character(round(eigenvalues[8]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,7], scores[,8], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,7][10102:20101], scores[,8][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][20102:30101], scores[,2][20102:30101],col="pink",pch=20,cex=0.2)
points(scores[,7][102:10101], scores[,8][102:10101],col="black",pch=20,cex=0.2)
points(scores[,7][2:101], scores[,8][2:101],col="blue",pch=20,cex=0.2)
points(scores[,7][1],scores[,8][1],col="red",pch=20)

xlab=paste("PCA 9 (",as.character(round(eigenvalues[9]/sum(eigenvalues)*100,2)),"%)", sep = "")
ylab=paste("PCA 10 (",as.character(round(eigenvalues[10]/sum(eigenvalues)*100,2)),"%)", sep = "")
plot(scores[,9], scores[,10], xlab=xlab, ylab=ylab,pch=20,cex=0.2,col="white")
points(scores[,9][10102:20101], scores[,10][10102:20101],col="grey",pch=20,cex=0.2)
points(scores[,1][20102:30101], scores[,2][20102:30101],col="pink",pch=20,cex=0.2)
points(scores[,9][102:10101], scores[,10][102:10101],col="black",pch=20,cex=0.2)
points(scores[,9][2:101], scores[,10][2:101],col="blue",pch=20,cex=0.2)
points(scores[,9][1],scores[,10][1],col="red",pch=20)
dev.off()


param=names(Tbl2)

pdf("Documents/GoogleDrive/LabProject/dog/simulation/Model1234_combo_SS.pdf")
for(j in 1:length(param)){
  
  hist(Tbl2[,j][10102:20101],main=param[j])
  abline(v = max(Tbl2[,j][20102:30101]), col = "pink")
  abline(v = min(Tbl2[,j][20102:30101]), col = "pink")
  abline(v = max(Tbl2[,j][102:10101]), col = "blue")
  abline(v = min(Tbl2[,j][102:10101]), col = "blue")
  abline(v = max(Tbl2[,j][2:101]), col = "green")
  abline(v = min(Tbl2[,j][2:101]), col = "green")
  abline(v = Tbl2[,j][1], col = "red")
  print(c(j,Tbl2[,j][1]-min(Tbl2[,j][20102:30101]),max(Tbl2[,j][20102:30101])-Tbl2[,j][1],Tbl2[,j][1]-min(Tbl2[,j][102:10101]),max(Tbl2[,j][102:10101])-Tbl2[,j][1]))
}
dev.off()

j=80
j=84
j=19
j=15

