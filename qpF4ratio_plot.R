A<-read.table("Documents/LabProject/dog/F3/qpF4ratio_dog_from_ChinaS_Coyote_result.txt")
A$V9=A$V6-A$V7
A$V10=A$V6+A$V7
A=A[-which(A$V6>0.7|A$V6<0),]
A=A[order(A[,6]),]
A=A[-c(1),]
par(mfrow=c(1,1))

Label=c("Taiwan","Papua New Guinea","North China","Siberian Husky","Pekingnese","India Tibetan Mastiff Mix","China Kunming","China Kazakhstan","HXH","Bosnia Tornjak","Boxer","Xolo","Chihuahua","China Mongolian Shepherd","Bulldog","Labrador Retriever","Beagle","Toy Poodle","Flat Coated Retriever","Istrian Shorthaired Hound","Kerry Blue Terrier","Standard Poodle","English Cocker Spaniel","Chinese Crest","Pembroke Welsh Corgi","Scottish Terrier","Great Dane","Mastiff")

pdf("Documents//LabProject/dog/F3/qpF4ratio_dog_from_ChinaS_coyote_label.pdf")
barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main="Proportion of admixture from south China village dogs")
segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
text( par("usr")[1],barcenters+0.8, labels = rev(Label), col=c(rep("blue",14),"black",rep("blue",4),"red","black","black","blue","blue","blue",rep("black",3)),srt = 0, pos = 1, xpd = TRUE,cex=0.6)
axis(1,at=seq(-1,1,0.1))
dev.off()

A<-read.table("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_dog_from_India_Coyote_result.txt")
A$V9=A$V6-A$V7
A$V10=A$V6+A$V7
A=A[-which(A$V6<(-1)|(A$V6>1)),]
A=A[order(A[,6]),]
A=A[-c(25),]
par(mfrow=c(1,1))

pdf("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_dog_from_India_coyote.pdf")
barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main="proportion of ancestry from India dogs")
segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
text( par("usr")[1],barcenters+0.8, labels = A$V3, srt = 0, pos = 1, xpd = TRUE,cex=0.6)
axis(1,at=seq(-1,1,0.1))
dev.off()

A<-read.table("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_IsraeliWolf_from_dog_result.txt")
A$V9=A$V6-A$V7
A$V10=A$V6+A$V7
A=A[order(A[,6]),]

pdf("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_IsraeliWolf_from_dog_coyote.pdf")
barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main="proportion of dog ancestry in IsraeliWolf")
segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
text( par("usr")[1],barcenters+0.8, labels = A$V2, srt = 0, pos = 1, xpd = TRUE,cex=0.6)
axis(1,at=seq(-1,1,0.1))
dev.off()

wolf_list=c('wolfSpain', 'IsraeliWolf', 'ChineseWolf','wolfItaly', 'wolfPortugal', 'wolfIberian', 'wolfIran','CroatianWolf')
ancestry=c()
for(i in 1:length(wolf_list)){
  A<-read.table(paste("Documents//LabProject/dog/F3/qpF4ratio_",wolf_list[i],"_from_dog_result.txt",sep=""))
  A$V9=A$V6-A$V7
  A$V10=A$V6+A$V7
  A=A[order(A[,6]),]
  ancestry=rbind(ancestry,summary(A$V6))
#  pdf(paste("Documents//LabProject/dog/F3/qpF4ratio_",wolf_list[i],"_from_dog_coyote.pdf",sep=""))
  barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main=paste("proportion of dog ancestry in ",wolf_list[i],sep=""))
  segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
  arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
  text( par("usr")[1],barcenters+0.8, labels = A$V2, srt = 0, pos = 1, xpd = TRUE,cex=0.6)
  axis(1,at=seq(-1,1,0.1))
#  dev.off()
}

#rownames(ancestry)=wolf_list
rownames(ancestry)=c("Spain","Israeli","China",'Italy','Portugal','Iberia','Iran','Croatia')
ancestry=ancestry[order(ancestry[,4]),]
pdf("Documents/LabProject/dog/F3/qpF4ratio_dog_ancestry_in_wolf_summary_label2.pdf")
barcenters=barplot(height=ancestry[,4],ylim=c(0,0.25),main="Inferred dog ancestry in each wolf",cex.names=0.6)
segments(barcenters,ancestry[,1],barcenters,ancestry[,6],lwd = 1.5)
arrows(barcenters,ancestry[,1],barcenters,ancestry[,6],lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()

wolf_list=c('wolfSpain', 'IsraeliWolf', 'ChineseWolf','wolfItaly', 'wolfPortugal', 'wolfIberian', 'wolfIran','CroatianWolf')
wolf_list=c('wolfIndia')
for(i in 1:length(wolf_list)){
  A<-read.table(paste("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_dog_from_",wolf_list[i],"_result.txt",sep=""))
  A$V9=A$V6-A$V7
  A$V10=A$V6+A$V7
  A=A[order(A[,6]),]
  pdf(paste("Documents/GoogleDrive/LabProject/dog/F3/qpF4ratio_dog_from_",wolf_list[i],"_coyote.pdf",sep=""))
  barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main=paste("proportion of ",wolf_list[i]," ancestry in dog",sep=""))
  segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
  arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
  text( par("usr")[1],barcenters+0.8, labels = A$V3, srt = 0, pos = 1, xpd = TRUE,cex=0.6)
  axis(1,at=seq(-1,1,0.1))
  dev.off()
}


A<-read.table("Documents/LabProject/dog/F4ratio/qpF4ratio_dog_from_ChinaS_result.txt")
A$V9=A$V6-A$V7
A$V10=A$V6+A$V7
A=A[-which(A$V6>0.7|A$V6<0),]
A=A[order(A[,6]),]
A=A[-c(1),]
par(mfrow=c(1,1))

Label=c("Taiwan","Papua New Guinea","North China","Siberian Husky","Pekingnese","HXH","China Kunming","India Tibetan Mastiff Mix","China Kazakhstan","NewGrange","Boxer","Xolo","Flat Coated Retriever","Bulldog","China Mongolian Shepherd","Labrador Retriever","Kerry Blue Terrier","Chihuahua","Pembroke Welsh Corgi","Chinese Crest","Standard Poodle","Bosnia Tornjak","Toy Poodle","Scottish Terrier","English Cocker Spaniel")

pdf("Documents//LabProject/dog/F4ratio/qpF4ratio_dog_from_ChinaS_fox_label.pdf")
barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main="Proportion of admixture from south China village dogs")
segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
text( par("usr")[1],barcenters+0.3, labels = rev(Label), col=c(rep("blue",15),"red","black","blue","black","red","blue","blue",rep("black",3)),srt = 0, pos = 1, xpd = TRUE,cex=0.6)
axis(1,at=seq(-1,1,0.1))
dev.off()

wolf_list=c('wolfSpain', 'IsraeliWolf', 'ChineseWolf','wolfItaly', 'wolfPortugal', 'wolfIberian', 'wolfIran','CroatianWolf')
ancestry=c()
for(i in 1:length(wolf_list)){
  A<-read.table(paste("Documents//LabProject/dog/F4ratio/qpF4ratio_",wolf_list[i],"_from_dog_result.txt",sep=""))
  A$V9=A$V6-A$V7
  A$V10=A$V6+A$V7
  A=A[order(A[,6]),]
  ancestry=rbind(ancestry,summary(A$V6))
    pdf(paste("Documents//LabProject/dog/F4ratio/qpF4ratio_",wolf_list[i],"_from_dog_fox.pdf",sep=""))
  barcenters=barplot(height=A[,6],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main=paste("proportion of dog ancestry in ",wolf_list[i],sep=""))
  segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
  arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
  text( par("usr")[1],barcenters+0.8, labels = A$V2, srt = 0, pos = 1, xpd = TRUE,cex=0.6)
  axis(1,at=seq(-1,1,0.1))
    dev.off()
}

#rownames(ancestry)=wolf_list
rownames(ancestry)=c("Spain","Israeli","China",'Italy','Portugal','Iberia','Iran','Croatia')
ancestry=ancestry[order(ancestry[,4]),]
pdf("Documents/LabProject/dog/F4ratio/qpF4ratio_dog_ancestry_in_wolf_summary_label_fox.pdf")
barcenters=barplot(height=ancestry[,4],ylim=c(0,0.25),main="Inferred dog ancestry in each wolf",cex.names=0.6)
segments(barcenters,ancestry[,1],barcenters,ancestry[,6],lwd = 1.5)
arrows(barcenters,ancestry[,1],barcenters,ancestry[,6],lwd = 1.5, angle = 90,code = 3, length = 0.05)
dev.off()



A<-read.table("Downloads//qpF4ratio_dog_from_ChinaS_result_v2.txt")
A$V9=1-(A$V6+A$V7)
A$V10=1-(A$V6-A$V7)
A$V11=1-A$V6
A=A[-which(A$V11>1|A$V11<0),]
A=A[-which(A$V3=="CTC"|A$V3=="Basenji"|A$V3=="Europe"|A$V3=="Mideast_EG"|A$V3=="SubSaharan"),]
A=A[order(A[,11]),]
A=A[-c(1),]
par(mfrow=c(1,1))

Label=c("Siberian Husky","Taiwan","North China","Papua New Guinea","NGD","HXH","India Tibetan Mastiff Mix","Pekingnese","China Kazakhstan","Xolo","Beagle","ChineseCrest","China Kunming","AfghanHound","Saluki","GreatDane")
#Label=c("Taiwan","Papua New Guinea","North China","Siberian Husky","Pekingnese","HXH","China Kunming","India Tibetan Mastiff Mix","China Kazakhstan","NewGrange","Boxer","Xolo","Flat Coated Retriever","Bulldog","China Mongolian Shepherd","Labrador Retriever","Kerry Blue Terrier","Chihuahua","Pembroke Welsh Corgi","Chinese Crest","Standard Poodle","Bosnia Tornjak","Toy Poodle","Scottish Terrier","English Cocker Spaniel")

pdf("Documents//LabProject/dog/F4ratio/qpF4ratio_dog_from_ChinaS_fox_label_v2.pdf")
barcenters=barplot(height=A[,11],xlim=c(-1,1),width=0.5,horiz=TRUE,axes=FALSE,main="Proportion of admixture from south China village dogs")
segments(A[,9],barcenters,A[,10],barcenters,lwd = 1.5)
arrows(A[,9],barcenters,A[,10],barcenters,lwd = 1.5, angle = 90,code = 3, length = 0.05)
text( par("usr")[1],barcenters+0.3, labels = A$V3, col=c(rep("blue",15),"red","black","blue","black","red","blue","blue",rep("black",3)),srt = 0, pos = 1, xpd = TRUE,cex=0.6)
axis(1,at=seq(-1,1,0.1))
dev.off()