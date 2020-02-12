read_Gphocs_calc_mig<- function(file,burn){
  #  mig=read.table(file,sep="\t",header=TRUE)
  mig=read.table(file,sep="\t",header=TRUE,nrows=1)
  header=colnames(mig)[-1]
  print(file)
  print(header)
  mig=read.table(file,sep="\t",header=FALSE,skip=burn)
  mig=mig[,-1]
  mig=mig[,-(dim(mig)[2]-2)]
  colnames(mig)=header
  mig$tau=mig$theta_root*0.5+mig$tau_root
  mig$m_AG=mig$m_A..G*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_GA=mig$m_G..A*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_BG=mig$m_B..G*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_GB=mig$m_G..B*min(mig$tau_GH,mig$tau_AB)*0.1
  if("m_C..F" %in% header){
    mig$m_CF=mig$m_C..F*min(mig$tau_FGH,mig$tau_ABCD)*0.1
    mig$m_FC=mig$m_F..C*min(mig$tau_FGH,mig$tau_ABCD)*0.1
    mig$m_CG=mig$m_C..F*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_GC=mig$m_F..C*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_CH=mig$m_C..F*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_HC=mig$m_F..C*min(mig$tau_GH,mig$tau_ABCD)*0.1 
  }
  if("m_C..D" %in% header){
    mig$m_CD=mig$m_C..D*min(mig$tau_ABD,mig$tau_ABCD)*0.1  #
    mig$m_DC=mig$m_D..C*min(mig$tau_ABD,mig$tau_ABCD)*0.1  #
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  } else if ("m_C..E" %in% header){
    mig$m_CE=mig$m_C..E*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EC=mig$m_E..C*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  } else{
    mig$m_EF=mig$m_E..F*min(mig$tau_ABDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABDEFGH..I*(mig$tau_root-mig$tau_ABDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABDEFGH*(mig$tau_root-mig$tau_ABDEFGH)*0.1  #
  }
  if("m_D..H" %in% header){
    if ("tau_ABCD" %in% header){
      mig$m_DH=mig$m_D..H*min(mig$tau_GH,mig$tau_ABCD)*0.1  #
      mig$m_HD=mig$m_H..D*min(mig$tau_GH,mig$tau_ABCD)*0.1  #
    }
    else{
      mig$m_DH=mig$m_D..H*min(mig$tau_GH,mig$tau_ABD)*0.1  #
      mig$m_HD=mig$m_H..D*min(mig$tau_GH,mig$tau_ABD)*0.1  #   
    }
  }
  mig$m_DG=mig$m_D..G*min(mig$tau_ABD,mig$tau_GH)*0.1  #
  mig$m_GD=mig$m_G..D*min(mig$tau_ABD,mig$tau_GH)*0.1  #
  mig$m_GI=mig$m_G..I*min(mig$tau_GH,mig$tau_root)*0.1  #
  mig$m_IG=mig$m_I..G*min(mig$tau_GH,mig$tau_root)*0.1  #
  vector=apply(mig,2,quantile_95)
  vector=data.frame(vector)
  colnames(vector)=colnames(mig)
  return(vector)
}

file="Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfIndian_ChinaN_Indian_Europe_Box.log"
vector=read_Gphocs_calc_mig(file,300000)
# tau*(1/10^4/1e-08*3)
divergence_time2=data.frame(mean=c(vector$tau_AB[1],vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time2=divergence_time2*(1/10^4/1e-08*3)
write.table(divergence_time2,file="Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfIndian_ChinaN_Indian_Europe_Box_divergence_time.txt",row.names=FALSE,sep="\t")

pdf('Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfIndian_ChinaN_Indian_Europe_Box_divergence_time.pdf')
plot(x=divergence_time$mean,y=seq(1,by=0.2,length=6),axes=F,xlim=c(0,3e04),ylab="",ylim=c(1,10),col=c(rep("orange",3),rep("blue",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",1,'x',10^-8,")")))
segments(divergence_time$lower,seq(1,by=0.2,length=6),divergence_time$upper,seq(1,by=0.2,length=6),lwd = 1.5,col=c(rep("orange",3),rep("blue",2),"darkgreen"))
arrows(divergence_time$lower,seq(1,by=0.2,length=6),divergence_time$upper,seq(1,by=0.2,length=6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c(rep("orange",3),rep("blue",2),"darkgreen"))
axis(1,at=seq(0,3e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


file="Documents/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box.log"
vector=read_Gphocs_calc_mig(file,300000)
# tau*(1/10^4/1e-08*3)
divergence_time=data.frame(mean=c(vector$tau_AB[1],vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/1e-08*3)
vector=data.frame(mean=vector[,1],lower=vector[,2],upper=vector[,3],name=name_tag)
barcenters=barplot(height=vector[,1],ylim=c(0,1),width=0.5,xlab=label[i],col=rainbow(6))
write.table(divergence_time,file="Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box_divergence_time.txt",row.names=FALSE,sep="\t")

pdf('Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box_divergence_time.pdf')
plot(x=divergence_time$mean,y=seq(1,by=0.2,length=6),axes=F,xlim=c(0,3e04),ylab="",ylim=c(1,10),col=c(rep("orange",3),rep("blue",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",1,'x',10^-8,")")))
segments(divergence_time$lower,seq(1,by=0.2,length=6),divergence_time$upper,seq(1,by=0.2,length=6),lwd = 1.5,col=c(rep("orange",3),rep("blue",2),"darkgreen"))
arrows(divergence_time$lower,seq(1,by=0.2,length=6),divergence_time$upper,seq(1,by=0.2,length=6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c(rep("orange",3),rep("blue",2),"darkgreen"))
axis(1,at=seq(0,3e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


# Use new mutation rate:
file="Documents/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box.log"
vector=read_Gphocs_calc_mig(file,300000)
# tau*(1/10^4/1e-08*3)
tau1_lower=1.349798e-05*1e04
tau1_upper=1.234545e-05*1e04
tau1=1.30065e-05*1e04
C=read.table("Documents/LabProject/dog/threeway_sim/sim1_uniform.txt")
divergence_time=data.frame(mean=c(vector$tau_AB[1],quantile_95(C$V1)[1]*1e04,vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],quantile_95(C$V1)[2]*1e04,vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],quantile_95(C$V1)[3]*1e04,vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)
write.table(divergence_time,file="Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box_divergence_time.txt",row.names=FALSE,sep="\t")

pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k/divergence_time_low_mu3.pdf')
par(mar=c(5,5,5,2))
plot(x=divergence_time$mean,y=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,xlim=c(0,5.5e04),ylab="",ylim=c(1,10),col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
arrows(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
axis(1,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


# Use new mutation rate:
file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_Europe_Box.log"
vector=read_Gphocs_calc_mig(file,300000)
C=read.table("Documents/LabProject/dog/threeway_sim/sim1_uniform_DP7.txt")
# tau*(1/10^4/1e-08*3)
tau1_lower=1.349798e-05*1e04
tau1_upper=1.234545e-05*1e04
tau1=1.30065e-05*1e04
C=read.table("Documents/LabProject/dog/threeway_sim/sim1_uniform.txt")
divergence_time=data.frame(mean=c(vector$tau_AB[1],quantile_95(C$V1)[1]*1e04,vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],quantile_95(C$V1)[2]*1e04,vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],quantile_95(C$V1)[3]*1e04,vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)
write.table(divergence_time,file="Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/WolfCroatian_ChinaN_Indian_Europe_Box_divergence_time.txt",row.names=FALSE,sep="\t")

pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/divergence_time_low_mu.pdf')
par(mar=c(5,5,5,2))
plot(x=divergence_time$mean,y=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,xlim=c(0,4.5e04),ylab="",ylim=c(1,10),col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
arrows(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
axis(1,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


# Use Dingo:
file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_Dingo_Indian_Europe_Box.log"
vector_dingo=read_Gphocs_calc_mig(file,300000)
vector_dingo=t(vector_dingo)
colnames(vector_dingo)=c("mean","lower","upper")
write.table(vector_dingo,file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_Dingo_Europe_Box_raw_estimate.txt",row.names=TRUE,sep="\t")

pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/divergence_time_low_mu2.pdf',10,7)
par(mar=c(5,5,5,2))
plot(x=divergence_time$mean,y=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,xlim=c(0,4.3e04),ylab="",ylim=c(1,10),col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
arrows(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("orange","red",rep("orange",2),rep("blue",2),"darkgreen"))
axis(1,at=seq(0,4.3e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


read_Gphocs_calc_mig_v2<- function(file,burn){
  #  mig=read.table(file,sep="\t",header=TRUE)
  mig=read.table(file,sep="\t",header=TRUE,nrows=1)
  header=colnames(mig)[-1]
  print(file)
  print(header)
  mig=read.table(file,sep="\t",header=FALSE,skip=burn)
  mig=mig[,-1]
  mig=mig[,-(dim(mig)[2]-2)]
  colnames(mig)=header
  mig$tau=mig$theta_root*0.5+mig$tau_root
  mig$m_AG=mig$m_A..G*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_GA=mig$m_G..A*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_BG=mig$m_B..G*min(mig$tau_GH,mig$tau_AB)*0.1
  mig$m_GB=mig$m_G..B*min(mig$tau_GH,mig$tau_AB)*0.1
  if("m_C..F" %in% header){
    mig$m_CF=mig$m_C..F*min(mig$tau_FGH,mig$tau_ABCD)*0.1
    mig$m_FC=mig$m_F..C*min(mig$tau_FGH,mig$tau_ABCD)*0.1
    mig$m_CG=mig$m_C..F*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_GC=mig$m_F..C*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_CH=mig$m_C..F*min(mig$tau_GH,mig$tau_ABCD)*0.1
    mig$m_HC=mig$m_F..C*min(mig$tau_GH,mig$tau_ABCD)*0.1 
  }
  if("m_C..D" %in% header){
    mig$m_CD=mig$m_C..D*min(mig$tau_ABD,mig$tau_ABCD)*0.1  #
    mig$m_DC=mig$m_D..C*min(mig$tau_ABD,mig$tau_ABCD)*0.1  #
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  } 
  if ("m_C..E" %in% header){
    mig$m_CE=mig$m_C..E*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EC=mig$m_E..C*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  } 
  if ("m_E..F" %in% header){
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  }
  if("m_D..H" %in% header){
    if ("tau_ABCD" %in% header){
      mig$m_DH=mig$m_D..H*min(mig$tau_GH,mig$tau_ABCD)*0.1  #
      mig$m_HD=mig$m_H..D*min(mig$tau_GH,mig$tau_ABCD)*0.1  #
    }
    else{
      mig$m_DH=mig$m_D..H*min(mig$tau_GH,mig$tau_ABD)*0.1  #
      mig$m_HD=mig$m_H..D*min(mig$tau_GH,mig$tau_ABD)*0.1  #   
    }
  }
  mig$m_DG=mig$m_D..G*min(mig$tau_ABD,mig$tau_GH)*0.1  #
  mig$m_GD=mig$m_G..D*min(mig$tau_ABD,mig$tau_GH)*0.1  #
  mig$m_GI=mig$m_G..I*min(mig$tau_GH,mig$tau_root)*0.1  #
  mig$m_IG=mig$m_I..G*min(mig$tau_GH,mig$tau_root)*0.1  #
  vector=apply(mig,2,quantile_95)
  vector=data.frame(vector)
  colnames(vector)=colnames(mig)
  return(vector)
}

vector=read_Gphocs_calc_mig_v2('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_NewGrange_Europe_Box.log',300000)
file='Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_NewGrange_Europe_Box.log'
divergence_time=data.frame(mean=c(vector$tau_AB[1],vector$tau_ABC[1],vector$tau_ABCD[1],vector$tau_ABCDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABCDEFGH[1]),lower=c(vector$tau_AB[2],vector$tau_ABC[2],vector$tau_ABCD[2],vector$tau_ABCDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABCDEFGH[2]),upper=c(vector$tau_AB[3],vector$tau_ABC[3],vector$tau_ABCD[3],vector$tau_ABCDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABCDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)

vector2=read_Gphocs_calc_mig_v2('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfIndian_ChinaN_Indian_NewGrange_Europe_Box.log',100000)
divergence_time=data.frame(mean=c(vector2$tau_AB[1],vector2$tau_ABC[1],vector2$tau_ABCD[1],vector2$tau_ABCDE[1],vector2$tau_GH[1],vector2$tau_FGH[1],vector2$tau_ABDEFGH[1]),lower=c(vector2$tau_AB[2],vector2$tau_ABC[2],vector2$tau_ABCD[2],vector2$tau_ABCDE[2],vector2$tau_GH[2],vector2$tau_FGH[2],vector2$tau_ABDEFGH[2]),upper=c(vector2$tau_AB[3],vector2$tau_ABC[3],vector2$tau_ABCD[3],vector2$tau_ABCDE[3],vector2$tau_GH[3],vector2$tau_FGH[3],vector2$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)

file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_Europe_Box_full.log"
vector=read_Gphocs_calc_mig(file,300000)
raw_estimate=t(vector)
colnames(raw_estimate)=c("mean","lower","upper")
write.table(raw_estimate,file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_Europe_Box_full_raw_estimate.txt",row.names=TRUE,sep="\t")

divergence_time=data.frame(mean=c(vector$tau_AB[1],vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)

file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_Dingo_Indian_Europe_Box_full.log"
vector2=read_Gphocs_calc_mig(file,300000)
raw_estimate=t(vector2)
colnames(raw_estimate)=c("mean","lower","upper")
write.table(raw_estimate,file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_Dingo_Indian_Europe_Box_full_raw_estimate.txt",row.names=TRUE,sep="\t")

divergence_time2=data.frame(mean=c(vector2$tau_AB[1],vector2$tau_ABD[1],vector2$tau_ABDE[1],vector2$tau_GH[1],vector2$tau_FGH[1],vector2$tau_ABDEFGH[1]),lower=c(vector2$tau_AB[2],vector2$tau_ABD[2],vector2$tau_ABDE[2],vector2$tau_GH[2],vector2$tau_FGH[2],vector2$tau_ABDEFGH[2]),upper=c(vector2$tau_AB[3],vector2$tau_ABD[3],vector2$tau_ABDE[3],vector2$tau_GH[3],vector2$tau_FGH[3],vector2$tau_ABDEFGH[3]))
divergence_time2=divergence_time2*(1/10^4/4e-09*3)


# Plot
file="Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_Europe_Box_full.log"
vector=read_Gphocs_calc_mig(file,300000)
C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_HXH_GJ_fox.txt",skip=1)
divergence_time=data.frame(mean=c(vector$tau_AB[1],quantile_95(C$V1)[1]*1e04,vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],quantile_95(C$V1)[2]*1e04,vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],quantile_95(C$V1)[3]*1e04,vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)

pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/divergence_time_low_mu_full.pdf')
par(mar=c(5,5,5,2))
plot(x=divergence_time$mean,y=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,xlim=c(0,4.5e04),ylab="",ylim=c(1,10),col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
arrows(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
axis(1,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()

# vertical
pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/divergence_time_low_mu_full_vertical.pdf')
par(mar=c(5,5,5,2))
plot(y=divergence_time$mean,x=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,ylim=c(0,4.5e04),xlab="",xlim=c(1,10),col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"),ylab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,lwd = 1.5,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
arrows(c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
axis(2,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()


## 02/20/17
quantile_95<-function(x){
  a=c(mean(x),quantile(x,0.05,names=FALSE),quantile(x,0.95,names=FALSE))
  return(a)
}

# Plot
file="/Volumes/FreeAgent/LabProject/dog/Gphocs/Gphocs_5k_v2/WolfCroatian_ChinaN_Indian_Europe_Box_full.log"
vector=read_Gphocs_calc_mig(file,300000)
C=read.table("Downloads/sim_full_result_HXH_GJ_fox_change_alpha.txt",skip=1)
divergence_time=data.frame(mean=c(vector$tau_AB[1],quantile_95(C$V1)[1]*1e04,vector$tau_ABD[1],vector$tau_ABDE[1],vector$tau_GH[1],vector$tau_FGH[1],vector$tau_ABDEFGH[1]),lower=c(vector$tau_AB[2],quantile_95(C$V1)[2]*1e04,vector$tau_ABD[2],vector$tau_ABDE[2],vector$tau_GH[2],vector$tau_FGH[2],vector$tau_ABDEFGH[2]),upper=c(vector$tau_AB[3],quantile_95(C$V1)[3]*1e04,vector$tau_ABD[3],vector$tau_ABDE[3],vector$tau_GH[3],vector$tau_FGH[3],vector$tau_ABDEFGH[3]))
divergence_time=divergence_time*(1/10^4/4e-09*3)

pdf('Documents/LabProject/dog/Gphocs/Gphocs_5k_v2/divergence_time_low_mu_full.pdf')
par(mar=c(5,5,5,2))
plot(x=divergence_time$mean,y=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,xlim=c(0,4.5e04),ylab="",ylim=c(1,10),col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"),xlab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
arrows(divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("blue","red",rep("blue",2),rep("orange",2),"darkgreen"))
axis(1,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()

# vertical
pdf('Downloads/divergence_time_low_mu_full_vertical.pdf')
par(mar=c(5,5,5,2))
plot(y=divergence_time$mean,x=c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),axes=F,ylim=c(0,4.5e04),xlab="",xlim=c(1,10),col=c("blue","red","darkgreen","purple",rep("orange",3)),ylab=expression(paste("Years (g=3,u=",4,'x',10^-9,")")))
segments(c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,lwd = 1.5,col=c("blue","red","darkgreen","purple",rep("orange",3)))
arrows(c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$lower,c(1.2,1.4,1.6,1.8,1.2,1.4,1.6),divergence_time$upper,lwd = 1.5, angle = 90,code = 3, length = 0.05,col=c("blue","red","darkgreen","purple",rep("orange",3)))
axis(2,at=seq(0,5.5e04,by=1e03), col.axis="black", las=1,cex.axis=0.6)
dev.off()