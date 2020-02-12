read_Gphocs<- function(file,burn){
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
  vector=apply(mig,2,quantile_95)
  vector=data.frame(vector)
  colnames(vector)=colnames(mig)
  return(vector)
}

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
  }
  else if ("m_C..E" %in% header){
    mig$m_CE=mig$m_C..E*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EC=mig$m_E..C*min(mig$tau_ABCD,mig$tau_ABCDE)*0.1  #
    mig$m_EF=mig$m_E..F*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_FE=mig$m_F..E*min(mig$tau_ABCDE,mig$tau_FGH)*0.1  #
    mig$m_ancI=mig$m_ABCDEFGH..I*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
    mig$m_Ianc=mig$m_I..ABCDEFGH*(mig$tau_root-mig$tau_ABCDEFGH)*0.1  #
  }
  else{
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

quantile_95<-function(x){
  a=c(mean(x),quantile(x,0.05,names=FALSE),quantile(x,0.95,names=FALSE))
  return(a)
}

folder='Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/'
wolf_list=c("WolfCroatian","WolfIndian")
dog_list=c("CTC_","HXH_","")
estimate_vector=list()
name_tag=rep(0,6)
index=0
for(i in 1:2){
  for(j in 1:3){
    file=paste(folder,wolf_list[i],"_ChinaN_",dog_list[j],"Indian_Europe_Box.log",sep="")
    vector=read_Gphocs_calc_mig(file,300000)
    index=index+1
    estimate_vector[[index]]=vector 
    name_tag[index]=paste(dog_list[j],wolf_list[i],sep="")
  }
}

# common header
header=unique(c(colnames(estimate_vector[[1]]),colnames(estimate_vector[[2]]),colnames(estimate_vector[[3]]),colnames(estimate_vector[[4]]),colnames(estimate_vector[[5]]),colnames(estimate_vector[[6]])))

for(i in 1:length(header)){
  print(header[i])
  vector=c()
  for(j in 1:6){
    if(header[i] %in% attributes(estimate_vector[[j]])$names){
      index=match(header[i],attributes(estimate_vector[[1]])$names)
      vector=rbind(vector,estimate_vector[[j]][,index])
    }
    else{
      vector=rbind(vector,rep(0,3))
    }
  }
  vector=data.frame(mean=vector[,1],lower=vector[,2],upper=vector[,3],name=name_tag)
}

# population size
pdf("Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/ChinaS_Ancient_Indian_Europe_Box_population_size_300k.pdf",10,7)
par(mfrow=c(3,6))
label=c("Boxer","V_Europe","Ancient","V_India","V_ChinaS","CHW","ISW","CRW/IDW","GLJ","ancDOG1","ancDOG2","ancDOG3*","ancDOG","ancWLF1","ancWLF","ancDW","root")
for(i in 1:17){
  print(header[i])
  vector=c()
  for(j in 1:6){
    if(header[i] %in% attributes(estimate_vector[[j]])$names){
      index=match(header[i],attributes(estimate_vector[[j]])$names)
      vector=rbind(vector,estimate_vector[[j]][,index])
    }
    else{
      if(j %in% c(3,6)){
        if(header[i]=="theta_ABCDE"){
          index=match("theta_ABDE",attributes(estimate_vector[[j]])$names)
          vector=rbind(vector,estimate_vector[[j]][,index])
        }
        else if(header[i]=="theta_ABCDEFGH"){
          index=match("theta_ABDEFGH",attributes(estimate_vector[[j]])$names)
          vector=rbind(vector,estimate_vector[[j]][,index])
        }
        else{
          vector=rbind(vector,rep(0,3))
        }
      }
      else{
        vector=rbind(vector,rep(0,3)) 
      }
    }
  }
  vector=data.frame(mean=vector[,1],lower=vector[,2],upper=vector[,3],name=name_tag)
  barcenters=barplot(height=vector[,1],ylim=c(0,max(vector[,3])*2),width=0.5,xlab=bquote(theta~.(label[i])),col=rainbow(6))
  segments(barcenters,vector[,2],barcenters,vector[,3],lwd = 1)
  arrows(barcenters,vector[,2],barcenters,vector[,3],lwd = 1, angle = 90,code = 3, length = 0.05)
  a=ceiling(max(vector[,3])*1.2*(1/10^4/4/1e-08)/1e03)
  if(a<2){
    increase=0.2
  }
  else if(a<=10){
    increase=1
  }
  else if(a<=20){
    increase=2
  }
  else if(a<=100){
    increase=5
  }
  else{
    increase=100
  }
  axis(4, at = seq(0,a,increase)*1e03/(1/10^4/4/1e-08), label=seq(0,a,increase),col="blue")  # in 
}
par(xpd=NA)
#legend(locator(1), legend=as.numeric(levels(factor(mtcars$cyl))), pch=19, col= as.numeric(levels(factor(mtcars$cyl))) )
legend(x=6, y=20, legend=name_tag,fill= rainbow(6),pch=15,bty="n",cex=0.7)
dev.off()

#ggplot(vector, aes(name_tag, mean,fill=name_tag)) +geom_bar(stat = "identity")+geom_errorbar(aes(ymax=upper,ymin=lower),width = 0.25)+labs(colour="",x=expression(paste(theta,"Boxer")),y="")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(legend.position = "none") 
# divergence time
pdf("Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/ChinaS_Ancient_Indian_Europe_Box_divergence_time_300k.pdf",10,7)
par(mfrow=c(2,5))
label=c("Boxer","V_Europe","Ancient","V_India","V_ChinaS","CHW","ISW","CRW/IDW","GLJ","ancDOG1","ancDOG2","ancDOG3*","ancDOG","ancWLF1","ancWLF2","ancWLF","root","DOG1(BOX,Europe)","DOG2(DOG1,India)","DOG3*(DOG2,Ancient)","DOG(DOG3,China_S)","WLF1(ISW,CRW/IDW)","WLF(WOLF1,CHW)","DW(DOG,WLF)","root")

for(i in 18:25){
  print(header[i])
  vector=c()
  for(j in 1:6){
    if(header[i] %in% attributes(estimate_vector[[j]])$names){
      index=match(header[i],attributes(estimate_vector[[j]])$names)
      vector=rbind(vector,estimate_vector[[j]][,index])
    }
    else{
      if(j %in% c(3,6)){
        if(header[i]=="tau_ABCDE"){
          index=match("tau_ABDE",attributes(estimate_vector[[j]])$names)
          vector=rbind(vector,estimate_vector[[j]][,index])
        }
        else if(header[i]=="tau_ABCDEFGH"){
          index=match("tau_ABDEFGH",attributes(estimate_vector[[j]])$names)
          vector=rbind(vector,estimate_vector[[j]][,index])
        }
        else{
          vector=rbind(vector,rep(0,3))
        }
      }
      else{
        vector=rbind(vector,rep(0,3)) 
      }
    }
  }
  vector=data.frame(mean=vector[,1],lower=vector[,2],upper=vector[,3],name=name_tag)
  barcenters=barplot(height=vector[,1],ylim=c(0,max(vector[,3])*2),width=0.5,xlab=bquote(tau~.(label[i])),col=rainbow(6))
  segments(barcenters,vector[,2],barcenters,vector[,3],lwd = 1.5)
  arrows(barcenters,vector[,2],barcenters,vector[,3],lwd = 1.5, angle = 90,code = 3, length = 0.05)
  a=ceiling(max(vector[,3])*1.2*(1/10^4/1e-08*3)/1e03)
  if(a<2){
    increase=0.2
  }
  else if(a<=10){
    increase=1
  }
  else if(a<=20){
    increase=2
  }
  else if(a<=100){
    increase=5
  }
  else{
    increase=100
  }
  axis(4, at = seq(0,a,increase)*1e03/(1/10^4/1e-08*3), label=seq(0,a,increase),col="blue")  # in 
}
par(xpd=NA)
#legend(locator(1), legend=as.numeric(levels(factor(mtcars$cyl))), pch=19, col= as.numeric(levels(factor(mtcars$cyl))) )
legend(x=6, y=20, legend=name_tag,fill= rainbow(6),pch=15,bty="n",cex=0.7)
dev.off()


# migration rate
folder='Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/'
wolf_list=c("WolfCroatian","WolfIndian")
dog_list=c("CTC_","HXH_","")
estimate_vector=list()
name_tag=rep(0,6)
index=0
for(i in 1:2){
  for(j in 1:3){
    file=paste(folder,wolf_list[i],"_ChinaN_",dog_list[j],"Indian_Europe_Box.log",sep="")
    vector=read_Gphocs_calc_mig(file,400000)
    index=index+1
    estimate_vector[[index]]=vector 
    name_tag[index]=paste(dog_list[j],wolf_list[i],sep="")
  }
}


pdf("Documents/GoogleDrive/LabProject/dog/Gphocs/Gphocs_5k/ChinaS_Ancient_Indian_Europe_Box_migration_rate_300k.pdf",10,7)
par(mfrow=c(2,6))
migration_header=c("m_AG","m_GA","m_BG","m_GB","m_CF","m_FC","m_CG","m_GC","m_CH","m_HC","m_CD","m_DC","m_CE","m_EC","m_DG","m_GD","m_DH","m_HD","m_EF","m_FE","m_GI","m_IG","m_ancI","m_Ianc")
label=c("BOX->\nISW","ISW->\nBOX","Europe->\nISW","ISW->\nEurope","Ancient->\nCHW","CHW->\nAncient","Ancient->\nISW","ISW->\nAncient","Ancient->\nCRW/IDW","CRW/IDW->\nAncient","Ancient->\nIndia","India->\nAncient","Ancient->\nChinaS","ChinaS->\nAncient","India->\nISW","ISW->\nIndia","India->\nIDW","IDW->\nIndia","ChinaS->\nCHW","CHW->\nChinaS","ISW->\nGLJ","GLJ->\nISW","ancDW->\nGLJ","GLJ->\nancDW")
for(i in 1:24){
  print(migration_header[i])
  vector=c()
  for(j in 1:6){
    if(migration_header[i] %in% attributes(estimate_vector[[j]])$names){
      index=match(migration_header[i],attributes(estimate_vector[[j]])$names)
      print(c(j,index))
      vector=rbind(vector,estimate_vector[[j]][,index])
    }
    else{
        vector=rbind(vector,rep(0,3)) 
    }
  }
  vector=data.frame(mean=vector[,1],lower=vector[,2],upper=vector[,3],name=name_tag)
  barcenters=barplot(height=vector[,1],ylim=c(0,1),width=0.5,xlab=label[i],col=rainbow(6))
  segments(barcenters,vector[,2],barcenters,vector[,3],lwd = 1.5)
  arrows(barcenters,vector[,2],barcenters,vector[,3],lwd = 1.5, angle = 90,code = 3, length = 0.05)
  if(i==12){
    par(xpd=NA)
    #legend(locator(1), legend=as.numeric(levels(factor(mtcars$cyl))), pch=19, col= as.numeric(levels(factor(mtcars$cyl))) )
    legend(x=0.8, y=1, legend=name_tag,fill= rainbow(6),pch=15,bty="n",cex=0.7)
  }
}
par(xpd=NA)
#legend(locator(1), legend=as.numeric(levels(factor(mtcars$cyl))), pch=19, col= as.numeric(levels(factor(mtcars$cyl))) )
legend(x=0.8, y=1, legend=name_tag,fill= rainbow(6),pch=15,bty="n",cex=0.7)
dev.off()
