quantile_95<-function(x){
  return(quantile(x,c(0.05,0.95),names=FALSE))
}

A=read.table("Documents/LabProject/dog/threeway_sim/sim1.txt")
plot(A$V1,A$V2,xlab="tau(HXH-Europe)",ylab="AB/AC ratio")
points(y=1.27,x=1e-05,col="red")
fit=lm(A$V2~A$V1)
Call:
  lm(formula = A$V2 ~ A$V1)

Coefficients:
  (Intercept)         A$V1  
1.638   -28293.535  

tau1=(1.27-1.638)/(-28293.535)
tau1
[1] 1.30065e-05

tau0=7.9e-06
tau2=2.42e-05
tau3=3.5e-05
mu=tau1/7000

B=read.table("Documents/LabProject/dog/threeway_sim/sim1.txt")
plot(A$V1,A$V2,xlab="tau(HXH-Europe)",ylab="AB/AC ratio")

library(gplots)
library(ggplot2)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
pdf("Documents/LabProject/dog/threeway_sim/tau1_vs_alpha.pdf")
colnames(B)[3]="AB/BC ratio"
ggplot(B, aes(V1,V2,color=`AB/BC ratio`)) + geom_point()+scale_color_gradient2(midpoint=1.27,low="blue", mid="white",high="red", space ="Lab" )+xlab("tau(HXH-Europe)")+ylab("AB/AC ratio")
dev.off()

#
tau1_lower=1.349798e-05
tau1_upper=1.234545e-05

C=read.table("Documents/LabProject/dog/threeway_sim/sim1_uniform_DP7.txt")

hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666

pdf("Documents/LabProject/dog/threeway_sim/tau1_DP7_hist.pdf")
hist(C$V1,xlab="divergence time estimates for HXH/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

mut=C$V1/7000*3
pdf("Documents/LabProject/dog/threeway_sim/mutation_DP7_hist.pdf")
hist(mut,xlab="mutation rate",main="",breaks=20)
abline(v=mean(mut),col="red")
abline(v=quantile_95(mut)[2],col="blue")
abline(v=quantile_95(mut)[3],col="blue")
dev.off()

D=read.table("Documents/LabProject/dog/threeway_sim/sim1_uniform_DP7_box.txt")
hist(D$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/tau0_test_hist.pdf")
hist(D$V1,xlab="divergence time estimates for Boxer/Europe",main="",breaks=20)
abline(v=mean(D$V1),col="red")
abline(v=quantile_95(D$V1)[2],col="blue")
abline(v=quantile_95(D$V1)[3],col="blue")
abline(v=vector$tau_AB[1]/1e04,col="red",lty=2)
abline(v=vector$tau_AB[2]/1e04,col="blue",lty=2)
abline(v=vector$tau_AB[3]/1e04,col="blue",lty=2)
dev.off()


# new analysis
D=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_box_GJ.txt",skip=1)
hist(D$V1/4e-09*3)
quantile_95(D$V1)/4e-09*3
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_box_GJ_hist.pdf")
hist(D$V1,xlab="divergence time estimates for Boxer/Europe",main="",breaks=20)
abline(v=mean(D$V1),col="red")
abline(v=quantile_95(D$V1)[2],col="blue")
abline(v=quantile_95(D$V1)[3],col="blue")
abline(v=vector$tau_AB[1]/1e04,col="red",lty=2)
abline(v=vector$tau_AB[2]/1e04,col="blue",lty=2)
abline(v=vector$tau_AB[3]/1e04,col="blue",lty=2)
dev.off()

D=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_box_fox.txt",skip=1)
hist(D$V1/4e-09*3)
quantile_95(D$V1)/4e-09*3
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_box_fox_hist.pdf")
hist(D$V1,xlab="divergence time estimates for Boxer/Europe",main="",breaks=20)
abline(v=mean(D$V1),col="red")
abline(v=quantile_95(D$V1)[2],col="blue")
abline(v=quantile_95(D$V1)[3],col="blue")
abline(v=vector$tau_AB[1]/1e04,col="red",lty=2)
abline(v=vector$tau_AB[2]/1e04,col="blue",lty=2)
abline(v=vector$tau_AB[3]/1e04,col="blue",lty=2)
dev.off()

C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_HXH_GJ.txt")
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_HXH_GJ_hist.pdf")
hist(C$V1,xlab="divergence time estimates for HXH/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_HXH_fox.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_HXH_fox_hist.pdf")
hist(C$V1,xlab="divergence time estimates for HXH/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_GJ.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_GJ_hist.pdf")
hist(C$V1,xlab="divergence time estimates for NGD/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_fox.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_fox_hist.pdf")
hist(C$V1,xlab="divergence time estimates for NGD/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

C=read.table("Downloads/sim_full_result_HXH_GJ_fox_change_alpha.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#6482.671 9719.19 12909.994
pdf("Downloads/sim_full_result_HXH_GJ_fox_change_alpha_hist.pdf")
hist(C$V1,xlab="divergence time estimates for HXH/Europe",main="",breaks=20,cex.axis=0.7)
#axis(side=1,pos=c(5e-06,1e-05,1.2e-05,1.e-05,2e-05))
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[1],col="blue")
abline(v=quantile_95(C$V1)[2],col="blue")
dev.off()

C=read.table("Downloads/sim_full_result_NGD_GJ_fox_change_alpha.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
mean(C$V1)/4e-09*3
#6364.492 9588.367 12592.639
pdf("Downloads/sim_full_result_NGD_GJ_fox_change_alpha_hist.pdf")
hist(C$V1,xlab="divergence time estimates for NGD/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[1],col="blue")
abline(v=quantile_95(C$V1)[2],col="blue")
dev.off()


D=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_box_GJ_fox.txt",skip=1)
hist(D$V1/4e-09*3)
quantile_95(D$V1)/4e-09*3
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_box_GJ_fox_hist.pdf")
hist(D$V1,xlab="divergence time estimates for Boxer/Europe",main="",breaks=20)
abline(v=mean(D$V1),col="red")
abline(v=quantile_95(D$V1)[2],col="blue")
abline(v=quantile_95(D$V1)[3],col="blue")
abline(v=vector$tau_AB[1]/1e04,col="red",lty=2)
abline(v=vector$tau_AB[2]/1e04,col="blue",lty=2)
abline(v=vector$tau_AB[3]/1e04,col="blue",lty=2)
dev.off()

C=read.table("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_GJ_fox.txt",skip=1)
hist(C$V1/4e-09*3)
quantile_95(C$V1)/4e-09*3
#7827.392  4723.704 12686.666
pdf("Documents/LabProject/dog/threeway_sim/sim_full_result_NGD_GJ_fox_hist.pdf")
hist(C$V1,xlab="divergence time estimates for NGD/Europe",main="",breaks=20)
abline(v=mean(C$V1),col="red")
abline(v=quantile_95(C$V1)[2],col="blue")
abline(v=quantile_95(C$V1)[3],col="blue")
dev.off()

C=read.table("Downloads/sim_full_result_HXH_GJ_fox_change_alpha.txt",skip=1)
mut=C$V1/7000*3
mean(mut)
] 5.553823e-09
> quantile_95(mut)
[1] 3.704384e-09 7.377140e-09
pdf("Downloads//mutation_HXH_GJ_fox_change_alpha_hist.pdf")
hist(mut,xlab="mutation rate per generation",main="",breaks=20)
abline(v=mean(mut),col="red")
abline(v=quantile_95(mut)[1],col="blue")
abline(v=quantile_95(mut)[2],col="blue")
dev.off()

