library(ape)
E<-read.table("Documents/LabProject/dog/NJ_tree/dist_matrix_canfam3.1_98genomes.txt",header=TRUE)
label<-read.table("Documents//LabProject/dog/Kidd_Freedman_aDNA_merge_label_nRW_detail_for_paper2.txt",sep="\t")
E=as.matrix(E)
colnames(E)=label$V5
colors=rep("red",length(E[1,]))
colors[which(label$V3=="Breed")]="blue"
colors[which(label$V3=="VillageDog")]="black"
colors[which(label$V3=="Wolf")]="orange"
colors[which(label$V3=="Coyote")]="darkgreen"
colors[which(label$V3=="Jackal")]="lightgreen"
colors[which(label$V3=="AncientDog")]="red"

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_label.pdf",10,10)
nj_values = nj(as.dist(E))
write.tree(root(nj_values, "Golden_Jackal"), file = "Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_newick_label.txt", append = TRUE,digits = 10, tree.names = TRUE)
plot(root(nj_values, "Golden_Jackal"), font = 0.7, cex = 0.6, edge.width = 0.7,tip.color=colors)
add.scale.bar(x=0.01, y=1,0.02)
dev.off()

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_bootstrap.pdf",10,10)
for(i in 0:99){
  E<-read.table(paste("Documents/LabProject/dog/NJ_tree/dist_matrix_canfam3.1_98genomes_bootstrap_",i,".txt",sep=""),skip=1,header=FALSE)
  E=as.matrix(E)
  colnames(E)=label$V5
  nj_values = nj(as.dist(E))
  write.tree(root(nj_values, "Golden_Jackal"), file = "Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_bootstrap_newick_label.txt", append = TRUE,digits = 10, tree.names = TRUE)
  plot(root(nj_values, "Golden_Jackal"), font = 0.7, cex = 0.6, edge.width = 0.7,tip.color=colors)
  add.scale.bar(x=0.01, y=1,0.02)
  print(i)
}

dev.off()

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_bootstrap_value_NJ_tree_label.pdf",8,10)
A=read.tree("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_bootstrap_score_label.nw")
color_match=colors[match(A$tip.label,colnames(E))]
plot.phylo(A,show.node.label=TRUE,tip.color=color_match,edge.width = 0.7,font = 0.7, cex = 0.6)
add.scale.bar(x=0.01, y=1,0.02)
dev.off()


## subset

library(ape)
E<-read.table("Documents/LabProject/dog/NJ_tree/dist_matrix_canfam3.1_98genomes.txt",header=TRUE)
label<-read.table("Documents//LabProject/dog/Kidd_Freedman_aDNA_merge_label_nRW_detail_for_paper2.txt",sep="\t")
E=as.matrix(E)
colnames(E)=label$V5

colors=rep("red",length(E[1,]))
colors[which(label$V3=="Breed")]="blue"
colors[which(label$V3=="VillageDog")]="black"
colors[which(label$V3=="Wolf")]="orange"
colors[which(label$V3=="Coyote")]="darkgreen"
colors[which(label$V3=="Jackal")]="lightgreen"
colors[which(label$V3=="AncientDog")]="red"

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_subset2_NJ_tree_label.pdf",10,10)
E_subset=E[c(2,4,7,8,24:29,31:34,37:52,56:68,70,71,74:77,79:81,83:84,87:89,91:98),c(2,4,7,8,24:29,31:34,37:52,56:68,70,71,74:77,79:81,83:84,87:89,91:98)]
colors_subset=colors[c(2,4,7,8,24:29,31:34,37:52,56:68,70,71,74:77,79:81,83:84,87:89,91:98)]
nj_values = nj(as.dist(E_subset))
write.tree(root(nj_values, "Golden_Jackal"), file = "Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_subset2_newick_label.txt", append = TRUE,digits = 10, tree.names = TRUE)
plot(root(nj_values, "Golden_Jackal"), font = 0.7, cex = 0.6, edge.width = 0.7,tip.color=colors_subset)
add.scale.bar(x=0.01, y=1,0.02)
dev.off()

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_subset2_NJ_tree_bootstrap_label.pdf",10,10)
for(i in 0:99){
  E<-read.table(paste("Documents/LabProject/dog/NJ_tree/dist_matrix_canfam3.1_98genomes_bootstrap_",i,".txt",sep=""),skip=1,header=FALSE)
  E=as.matrix(E)
  colnames(E)=label$V5
  E_subset=E[c(2,4,7,8,24:29,31:34,37:52,56:68,70,71,74:77,79:81,83:84,87:89,91:98),c(2,4,7,8,24:29,31:34,37:52,56:68,70,71,74:77,79:81,83:84,87:89,91:98)]
  nj_values = nj(as.dist(E_subset))
  write.tree(root(nj_values, "Golden_Jackal"), file = "Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_subset2_bootstrap_newick_label.txt", append = TRUE,digits = 10, tree.names = TRUE)
  plot(root(nj_values, "Golden_Jackal"), font = 0.7, cex = 0.6, edge.width = 0.7,tip.color=colors_subset)
  add.scale.bar(x=0.01, y=1,0.02)
  print(i)
}

dev.off()

pdf("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_subset2_bootstrap_value_NJ_tree_label.pdf",10,10)
A=read.tree("Documents/LabProject/dog/NJ_tree/Kidd_Freedman_aDNA_wgs_snp_NJ_tree_subset2_bootstrap_score_label.nw")
color_match=colors[match(A$tip.label,colnames(E))]
plot.phylo(A,show.node.label=TRUE,tip.color=color_match,edge.width = 0.7,font = 1, cex = 0.6)
add.scale.bar(x=0.01, y=1,0.02)
dev.off()