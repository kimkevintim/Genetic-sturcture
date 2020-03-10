# require packages
library(adegenet)
library(ggplot2)


# input txt file
Rpcf_aflp<-read.table("Rpc_Rf_AFLP_1050_663.txt", sep="\t", header=T)
# load to genind objects
Rpcf.gind<-df2genind(X=Rpcf_aflp[,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[,4], pop=Rpcf_aflp[,3], type="PA")
Rpc.gind<-df2genind(X=Rpcf_aflp[61:232,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[61:232,4], pop=Rpcf_aflp[61:232,3], type="PA")
Rf.gind<-df2genind(X=Rpcf_aflp[1:60,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[1:60,4], pop=Rpcf_aflp[1:60,3], type="PA")


# DAPC for Rpc and Rf load 90% PCs
Rpcf.pca <- dudi.pca(Rpcf.gind, cent=FALSE,scale=FALSE,scannf=FALSE,nf=2)
sum(Rpcf.pca$eig[1:66])/sum(Rpcf.pca$eig)
Rpcf.dapc<-dapc(Rpcf.gind, Rpcf.gind$pop, var.contrib=TRUE, scale=FALSE, n.pca=66, n.da=2)
Rpcf.dapc.li <- data.frame(Rpcf.dapc$ind.coord)
Rpcf.dapc.li$Populations<-Rpcf.gind$pop
head(Rpcf.dapc.li)
Rpcf.dapc$eig[1:2]/sum(Rpcf.dapc$eig)
# ggplot
Rpcf.ggdapc<-ggplot(Rpcf.dapc.li)+geom_point(aes(x=LD1, y=LD2, color=Populations,
                                                 fill=Populations, shape=Populations), size=4)
Rpcf.ggdapc<-Rpcf.ggdapc+labs(x="LD1 (64.18%)", y="LD2 (10.20%)")+
  theme_bw()+scale_shape_manual(values=c(1:5, 16:24))+
  scale_colour_manual(values=rainbow(14))+
  scale_fill_manual(values=rainbow(14))
Rpcf.ggdapc<-Rpcf.ggdapc+theme(legend.position=c(0.5, 0.74),panel.grid.major=element_blank(),
                               panel.grid.minor=element_blank(), legend.text=element_text(size=12), legend.title=element_text(size=14))
Rpcf.ggdapc<-Rpcf.ggdapc+theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# dapc out
tiff(file="Rpcf_AFLP_1050_ggDAPC.tif",width=5000,height=5000,res=600, compression = "lzw")
Rpcf.ggdapc
dev.off()



# DAPC for only Rpc load 90% PCs
Rpc.pca <- dudi.pca(Rpc.gind, cent=FALSE,scale=FALSE,scannf=FALSE,nf=2)
sum(Rpc.pca$eig[1:48])/sum(Rpc.pca$eig)
Rpc.dapc<-dapc(Rpc.gind, Rpc.gind$pop, var.contrib=TRUE, scale=FALSE, n.pca=48, n.da=2)
Rpc.dapc.li <- data.frame(Rpc.dapc$ind.coord)
Rpc.dapc.li$Populations<-Rpc.gind$pop
head(Rpc.dapc.li)
Rpc.dapc$eig[1:2]/sum(Rpc.dapc$eig)
# ggplot
Rpc.ggdapc<-ggplot(Rpc.dapc.li)+geom_point(aes(x=LD1, y=LD2, color=Populations,
                                               fill=Populations, shape=Populations), size=4)
Rpc.ggdapc<-Rpc.ggdapc+labs(x="LD1 (38.78%)", y="LD2 (22.54%)")+
  theme_bw()+scale_shape_manual(values=c(16:24))+
  scale_colour_manual(values=rainbow(14)[6:14])+
  scale_fill_manual(values=rainbow(14)[6:14])
Rpc.ggdapc<-Rpc.ggdapc+theme(legend.position=c(0.89, 0.17),panel.grid.major=element_blank(),
                             panel.grid.minor=element_blank(), legend.text=element_text(size=12), legend.title=element_text(size=14))
Rpc.ggdapc<-Rpc.ggdapc+theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
Rpc.ggpca<-Rpc.ggpca+theme(axis.text=element_text(size=6), axis.title=element_text(size=8))
# dapc out
tiff(file="Rpc_AFLP_1050_ggDAPC.tif",width=5000,height=5000,res=600, compression = "lzw")
Rpc.ggdapc
dev.off()



