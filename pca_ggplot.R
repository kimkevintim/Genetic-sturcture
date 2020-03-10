# require packages
library(adegenet)
library(ggplot2)


# input txt file
Rpcf_aflp<-read.table("Rpc_Rf_AFLP_1050_663.txt", sep="\t", header=T)
# load to genind objects
Rpcf.gind<-df2genind(X=Rpcf_aflp[,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[,4], pop=Rpcf_aflp[,3], type="PA")
Rpc.gind<-df2genind(X=Rpcf_aflp[61:233,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[61:233,4], pop=Rpcf_aflp[61:233,3], type="PA")
Rf.gind<-df2genind(X=Rpcf_aflp[1:60,5:667], sep=NULL, ncode=1, ind.names=Rpcf_aflp[1:60,4], pop=Rpcf_aflp[1:60,3], type="PA")


# pca for Rpc and Rf
Rpcf.pca <- dudi.pca(Rpcf.gind, cent=FALSE,scale=FALSE,scannf=FALSE,nf=2)
summary(Rpcf.pca)
Rpcf.pca.li <- data.frame(Rpcf.pca$li)
Rpcf.pca.li$Populations<-Rpcf.gind$pop
head(Rpcf.pca.li)
# gplot
Rpcf.ggpca<-ggplot(Rpcf.pca.li)+geom_point(aes(x=Axis1, y=Axis2, color=Populations,
                                               fill=Populations, shape=Populations), size=4)
Rpcf.ggpca<-Rpcf.ggpca+labs(x="PC1 (67.54%)", y="PC2 (4.75%)")+
  theme_bw()+scale_shape_manual(values=c(1:5, 16:24))+
  scale_colour_manual(values=rainbow(14))+
  scale_fill_manual(values=rainbow(14))
Rpcf.ggpca<-Rpcf.ggpca+theme(legend.position=c(0.11, 0.26),panel.grid.major=element_blank(),
                             panel.grid.minor=element_blank(), legend.text=element_text(size=12), legend.title=element_text(size=14))
Rpcf.ggpca<-Rpcf.ggpca+theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# pca out
tiff(file="Rpcf_AFLP_ggPCA.tif",width=5000,height=5000,res=600, compression = "lzw")
Rpcf.ggpca
dev.off()



# pca for only Rpc
Rpc.pca <- dudi.pca(Rpc.gind, cent=FALSE,scale=FALSE,scannf=FALSE,nf=2)
summary(Rpc.pca)
Rpc.pca.li <- data.frame(Rpc.pca$li)
Rpc.pca.li$Populations<-Rpc.gind$pop
head(Rpc.pca.li)
# gplot
Rpc.ggpca<-ggplot(Rpc.pca.li)+geom_point(aes(x=Axis1, y=Axis2, color=Populations,
                                             fill=Populations, shape=Populations), size=4)
Rpc.ggpca<-Rpc.ggpca+labs(x="PC1 (72.76%)", y="PC2 (1.86%)")+
  theme_bw()+scale_shape_manual(values=c(16:24))+
  scale_colour_manual(values=rainbow(14)[6:14])+
  scale_fill_manual(values=rainbow(14)[6:14])
Rpc.ggpca<-Rpc.ggpca+theme(legend.position=c(0.89, 0.18),panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(), legend.text=element_text(size=12), legend.title=element_text(size=14))
Rpc.ggpca<-Rpc.ggpca+theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# pca out
tiff(file="Rpc_AFLP_ggPCA.tif",width=5000,height=5000,res=600, compression = "lzw")
Rpc.ggpca
dev.off()



# pca for only Rf
Rf.pca <- dudi.pca(Rf.gind, cent=FALSE,scale=FALSE,scannf=FALSE,nf=2)
summary(Rf.pca)
Rf.pca.li <- data.frame(Rf.pca$li)
Rf.pca.li$Populations<-Rf.gind$pop
head(Rf.pca.li)
# gplot
Rf.ggpca<-ggplot(Rf.pca.li)+geom_point(aes(x=Axis1, y=Axis2, color=Populations,
                                           fill=Populations, shape=Populations), size=4)
Rf.ggpca<-Rf.ggpca+labs(x="PC1 (64.54%)", y="PC2 (2.69%)")+
  theme_bw()+scale_shape_manual(values=c(1:5))+
  scale_colour_manual(values=rainbow(14)[1:5])+
  scale_fill_manual(values=rainbow(14)[1:5])
Rf.ggpca<-Rf.ggpca+theme(legend.position=c(0.89, 0.87),panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(), legend.text=element_text(size=12), legend.title=element_text(size=14))
Rf.ggpca<-Rf.ggpca+theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
# pca out
tiff(file="Rf_AFLP_ggPCA.tif",width=5000,height=5000,res=600, compression = "lzw")
Rf.ggpca
dev.off()
