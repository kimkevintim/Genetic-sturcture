# require packages
library(adegenet)
library(poppr)
library(ape)
library(pegas)

# input txt file
#Rpcf_aflp<-read.table("Rpc_Rf_AFLP_1550_535.txt", sep="\t", header=T)
#Rpcf_aflp<-read.table("Rpc_Rf_AFLP_CJJ_1050_337.txt", sep="\t", header=T)
Rpcf_aflp<-read.table("Rpc_Rf_AFLP_CJJ_1050_337_byoperator.txt", sep="\t", header=T)
Rpc_aflp<-read.table("Rpc_AFLP_CJJ_1050_337_byoperator.txt", sep="\t", header=T)
# load to genind objects
#Rpcf.gind<-df2genind(X=Rpcf_aflp[,5:539], sep=NULL, ncode=1, ind.names=Rpcf_aflp[,4], pop=Rpcf_aflp[,3], type="PA")
#Rpc.gind<-df2genind(X=Rpcf_aflp[61:232,5:539], sep=NULL, ncode=1, ind.names=Rpcf_aflp[61:232,4], pop=Rpcf_aflp[61:232,3], type="PA")
Rpcf.gind<-df2genind(X=Rpcf_aflp[,5:341], sep=NULL, ncode=1, ind.names=Rpcf_aflp[,4], pop=Rpcf_aflp[,2], type="PA")
Rpc.gind<-df2genind(X=Rpc_aflp[,5:341], sep=NULL, ncode=1, ind.names=Rpc_aflp[,4], pop=Rpc_aflp[,2], type="PA")


# build tree inclued Rpc and Rf
Rpcf.dist<-nei.dist(Rpcf.gind)
Rpcf.nj<-nj(Rpcf.dist)
Rpcf.loci<-as.loci(Rpcf.gind)
Rpcf.loci<-Rpcf.loci[,-1]
Rpcf_boot1000.nei<-boot.phylo(Rpcf.nj, Rpcf.loci, function(x) nj(nei.dist(loci2genind(x))), B=1000)
Rpcf_p <-character(length(Rpcf_boot1000.nei))
c <-c("green", "red", NULL)
Rpcf_p[Rpcf_boot1000.nei >= 700] <-c[1]
Rpcf_p[Rpcf_boot1000.nei < 700 & Rpcf_boot1000.nei >= 500 ] <-c[2]
Rpcf_p[Rpcf_boot1000.nei < 500] <-c[3]
#Rpcf.tip.color<-rep(rainbow(14), times=as.character(table(Rpcf.gind@pop)))
Rpcf.tip.color<-rep(rainbow(2), times=as.character(table(Rpcf_aflp$operator)))
# tree out
#tiff(file="Rpc_Rf_AFLP_1550_NJ.tif",width=6000,height=6000,res=600, compression = "lzw")
#tiff(file="Rpc_Rf_AFLP_CJJ_1050_NJ.tif",width=6000,height=6000,res=600, compression = "lzw")
tiff(file="Rpc_Rf_AFLP_CJJ_1050_byoperator_NJ.tif",width=6000,height=6000,res=600, compression = "lzw")
plot(Rpcf.nj, edge.width= 1, cex= 0.35, adj= 0.55, no.margin= TRUE, tip.color=Rpcf.tip.color)
nodelabels(pch= 20, col = Rpcf_p, cex=0.5)
dev.off()
# unrooted tree out
#tiff(file="Rpc_Rf_AFLP_1550_urNJ.tif",width=5000,height=4000,res=600, compression = "lzw")
#tiff(file="Rpc_Rf_AFLP_CJJ_1050_urNJ.tif",width=5000,height=4000,res=600, compression = "lzw")
tiff(file="Rpc_Rf_AFLP_CJJ_1050_byoperator_urNJ.tif",width=5000,height=4000,res=600, compression = "lzw")
plot(Rpcf.nj, type="unrooted", edge.width= 1, cex= 0.45, adj= 0.55, no.margin= TRUE, tip.color=Rpcf.tip.color)
nodelabels(pch= 20, col = Rpcf_p, cex=0.35)
dev.off()



# build tree only Rpc
Rpc.dist<-nei.dist(Rpc.gind)
Rpc.nj<-nj(Rpc.dist)
Rpc.loci<-as.loci(Rpc.gind)
Rpc.loci<-Rpc.loci[,-1]
Rpc_boot1000.nei<-boot.phylo(Rpc.nj, Rpc.loci, function(x) nj(nei.dist(loci2genind(x))), B=1000)
Rpc_p <-character(length(Rpc_boot1000.nei))
c <-c("green", "red", NULL)
Rpc_p[Rpc_boot1000.nei >= 700] <-c[1]
Rpc_p[Rpc_boot1000.nei < 700 & Rpc_boot1000.nei >= 500 ] <-c[2]
Rpc_p[Rpc_boot1000.nei < 500] <-c[3]
#Rpc.tip.color<-rep(rainbow(14), times=as.character(table(Rpc.gind@pop)))
Rpc.tip.color<-rep(rainbow(2), times=as.character(table(Rpc_aflp$operator)))
# tree out
#tiff(file="Rpc_AFLP_CJJ_1550_NJ.tif",width=4000,height=5000,res=600, compression = "lzw")
tiff(file="Rpc_AFLP_CJJ_1050_operator_NJ.tif",width=4000,height=5000,res=600, compression = "lzw")
plot(Rpc.nj, edge.width= 1, cex= 0.35, adj= 0.55, no.margin= TRUE, tip.color=Rpc.tip.color)
nodelabels(pch= 20, col = Rpc_p, cex=0.5)
dev.off()
# unrooted tree out
#tiff(file="Rpc_AFLP_1550_urNJ.tif",width=5000,height=5000,res=600, compression = "lzw")
tiff(file="Rpc_AFLP_CJJ_1050_operator_urNJ.tif",width=5000,height=5000,res=600, compression = "lzw")
plot(Rpc.nj, type="unrooted", edge.width= 1, cex= 0.45, adj= 0.55, no.margin= TRUE,tip.color=Rpc.tip.color)
nodelabels(pch= 20, col = Rpc_p, cex=0.35)
dev.off()
