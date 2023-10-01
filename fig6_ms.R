load("1kg_RMSE_FsM.RData")
library(corrplot)
source("fig_panel.R")

png("1kgFsts_new_500x500.png",width=500,height=500)
par(mfrow=c(2,2))
corrplot(fs.all$FsM[op,op],is.corr=FALSE,cl.pos="b",cl.length=7)
#corrplot(FstOS,is.corr=FALSE,tl.cex=0.5)

fig_label("A",cex=2)

fst2x2<-fs.all$Fst2x2
diag(fst2x2)<-0.0
corrplot(fst2x2[op,op],is.corr=FALSE,cl.pos="b",cl.length=7)

fig_label("B",cex=2)

#something for pop specific Fst here
boxplot(lapply(ls(pattern="Rmse")[c(3:5,1:2,9,7,6,8)],get),xlab="sampling scheme",ylab="RMSE Fst",names=c("2.2m L","220k L","22k L","11k L","2.2k L ","50 i","20 i","10 i","5 i"),las=2,col=rep(c("blue","red"),c(5,4)))

abline(v=5.5,lwd=2)

fig_label("C",cex=2)

# boxplot(lapply(ls(pattern="Rmse")[9+c(3:5,1:2,9,7,6,8)],get),xlab="sampling scheme",ylab="RMSE Fstp",names=c("2.2m","220k","22k","11k","2.2k","50","20","10","5"),las=2,col=rep(c("blue","red"),c(5,4)),ylim=c(1e-5,1e-1))


#par(mfrow=c(2,2))
#plot(fs.all$FsM,FsM.22k[[1]],pch=16,col="red",xlim=c(-0.2,0.3),ylim=c(-0.2,0.3),xlab="FST all",ylab="FST 22k SNPs")
#abline(c(0,1))
#for (i in 2:100) points(fs.all$FsM,FsM.22k[[i]],pch=16,col="red")

#plot(fs.all$Fst2x2,get2x2fromFsM(FsM.22k[[1]]),pch=16,col="red",xlim=c(0,0.18),ylim=c(0,0.18),xlab="FSTp all",ylab="FSTp 22k SNPs")
#abline(c(0,1))
#for (i in 2:100) points(fs.all$Fst2x2,get2x2fromFsM(FsM.22k[[i]]),pch=16,col="red")

plot(fs.all$FsM,Fsm.ni10[[1]],pch=16,col="red",xlim=c(-0.2,0.3),ylim=c(-0.2,0.3),xlab="FST all",ylab="FST 10 inds")
abline(c(0,1))
for (i in 2:100) points(fs.all$FsM,Fsm.ni10[[i]],pch=16,col="red")

fig_label("D",cex=2)

#plot(fs.all$Fst2x2,get2x2fromFsM(Fsm.ni10[[1]]),pch=16,col="red",xlim=c(0,0.18),ylim=c(0,0.18),xlab="FSTp all",ylab="FSTp 10 inds")
#abline(c(0,1))
#for (i in 2:100) points(fs.all$Fst2x2,get2x2fromFsM(Fsm.ni10[[i]]),pch=16,col="red")

par(mfrow=c(1,1))


dev.off()

