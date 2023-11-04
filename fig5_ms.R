load("rivsys.RData")
library(corrplot)
library(gaston)
source("fig_label.R")

#png("rmse_riversystem_new2_500x500.png",width=500,height=500)
#Obs vs expected for 10k, 1k and 100 loci for the whole matrix
par(mfrow=c(3,2))
corrplot(ebeta[[20]],is.corr=FALSE,cl.pos="b",cl.length=7)

fig_label("A",cex=2)

#hist of RMSEs as a function of nb loci
boxplot(lapply(ls(pattern="Rmse")[c(12,1,2,3,6,11,4,5,9,7,10,8)],get),xlab="sampling scheme",ylab="RMSE",names=c("all.ms","100k L.ms","10k L.ms","1k L.ms","10k L","1k L","100 L","AFM","20 i","10 i","5 i","2 i"),las=2,main="",col=rep(c("orange","red","black","blue"),c(4,3,1,4)),ylim=c(0,0.05));abline(v=c(4.5,7.5,8.5),lwd=2)

fig_label("B",cex=2)

#hist(Rmses.10kl,breaks=0:150/1000,col="green",xlab="RMSE",main="")
#hist(Rmses.1kl,breaks=0:150/1000,col="blue",add=TRUE)
#hist(Rmses.100l,breaks=0:150/1000,col="orange",add=TRUE)
#hist(Rmses.100l.raf,breaks=0:150/1000,col="red",add=TRUE)
#legend("topright",legend=c("10k","1k","100","100AFM"),col=c("green","blue","orange","red"),lwd=2)

plot(ebeta[[22]],FsM.10kl[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[22]],FsM.10kl[[i]],pch=16,col="red")

fig_label("C",cex=2)


plot(ebeta[[20]],FsM.ms.10kl[[1]],pch=16,col="orange",ylim=range(unlist(FsM.ms.10kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.ms.10kl[[i]],pch=16,col="orange")

fig_label("D",cex=2)

plot(ebeta[[20]],FsM.100l[[1]],pch=16,col="red",ylim=range(unlist(FsM.100l)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.100l[[i]],pch=16,col="red")

fig_label("E",cex=2)

plot(ebeta[[20]],FsM.100l.raf[[1]],pch=16,col="black",ylim=range(unlist(FsM.100l)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.100l.raf[[i]],pch=16,col="black")

fig_label("F",cex=2)


par(mfrow=c(1,1))
#dev.off()

