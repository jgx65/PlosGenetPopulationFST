#fig 3 ms popFST
load("FIM_M4.RData")
library(corrplot)
library(gaston)
source("fig_label.R")

png("rmse_FIM_M4_new2.png",height=500,width=500)
par(mfrow=c(3,2))
corrplot(ebeta[[20]],is.corr=FALSE)

fig_label("A",cex=2)


#hist of RMSEs as a function of nb loci
boxplot(lapply(ls(pattern="Rmse")[c(3,8,9,6,4,7,5)],get),xlab="sampling scheme",ylab="RMSE",names=c("10k L","1k L","AFM","20 i","10 i","5 i","2 i"),las=2,main="",col=rep(c("red","black","blue"),c(2,1,4)),ylim=c(0,0.04))
abline(v=c(2.5,3.5),lwd=2)

fig_label("B",cex=2)



#hist(Rmses.10kl,breaks=0:150/1000,col="green",xlab="RMSE",main="")
#hist(Rmses.1kl,breaks=0:150/1000,col="blue",add=TRUE)
#hist(Rmses.100l,breaks=0:150/1000,col="orange",add=TRUE)
#hist(Rmses.100l.raf,breaks=0:150/1000,col="red",add=TRUE)
#legend("topright",legend=c("10k","1k","100","100AFM"),col=c("green","blue","orange","red"),lwd=2)
plot(ebeta[[20]],FsM.10kl[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.10kl[[i]],pch=16,col="red")

fig_label("C",cex=2)



plot(ebeta[[20]],FsM.10kl.ni2[[1]],pch=16,col="blue",ylim=range(unlist(FsM.10kl.ni2)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.10kl.ni2[[i]],pch=16,col="blue")


fig_label("D",cex=2)



plot(ebeta[[20]],FsM.1kl[[1]],pch=16,col="red",ylim=range(unlist(FsM.1kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.1kl[[i]],pch=16,col="red")

fig_label("E",cex=2)




plot(ebeta[[20]],FsM.1kl.raf[[1]],pch=16,col="black",ylim=range(unlist(FsM.1kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.1kl.raf[[i]],pch=16,col="black")

fig_label("F",cex=2)



par(mfrow=c(1,1))
dev.off()

