#Fig 4 MS PopFst

load("ss1dM20.RData")
source("fig_label.R")
library(corrplot)
library(gaston)

png("rmse_ss1dM20_new2_1000x1000.png",height=1000,width=1000)

par(mfrow=c(3,2))
corrplot(ebeta[[23]],is.corr=FALSE,cl.pos="b",cl.length=7)

fig_label("A",cex=2)

#par(xpd=NA)
#text(-0.3,10,"A",cex=2)
#par(xpd=FALSE)


#hist of RMSEs as a function of nb loci
boxplot(lapply(ls(pattern="Rmse")[c(3,8,1,2,6,4,7,5)],get),xlab="sampling scheme",ylab="RMSE",names=c("10k L",",1k L","100 L","AFM","20 i","10 i","5 i","2 i"),las=2,main="",col=rep(c("red","black","blue"),c(3,1,4)),ylim=c(0,0.09))
abline(v=c(3.5,4.5),lwd=2)

fig_label("B",cex=2)
#par(xpd=NA)
#text(0,0.10,"B",cex=2)
#par(xpd=FALSE)

plot(ebeta[[23]],FsM.10kl[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[23]],FsM.10kl[[i]],pch=16,col="red")

fig_label("C",cex=2)

#par(xpd=NA)
#text(-0.08,0.2,"C",cex=2)
#par(xpd=FALSE)


plot(ebeta[[23]],FsM.10kl.ni2[[1]],pch=16,col="blue",ylim=range(unlist(FsM.10kl.ni2)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[23]],FsM.10kl.ni2[[i]],pch=16,col="blue")

fig_label("D",cex=2)

#par(xpd=NA)
#text(-0.08,0.2,"D",cex=2)
#par(xpd=FALSE)


plot(ebeta[[23]],FsM.100l[[1]],pch=16,col="red",ylim=range(unlist(FsM.100l)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[23]],FsM.100l[[i]],pch=16,col="red")

fig_label("E",cex=2)

#par(xpd=NA)
#text(-0.08,0.3,"E",cex=2)
#par(xpd=FALSE)



plot(ebeta[[23]],FsM.100l.raf[[1]],pch=16,col="black",ylim=range(unlist(FsM.100l)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[23]],FsM.100l.raf[[i]],pch=16,col="black")

fig_label("F",cex=2)

#par(xpd=NA)
#text(-0.08,0.09,"F",cex=2)
#par(xpd=FALSE)


par(mfrow=c(1,1))
dev.off()

