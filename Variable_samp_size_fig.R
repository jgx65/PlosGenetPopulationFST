##assumes dos.g2k, 20 replicates of dosage data for 10'000 SNPs and 50 individuals in each of the 14 river-system populations obtained from sim.genot.metapop.t loaded in the environment

##unequal sample sizes
#repl 1
ss<-sample(rep(c(1:5*2,20,40),each=2),replace=FALSE)

FsM.10kl.niv1<-lapply(dos.g2k,function(x) {a<-rep(1:np-1,c(ss))*50+unlist(lapply(ss,function(x) sample(50,size=x)));hierfstat::fs.dosage(x[a,],pop=rep(1:np,ss))$FsM})

tmp<-lapply(FsM.10kl.niv1,function(x) (x-ebeta[[20]])^2)
Rmses.10kl.niv1<-matrix(numeric(np^2),ncol=np)
for (i in 1:20) Rmses.10kl.niv1<-Rmses.10kl.niv1+tmp[[i]]
Rmses.10kl.niv1<-(Rmses.10kl.niv1/20)^.5

#repl 2
ss<-sample(rep(c(1:5*2,20,40),each=2),replace=FALSE)

FsM.10kl.niv2<-lapply(dos.g2k,function(x) {a<-rep(1:np-1,c(ss))*50+unlist(lapply(ss,function(x) sample(50,size=x)));hierfstat::fs.dosage(x[a,],pop=rep(1:np,ss))$FsM})

tmp<-lapply(FsM.10kl.niv2,function(x) (x-ebeta[[20]])^2)
Rmses.10kl.niv2<-matrix(numeric(np^2),ncol=np)
for (i in 1:20) Rmses.10kl.niv2<-Rmses.10kl.niv2+tmp[[i]]
Rmses.10kl.niv2<-(Rmses.10kl.niv2/20)^.5

##subsampling populations

#rep1
rs1<-sample(14,size=7)
FsM.10kl.rs1<-lapply(FsM.10kl,function(x) {tmp<-mean(hierfstat::mat2vec(x[rs1,rs1]));(x[rs1,rs1]-tmp)/(1-tmp)})

#rep 2
rs2<-sample(14,size=7)
FsM.10kl.rs2<-lapply(FsM.10kl,function(x) {tmp<-mean(hierfstat::mat2vec(x[rs2,rs2]));(x[rs2,rs2]-tmp)/(1-tmp)})

#rep 3
rs3<-sample(14,size=7)
FsM.10kl.rs3<-lapply(FsM.10kl,function(x) {tmp<-mean(hierfstat::mat2vec(x[rs3,rs3]));(x[rs3,rs3]-tmp)/(1-tmp)})



##################################################
##assumes ebeta[[20]], the expected values for FST in the river system after 2000 generations loaded in the environment

png("Subsamp_ind_pop_500x500.png",width=500,height=500)
par(mfrow=c(2,3))
plot(ebeta[[20]],FsM.10kl.niv1[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl.niv1)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.10kl.niv1[[i]],pch=16,col="red")

fig_label("A",cex=2)

plot(ebeta[[20]],FsM.10kl.niv2[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl.niv2)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.10kl.niv2[[i]],pch=16,col="red")

fig_label("B",cex=2)

##assumes Rmses for 10k loci and 50, 20, 10, 5 and 2 individuals loaded in the environment

boxplot(as.vector(Rmses.10kl),as.vector(Rmses.10kl.ni20),as.vector(Rmses.10kl.ni10),as.vector(Rmses.10kl.ni5),as.vector(Rmses.10kl.ni2),as.vector(Rmses.10kl.niv1),as.vector(Rmses.10kl.niv2),ylab="RMSEs",names=c("50 i","20 i","10 i","5 i","2 i","var i","var i"),las=2,xlab="sampling scheme")

fig_label("C",cex=2)

tmp1<-ebeta[[20]][rs1,rs1]
tmp2<-mean(hierfstat::mat2vec(tmp1))
tmp1<-(tmp1-tmp2)/(1-tmp2)


plot(tmp1,FsM.10kl.rs1[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl.rs1)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(tmp1,FsM.10kl.rs1[[i]],pch=16,col="red")

fig_label("D",cex=2)


tmp1<-ebeta[[20]][rs2,rs2]
tmp2<-mean(hierfstat::mat2vec(tmp1))
tmp1<-(tmp1-tmp2)/(1-tmp2)


plot(tmp1,FsM.10kl.rs2[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl.rs2)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(tmp1,FsM.10kl.rs2[[i]],pch=16,col="red")

fig_label("E",cex=2)


tmp1<-ebeta[[20]][rs3,rs3]
tmp2<-mean(hierfstat::mat2vec(tmp1))
tmp1<-(tmp1-tmp2)/(1-tmp2)

plot(tmp1,FsM.10kl.rs3[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl.rs3)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(tmp1,FsM.10kl.rs3[[i]],pch=16,col="red")

fig_label("F",cex=2)


par(mfrow=c(1,1))

dev.off()
#################################
plot(ebeta[[20]],FsM.10kl[[1]],pch=16,col="red",ylim=range(unlist(FsM.10kl)),xlab="E[FST]",ylab="FST",main="");abline(c(0,1))
for (i in 2:20) points(ebeta[[20]],FsM.10kl[[i]],pch=16,col="red")

corrplot(FsM.10kl[[1]],is.corr=FALSE)

corrplot(FsM.10kl.niv[[1]],is.corr=FALSE)


