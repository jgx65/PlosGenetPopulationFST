####Fig 1 MS

get.coan<-function(N=rep(1000,10),mut=1e-06,M=diag(10),ngen=100){
np<-length(N)
J<-matrix(rep(1,np^2),nrow=np)
#gens<-c(1:9*10,1:9*100,1:5*1000)
etheta<-NULL
a<-(1-mut)^2
ig<-0
theta <- matrix(numeric(np^2), nrow = np)
for (i in 1:ngen){
theta<-a*(M%*%((J-diag(1/2/N))*theta+diag(1/2/N))%*%t(M))
}
#if (i %in% gens) {ig<-ig+1;etheta[[ig]]<-theta}
#theta<-b*(M%*%theta%*%t(M)+diag(1/2/N))
theta
}

get.beta<-function(x) {tmp<-x;diag(tmp)<-NA;Mb<-mean(tmp,na.rm=TRUE);(x-Mb)/(1-Mb)}


breaks<--1:2+.5
#coul<-c("grey","red","blue")
#im_F<-function(x,main=""){
#image(1:16,1:16,x,xlab="",ylab="",
#      main=main,col=coul,breaks=breaks)
#abline(h=0:4*4+0.5,v=0:4*4+0.5)
#abline(h=0:16+0.5,v=0:16+0.5,lty=2)  
#}

#im_Fp<-function(x,main=""){
#image(1:4,1:4,x,xlab="",ylab="",
#      main=main,col=coul,breaks=breaks)
#abline(h=0:4+0.5,v=0:4+0.5)
#}

#coul<-c("grey","red","blue")
###F[ST]^l
#x<-matrix(2,nrow=16,ncol=16)
#x[5:8,5:8]<-1
#diag(x)<-0
#x[1:4,1:4]<-x[9:12,9:12]<-x[13:16,13:16]<-0
#im_F(x,main=expression(F[ST]^2))



#x<-matrix(2,nrow=4,ncol=4)
#diag(x)<-0
#x[2,2]<-1
#im_Fp(x,main=expression(F[ST]^2))

breaks<--1:2+.5

im_Fp<-function(x,main="",a=6){
image(1:a,1:a,t(x[a:1,]),xlab="",ylab="",
      main=main,col=coul,breaks=breaks,axes=FALSE)
abline(h=0:a+0.5,v=0:a+0.5)
}
########################################################
source("fig_label.R")
png("mig_coan_fst_new2.png",width=500,height=500)

mycols<-c("orange","red","brown","purple","blue","green")
par(mfrow=c(3,3))
#Cont-island
coul<-c("white","lightgrey","darkgrey")
x<-matrix(0,nrow=7,ncol=7)
diag(x)<-2
x[-1,1]<-1
im_Fp(x,main="",a=7)#Continent-island")

fig_label("A",cex=2)

#Fim
x<-matrix(1,nrow=6,ncol=6)
diag(x)<-2
#x[-1,1]<-1
im_Fp(x,main="")#"Finite island")

fig_label("B",cex=2)

x<-matrix(0,nrow=6,ncol=6)
diag(x)<-2
diag(x[-1,-6])<-diag(x[-6,-1])<-1
im_Fp(x,main="")#Stepping-stone")

fig_label("C",cex=2)

m<-0.001
gens<-c(1:9*10,1:9*100,1:5*1000)
N<-c(1e+9,10,50,100,500,1000,5000)
CIM<-diag(7)
CIM[-1,1]<-m
diag(CIM)[-1]<-1-CIM[-1,1]
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=CIM,ngen=x))
ebeta<-lapply(etheta,function(y) get.beta(y[-1,-1]))
etd<-matrix(unlist(lapply(etheta,function(x) diag(x)[-1])),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etd[,1],xlab="gens",ylab="Coan.",type="l",col="orange",lwd=2,ylim=c(0,1),
main="")#c(paste0("m: ",m)," N: ",paste0(N[-1],collapse="/")))
for (i in 2:6) lines(gens,etd[,i],col=mycols[i],lwd=2)

fig_label("D",cex=2)

N<-c(10,50,100,500,1000,5000)
FIM<-matrix(m/5,ncol=6,nrow=6)
diag(FIM)<-NA
diag(FIM)<-1-rowSums(FIM,na.rm=TRUE)
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=FIM,ngen=x))
ebeta<-lapply(etheta,get.beta)
etd<-matrix(unlist(lapply(etheta,function(x) diag(x))),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etd[,1],xlab="gens",ylab="Coan.",type="l",col="orange",lwd=2,ylim=c(0,1),
main="")#c(paste0("m: ",m)," N: ",paste0(N,collapse="/")))
for (i in 2:6) lines(gens,etd[,i],col=mycols[i],lwd=2)

fig_label("E",cex=2)

N<-c(10,50,1000,1000,100,100)
m<-0.01
SSM<-matrix(0,ncol=6,nrow=6)
diag(SSM[-1,-6])<-diag(SSM[-6,-1])<-m/2
diag(SSM)<-1-rowSums(SSM)
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=SSM,ngen=x))
ebeta<-lapply(etheta,get.beta)
etd<-matrix(unlist(lapply(etheta,function(x) diag(x))),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etd[,1],xlab="gens",ylab="Coan.",type="l",col="orange",lwd=2,ylim=c(0,1),main="")#c(paste0("m: ",m)," N: ",paste0(N,collapse="/")))
for (i in 2:6) lines(gens,etd[,i],col=mycols[i],lwd=2)

fig_label("F",cex=2)





m<-0.001
gens<-c(1:9*10,1:9*100,1:5*1000)
N<-c(1e+9,10,50,100,500,1000,5000)
CIM<-diag(7)
CIM[-1,1]<-m
diag(CIM)[-1]<-1-CIM[-1,1]
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=CIM,ngen=x))
ebeta<-lapply(etheta,function(y) get.beta(y[-1,-1]))
etd<-matrix(unlist(lapply(etheta,function(x) diag(x)[-1])),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etb[,1],xlab="gens",ylab="FST",type="l",col="orange",lwd=2,ylim=c(-0.1,1));for (i in 2:6) lines(gens,etb[,i],col=mycols[i],lwd=2)

fig_label("G",cex=2)



N<-c(10,50,100,500,1000,5000)
FIM<-matrix(m/5,ncol=6,nrow=6)
diag(FIM)<-NA
diag(FIM)<-1-rowSums(FIM,na.rm=TRUE)
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=FIM,ngen=x))
ebeta<-lapply(etheta,get.beta)
etd<-matrix(unlist(lapply(etheta,function(x) diag(x))),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etb[,1],xlab="gens",ylab="FST",type="l",col="orange",lwd=2,ylim=c(-0.1,1));for (i in 2:6) lines(gens,etb[,i],col=mycols[i],lwd=2)

fig_label("H",cex=2)



N<-c(10,50,1000,1000,100,100)
m<-0.01
SSM<-matrix(0,ncol=6,nrow=6)
diag(SSM[-1,-6])<-diag(SSM[-6,-1])<-m/2
diag(SSM)<-1-rowSums(SSM)
etheta<-lapply(gens,function(x) get.coan(N=N,mut=1e-8,M=SSM,ngen=x))
ebeta<-lapply(etheta,get.beta)
etd<-matrix(unlist(lapply(etheta,function(x) diag(x))),ncol=6,byrow=TRUE)
etb<-matrix(unlist(lapply(ebeta,function(x) diag(x))),ncol=6,byrow=TRUE)
plot(gens,etb[,1],xlab="gens",ylab="FST",type="l",col="orange",lwd=2,ylim=c(-0.1,1));for (i in 2:6) lines(gens,etb[,i],col=mycols[i],lwd=2)

fig_label("I",cex=2)
par(mfrow=c(1,1))
dev.off()



###########################
#par(mfrow=c(2,2))
#image(1:6,1:6,ebeta[[1]],xlab="",ylab="",main="FST, gen: 10");abline(h=0:6+0.5,v=0:6+0.5)
#image(1:6,1:6,ebeta[[10]],xlab="",ylab="",main="FST, gen: 100");abline(h=0:6+0.5,v=0:6+0.5)
#image(1:6,1:6,ebeta[[19]],xlab="",ylab="",main="FST, gen: 1000");abline(h=0:6+0.5,v=0:6+0.5)
#image(1:6,1:6,ebeta[[23]],xlab="",ylab="",main="FST, gen: 5000");abline(h=0:6+0.5,v=0:6+0.5)
#par(mfrow=c(2,2))


