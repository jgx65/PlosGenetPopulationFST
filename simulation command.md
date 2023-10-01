The code necessary to generate the simulated data set and to estimate the FST matrix is given below


## Island model

`sim.genot.metapop` commands to simulate the island model

```
m<-0.001
mut<-1e-6
np<-10
N<-c(1000,1000,10,10,100,100,500,500,2000,2000)
nl<-10000
FIM<-matrix(rep(m/(np-1),np^2),ncol=np)
diag(FIM)<-NA
diag(FIM)<-1-rowSums(FIM,na.rm=TRUE)
dat<-hierfstat::sim.genot.metapop.t(nbal=2,nbloc=10000,nbpop=np,N=N,mig=M,mut=1e-6,t=2000)

```

## Stepping-stone

`sim.genot.metapop` commands to simulate the stepping-stone

```
np<-10
N<-rep(1e+3,10)
m<-20/4/N[1]
ss1d20<-matrix(numeric(np^2),ncol=np)
diag(ss1d20[-1,-np])<-diag(ss1d20[-np,-1])<-m
diag(ss1d20)<-1-2*m
ss1d20[1,1]<-ss1d20[np,np]<-1-m
M<-ss1d20

dat<-hierfstat::sim.genot.metapop.t(nbal=2,nbloc=10000,nbpop=np,N=N,mig=M,mut=1e-6,t=4000)
```


## River system

`sim.genot.metapop` commands to simulate the river system:

```
md=0.02
mu=0.005
N<-100*c(1,1,2,2,4,1,1,1,2,2,2,8,10,50)
M<-matrix(numeric(14^2),ncol=14)
M[1,2]<-M[2,3]<-M[3,4]<-M[4,5]<-M[5,12]<-M[12,13]<-M[13,14]<-md
M[6,7]<-M[7,8]<-M[8,3]<-md
M[9,10]<-M[10,11]<-M[11,5]<-md

M[2,1]<-M[3,2]<-M[4,3]<-M[5,4]<-M[12,5]<-M[13,12]<-M[14,13]<-mu
M[7,6]<-M[8,7]<-M[3,8]<-mu
M[10,9]<-M[11,10]<-M[5,11]<-mu
diag(M)<-1-rowSums(M)

np<-14

 hierfstat::sim.genot.metapop.t(nbal=2,nbloc=10000,nbpop=np,N=N,mig=M,mut=1e-6,t=2000)
```

`mspms` commands to simulate the river system:

```
mspms 1400 20 -t 4000 -r 4000 100000000 -I 15 100 100 100 100 100 100 100 100 100 100 100 100 100 100 0 -m 1 2 80 -m 2 3 80 -m 3 4 80 -m 4 5 80 -m 5 12 80 -m 12 13 80 -m 13 14 80 -m 6 7 80 -m 7 8 80 -m 8 3 80 -m 9 10 80 -m 10 11 80 -m 11 5 80 -m 2 1 20 -m 3 2 20 -m 4 3 20 -m 5 4 20 -m 12 5 20 -m 13 12 20 -m 14 13 20 -m 7 6 20 -m 8 7 20 -m 3 8 20 -m 10 9 20 -m 11 10 20 -m 5 11 20 -n 1 0.1 -n 2 0.1 -n 3 0.2 -n 4 0.2 -n 5 0.4 -n 6 0.1 -n 7 0.1 -n 8 0.1 -n 9 0.2 -n 10 0.2 -n 11 0.2 -n 12 0.8 -n 14 5 -ej 2 1 15 -ej 2 2 15 -ej 2 3 15 -ej 2 4 15 -ej 2 5 15 -ej 2 6 15 -ej 2 7 15 -ej 2 8 15 -ej 2 9 15 -ej 2 10 15 -ej 2 11 15 -ej 2 12 15 -ej 2 13 15 -ej 2 14 15 -p 9
```

## R commands to obtain the FST matrix

From the output of `sim.genot.metapop.t`, stored in object `dat`

```
dos<-hierfstat::biall2dos(dat[,-1])
FST.1<-hierfstat::fs.dosage(dos,pop=dat[,1])$FsM
```
From the output of `mspms` stored in file `msres.txt`

```
bed<-hierfstat::ms2bed("msres.txt")
M<-hierfstat::matching(bed)
FST.2<-hierfstat::fs.dosage(M,pop=dat[,1],matching=TRUE)$FsM
```

## R commands to estimate FST using AFM:

```
das<-data.frame(pop=dat[,1],cbind(dat[,-1]%/%10,dat[,-1]%%10)
res.rafm<-rafm::do.all(das,15000,5000,5)

#mean of the 200 last elements of each MCMC chains for coancestries
mean.coan.rafm<-apply(res.rafm$theta[,,1801:2000],c(1,2),mean)
MB<-mean(hierfstat::mat2vec(mean.coan.rafm))
FST.RAF<-(mean.coan.rafm-MB)/(1-MB)
```