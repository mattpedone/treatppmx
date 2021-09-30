##################################################################################
## A total of 92 genes are selected for 152 simuluated patients               ####
## Also simulate 100 realizations data sets                                   ####
## Scenario 1                                                                 ####
## Junsheng Ma; 13 April, 2016                                                ####
##################################################################################
rm(list=ls()); set.seed(123456);
library("gtools");

##################################################################################
gene.norm<-t(scale(t(read.table(file="data/Simu152pats.txt"))));

##################################################################################
mypca<-prcomp(t(gene.norm[-c(11:92),]));
noise <- matrix(rep(rnorm(152), 10), nrow = 152, 10)
noisy_pbm <- cbind(t(gene.norm[-c(11:92),]), noise)
nsub<-length(mypca$x[,1]);
###################################################################################
#### use the combination of the pca to generate the exponetial
metx1<-(scale(mypca$x[,1])+scale(mypca$x[,2])+0.50);
metx<-sign(metx1)*(sign(metx1)*metx1)^(0.2)*0.45

###################################################################################
##### function of generating outcome variables
genoutcome<-function(nsub,alpha1,beta1,beta2,beta3,metx,x2,x3){
  nsub<-nsub;   alpha1<-alpha1;   beta1<-beta1;
  metx<-metx;  beta2<-beta2;   beta3<-beta3;
  prob<-matrix(0,nrow=nsub,ncol=3);
  myy<-rep(0,nsub);
  for (i in 1:nsub){
    temp<-NULL;
    temp<-beta1*metx[i]+beta2*x2[i]+beta3*x3[i];
    eta1<-alpha1[1]+temp[1];
    eta2<-alpha1[2]+temp[2];

    prob[i,1]<-1/((exp(eta1)+1)*(exp(eta2)+1));
    prob[i,2]<-exp(eta1)/((exp(eta1)+1)*(exp(eta2)+1));
    prob[i,3]<-exp(eta2)/(exp(eta2)+1);
    myy[i]<-sample(c(0,1,2),size=1,prob=prob[i,]);
  }
  return<-round(cbind(myy,prob),3);return;
}
##################################################################################
#### parametrs; ## prob with preditive features only for treatment 1;
x2<-gene.norm[91,]; x3<-gene.norm[92,];
myx2<-sign(x2)*(sign(x2)*x2)^(0.5);
myx3<-sign(x3)*(sign(x3)*x3)^(0.2);
newx<-cbind(rnorm(n=nsub,mean=0,sd=1),rnorm(n=nsub,mean=0,sd=3));
##################################################################################
alpha11<--1*c(0.5,1);  beta11<-1*c(2,2.6);
prob1<-genoutcome(nsub,alpha11,beta11,c(0,0),c(0,0),metx,myx2,myx3);
myw<-c(0,40,100);

## prob with preditive features only for treatment
alpha21<--1*c(-0.7,1);  beta21<--1*c(1,3);
prob2<-genoutcome(nsub,alpha21,beta21,c(0,0),c(0,0),metx,myx2,myx3);

## prob with prognostic only
alpha30<--1*c(-1,0.5);beta2<-1*c(1,0.5);beta3<-1*c(0.7,1)
prob_prog<-genoutcome(nsub,alpha30,c(0,0),beta2,beta3,metx,myx2,myx3);
## prob with prog and pred features
myprob1<-myprob2<-matrix( 0,nrow=nsub,ncol=3);

for (i in 1:nsub){
  myprob1[i,]<-prob1[i,2:4]*prob_prog[i,2:4]/sum(prob1[i,2:4]*prob_prog[i,2:4]);
  myprob2[i,]<-prob2[i,2:4]*prob_prog[i,2:4]/sum(prob2[i,2:4]*prob_prog[i,2:4]);
}
##################################################################################
### simulate the outcomes with the above parameters.
nset<-100;
myprob<-list( myprob1, myprob2);    trtsgn<-rep(c(1,2),nsub/2);#trtsgn<-sample(c(1,2),nsub, replace = T);#
###############Generat outcomes for replicated dataset ################################
mytot<-array(0,dim=c(nsub,3,nset));myoutot<-matrix(0,nrow=nsub,ncol= nset);
for (i in 1: nset){
  myy<-matrix(0,nrow=nsub,ncol=3);myyout<-matrix(0,nrow=nsub,ncol=1);
  for (k in 1:nsub){
    trtemp<-trtsgn[k];
    if (trtemp==1) {myy[k,1:3]<-t(rmultinom(n=1,size=1,prob= myprob[[1]][k,]));}
    if (trtemp==2) {myy[k,1:3]<-t(rmultinom(n=1,size=1,prob=myprob[[2]][k,]));}
    myyout[k]<-match(1,myy[k,]);
    trtemp<-NULL;
  }
  mytot[,,i]<-myy; myoutot[,i]<-myyout;
}
###################################################################################

mydata<-gene.norm;
orgx<-cbind(x2,x3);
    save(myoutot,mytot,mydata, noisy_pbm, trtsgn, myprob,orgx,myx2,myx3,newx,file="data/SimuOutsce2_noise.rda")






