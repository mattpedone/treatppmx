#############################################################################
## This is used to simulate 38*4=152 subjects with 92 features
## using the top 1000 varied genes;

#############################################################################
rm(list=ls());
set.seed(20101027);
library("gtools");
#    sessionInfo()

# read raw data
brunet<-read.table("data/data-raw/bionmfBrunet.txt")[,c(1:38)];
rowtrack<-c(1:5000);

# pick the most varied genes
mmmax<-apply(brunet,1,max);  mmmin<-apply(brunet,1,min)
diff<-mmmax-mmmin;

brunet<-cbind(rowtrack,brunet) ### track genes with rowtrack variable
brunet<-brunet[order(diff),];
test1data<-brunet[4001:5000,];

# standarize the gene expressions
test2data<-cbind(test1data[,1],t(scale(t(test1data[,2:39]))));

# transform the data into positive matrix, +0.1 for 0s
mmmin2<-apply(test2data[,2:39],1,min)
for (i in 1:length(mmmin2)){
  test2data[i,2:39]<-test2data[i,2:39]-mmmin2[i]+0.1
}

# hcluster genes using correlation matrix
genecor<-round(cor(t(test2data[,2:39]),method="spearman"),2);
dissimilarity<-1-genecor;
distance <- as.dist(dissimilarity);
hcgen<-hclust(dist( distance),method="average");
#plot(hcgen,main="Dissimilarity = 1 - Correlation",xlab="gene");

# find the gene groups with high correlations.
memb<-cutree(hcgen,k=200);memb;
table(memb)/sum(table(memb));
seltemp<-subset(table(memb),table(memb)>=4);rownames(seltemp)
length(seltemp);
grpnumb<-as.numeric(rownames(seltemp));

# find groups that each has at least 4 genes
tempdata<-cbind(memb,test2data) ##sign the groups;
mydatatemp<-matrix(0,ncol= length(seltemp),nrow=38*4);

for(k in 1:length(seltemp)){
  ## find the subgroups
  tempsubdata<-subset(tempdata,tempdata[,1]==grpnumb[k]);
  subgencor<-round(cor(t(tempsubdata[,3:40]),method="spearman"),2);
  apply(subgencor,1,sum);
  ## sort the data with spearman's correlations
  tempcor<-rbind(apply(subgencor,1,sum),subgencor);
  subgencor2<-tempcor[,order(-tempcor[1,])];
  ## find the data with the first four highly correlated genes
  subtemp<-t(tempsubdata[,3:40])[,order(-tempcor[1,])];
  subtemp2<-subtemp[,1:4]
  ## extract the data
  mydatatemp[,k]<-matrix(subtemp2,ncol=1,byrow=TRUE);
}

mydatagen<-t(mydatatemp)

#write.table(mydatagen,file="data/simu152pats.txt")
