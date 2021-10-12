### Function to calculate the NPC
    countUT<-function(resultsum, trtsgn, myoutot){
      myctut<-array(0,dim=c(3,3,100));
      myctutSum<- NULL;
      for(i in 1:(dim(resultsum)[3])){
            mycurdata<-resultsum[,,i];
            mypre<-NULL;
            pretrt1<-apply(mycurdata[,4:6],1,which.max);
            pretrt2<-apply(mycurdata[,7:9],1,which.max);
            mypreTall<-cbind(pretrt1,pretrt2);
            for (j in 1:(dim(resultsum)[1])){mypre[j]<- mypreTall[j,trtsgn[j]]};
####################################################################################
## in case of sparse tables we do the following steps.
           sts<-table(mypre,myoutot[,i]);
           mysdls<-as.numeric(rownames(sts));
           str1<-matrix(0,nrow=3,ncol=3);
           str1[mysdls,]<-sts;
####################################################################################
            myctut[,,i]<-str1*diag(3);
            myctutSum[i]<-sum(str1*diag(3));

   }
   return<-cbind(myctutSum);
 }
