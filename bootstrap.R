#!/usr/bin/env Rscript

library(boot)


args=(commandArgs(TRUE))

x<-c(rep(1,args[1]),rep(0,args[2]));
indices=1:length(x)
n=length(x)
xx=unlist(lapply(1:1000, function(y){
  b.x=x[sample(indices, n, replace=T)]
  return (mean(b.x))
}))
q<-quantile(xx, probs=c(0.025, 0.975),names=FALSE);

cat(q[1],"\n");
cat(q[2],"\n");

#div<-function(x, i) {
#  data<-x[i];
#  return (sum(data==0)/ (length(data) ));
#} 



#res<-boot(data=x, statistic=div, R=1000);
#res2<-boot.ci(res, type="bca");




#cat(res2$bca[4],"\n");
#cat(res2$bca[5],"\n");


