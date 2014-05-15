#!/usr/bin/env Rscript

library(boot)


args=(commandArgs(TRUE))


success<-strtoi(args[1]);
total<-(strtoi(args[1])+strtoi(args[2]));


b<-binom.test(success,total);
cat(b$conf[1],"\n");
cat(b$conf[2],"\n");


