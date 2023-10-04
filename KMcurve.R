
#km curve
remove(list=ls())
getwd()
library(survival)
library(survminer)
input=read.csv("ospfs.csv")
colnames(input)
os$Status=as.factor(os$Status)
os=survfit(Surv(os_time,os_status)~1,data=input)
ggsurvplot(os,data=input,palette = c("red"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="OS",
           pval = F,
           xlim=c(0,24),
           break.time.by=2)
pfs=survfit(Surv(pfs_time,pfs_status)~1,data=input)
ggsurvplot(pfs,data = input,palette =c("#27408B"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="PFS",
           xlim=c(0,14),
           break.time.by=2)


#READ
library(survival)
library(survminer)
data=read.csv("READ.csv")
colnames(data)
ostmb=survfit(Surv(os,os_status)~TMB10,data=data)
ggsurvplot(ostmb,data = data)
pfstmb=survfit(Surv(pfs,pfs_status)~TMB10,data=data)
ggsurvplot(pfstmb,data = data)
