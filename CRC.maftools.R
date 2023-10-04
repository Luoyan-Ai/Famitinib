
library(maftools)       
setwd("D:\\桌面\\突变瀑布图")      

#读取分组文件
score=read.table("group.txt", header=T, sep="\t", check.names=F)


clinicaldata=score[,c(1, 2,3,4,5,6)]
colnames(clinicaldata)=c("Tumor_Sample_Barcode", "location","PFS_status",	"OS_status","BOR")
write.table(clinicaldata, file="ann.txt", sep="\t", quote=F, row.names=F)

#绘制

pdf(file="CRC_mution3333333333333333333333.pdf", width=6, height=6)
maf=read.maf(maf="CRC.maf", clinicalData="ann2.txt")
oncoplot(maf=maf, showTumorSampleBarcodes = T ,clinicalFeatures=c("BOR","location","PFS_status","OS_status"),sortByAnnotation = TRUE)
dev.off()













pdf(file="CRC_mution2233333333333.pdf", width=6, height=6)
maf=read.maf(maf="CRC.maf", clinicalData="ann1.txt")
oncoplot(maf=maf, showTumorSampleBarcodes = T ,clinicalFeatures="BOR",sortByAnnotation = TRUE)
dev.off()









pdf(file="CRC_mution22.pdf", width=6, height=6)
maf=read.maf(maf="CRC.maf", clinicalData="ann.txt")



col_location=pal_jama("default",alpha=0.8)(3)
assign_location=setNames(col_location,unique(clinicaldata$location))

col_BOR=pal_jama("default",alpha=0.8)(5)
assingn_BOR=setNames(col_BOR,unique(clinicaldata$BOR))
col_PFS=c("#743A3A","#484891")
assign_PFS=setNames(col_PFS,unique(clinicaldata$PFS))
col_OS=c("#743A3A","#484891")
assign_OS=setNames(col_OS,unique(clinicaldata$OS))
oncoplot(maf=maf, showTumorSampleBarcodes = T ,clinicalFeatures=c("BOR","location","PFS_status","OS_status"),sortByAnnotation = TRUE,
annotationColor = list(location=assign_location,BOR=assingn_BOR,PFS=assign_PFS,OS=assign_OS))
dev.off()



















location=pal_jama("default",alpha=0.8)(3)
assign_location=setNames(col_location,unique(clinicaldata$location))

col_BOR=pal_jama("default",alpha=0.8)(5)
assingn_BOR=setNames(col_BOR,unique(clinicaldata$BOR))
col_PFS=c("#743A3A","#484891")
assign_PFS=setNames(col_PFS,unique(clinicaldata$PFS))
col_OS=c("#743A3A","#484891")
assign_OS=setNames(col_OS,unique(clinicaldata$OS))



oncoplot(maf =maf, fontSize = 0.35 ,showTumorSampleBarcodes = T ,
clinicalFeatures = c('PD.L1_expression','bTMB','BOR','PFS','Smoking_status','OS'),
titleFontSize=1.2,legendFontSize=0.7,,SampleNamefontSize=0.6,removeNonMutated=F,writeMatrix=T,annotationFontSize=0.5,sortByAnnotation = TRUE,genes=gene,keepGeneOrder=TRUE,
sampleOrder=sampeOrder,
annotationColor = list(PD.L1_expression=assign_PD.L1_expression,bTMB=assign_bTMB,BOR=assingn_BOR,PFS=assign_PFS,OS=assign_OS,
Smoking_status=assign_Smoking_status))












#读取分组文件
score=read.table("group.txt", header=T, sep="\t", check.names=F)


outTab=score[,c(1, 2,3,4,5,6)]
colnames(outTab)=c("Tumor_Sample_Barcode", "location","PFS_status",	"OS_status","BOR","Value")
write.table(outTab, file="ann1.txt", sep="\t", quote=F, row.names=F)

#绘制

pdf(file="CRC_mution2111112.pdf", width=6, height=6)
maf=read.maf(maf="CRC.maf", clinicalData="ann1.txt")
oncoplot(maf=maf, showTumorSampleBarcodes = T ,clinicalFeatures="Value",,sortByAnnotation = TRUE)
dev.off()


















clin <- read.table("response_feature.txt",header=TRUE)
maf <- read.maf(maf = "CRC.maf", clinicalData = clin)#注意读取maf时需要增加临床表型信息




sampeOrder<-c("30005","15002","09003","33006","20007","22001","01005","23009","01039","15001","01008","20013","31007","01013",
"31006","01015","22005","15012","01036","01041","20009","09002","01016","01001",
"31008","01031","09005","28005","01009","31004","01022","01042","06007","01006","20005","22008","30009","20010","20018",
"31001","31005","36001","43002","43004")
pdf("All.oncoplot.clin.pdf", width=6, height=6)











clin <- read.table("response_feature.txt",header=TRUE)
pdf("All.oncoplot.clin.pdf", width=6, height=6)

maf <- read.maf(maf = "CRC.maf", clinicalData = clin)#注意读取maf时需要增加临床表型信息

oncoplot(maf =maf, fontSize = 0.45 ,showTumorSampleBarcodes = T ,clinicalFeatures = 'BOR',titleFontSize=1.2,legendFontSize=0.7,removeNonMutated=F,writeMatrix=T,annotationFontSize=0.7,sampleOrder=sampeOrder)

dev.off()



oncoplot(maf =maf, fontSize = 0.45 ,showTumorSampleBarcodes = T ,clinicalFeatures = 'Value',titleFontSize=1.2,legendFontSize=0.7,
removeNonMutated=F,writeMatrix=T,annotationFontSize=0.7,sortByAnnotation = TRUE)









3）样本排序
如果你想按照表达量高低对样本进行排序，也可以指定样本顺序：

oncoplot(maf =maf, fontSize = 0.45 ,showTumorSampleBarcodes = T ,clinicalFeatures = 'ORR',titleFontSize=1.2,legendFontSize=0.7,
removeNonMutated=F,writeMatrix=T,annotationFontSize=0.7,sortByAnnotation = TRUE)




