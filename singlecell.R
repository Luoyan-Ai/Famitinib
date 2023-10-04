library(Seurat)
library(readxl)
library(tidyverse)
library(devtools)
library(monocle)
library(ggplot2)
library(pheatmap)

a <- readRDS("41467_2022_29366_MOESM6_ESM.rds")
a <- as.matrix(a)

b <-  read_excel("41467_2022_29366_MOESM7_ESM.xlsx",1)
table(b$PatientID) 
table(b$Tissues) 

b.group <- b[,c("PatientID","Tissues","MainTypes","Barcode")]
head(b.group)
identical(colnames(a),b.group$Barcode)#a中样本是否都在b中表达

index <- which(b.group$Tissues =="T")
length(index)
group <- b.group[index,] #筛选的是行
head(group)
a.filt <- a[,index] #筛选的是列
dim(a.filt)
identical(colnames(a.filt),group$Barcode)

sce.meta <- data.frame(Patient_ID=group$PatientID,
                       row.names = group$Barcode,
                       celltype = group$MainTypes)
head(sce.meta)
table(sce.meta$Patient_ID)

sce = CreateSeuratObject(counts=a.filt,
                         meta.data = sce.meta,
                         min.cells = 3, 
                         min.features = 50)
head(sce@meta.data)
#nCount_RNA：the number of cell total counts
#nFeature_RNA：the number of cell's detected gene
summary(sce@meta.data)
sce@assays$RNA@counts[1:4,1:4]

dim(sce)

table(grepl("^MT-",rownames(sce)))
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
head(sce@meta.data)
summary(sce@meta.data)
pctMT=15
sce <- subset(sce, subset = percent.mt < pctMT)
dim(sce)

table(grepl("^ERCC-",rownames(sce)))
sce[["percent.ERCC"]] <- PercentageFeatureSet(sce, pattern = "^ERCC-")
head(sce@meta.data)
summary(sce@meta.data)
rownames(sce)[grep("^ERCC-",rownames(sce))]

sum(sce$percent.ERCC< 10)  
pctERCC=10
sce <- subset(sce, subset = percent.ERCC < pctERCC)
dim(sce)
dim(a.filt)

col.num <- length(unique(sce@meta.data$Patient_ID))
library(ggplot2)

p1_1.1 <- VlnPlot(sce,
                  features = c("nFeature_RNA"),
                  group.by = "Patient_ID",
                  cols =rainbow(col.num)) +
  theme(legend.position = "none") +
  labs(tag = "A")
p1_1.1
p1_1.2 <- VlnPlot(sce,
                  features = c("nCount_RNA"),
                  group.by = "Patient_ID",
                  cols =rainbow(col.num)) +
  theme(legend.position = "none") 
p1_1.2
p1_1 <- p1_1.1 | p1_1.2
p1_1
VlnPlot(sce,
        features = c("nFeature_RNA","nCount_RNA","percent.ERCC"))
p1_2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by = "Patient_ID",pt.size = 1.3) +
  labs(tag = "B")
p1_2
FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.ERCC")
sessionInfo()
sce$celltype <- sce$orig.ident
sce$orig.ident <- sce$Patient_ID
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 1000) 
top10 <- head(VariableFeatures(sce), 10) 
top10
plot1 <- VariableFeaturePlot(sce) 
p1_3 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) +
  theme(legend.position = c(0.1,0.8)) +
  labs(tag = "C")
p1_3 

library(scater)
sce  <- sce %>% NormalizeData %>% ScaleData %>% RunPCA
DimPlot(sce, reduction = "pca",group.by = "orig.ident")
pc.num=1:20
library(harmony)
sce <-sce %>% RunHarmony("orig.ident", plot_convergence = TRUE)
sce <- JackStraw(sce,reduction = "pca", dims=20)
sce <- ScoreJackStraw(sce,dims = 1:20)

p2_2 <- JackStrawPlot(sce,dims = 1:20, reduction = "pca") +
  theme(legend.position="bottom") +
  labs(tag = "E")
p2_2
p2_3 <- ElbowPlot(sce, ndims=20, reduction="pca") 
p2_3

sce$patienttype<-NA
sce$patienttype[which(sce$Patient_ID == "p1")] = "COAD"
sce$patienttype[which(sce$Patient_ID == "p2")] = "COAD"
sce$patienttype[which(sce$Patient_ID == "p3")] = "READ"
sce$patienttype[which(sce$Patient_ID == "p4")] = "READ"
sce$patienttype[which(sce$Patient_ID == "p5")] = "READ"

sce <- FindNeighbors(sce, dims = pc.num,reduction = "harmony") 
sce <- FindClusters(sce, resolution = 0.5)

table(sce@meta.data$seurat_clusters)
sce <- RunUMAP(sce,reduction = "harmony",dims = pc.num)
UMAPPlot(sce, reduction = "harmony",label=F,group.by = "celltype")

Tcell <- subset(sce,subset = celltype == c("T","ProliferatingT"))
Tcell  <- Tcell %>% FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% NormalizeData %>% ScaleData %>% RunPCA

Tcell <- FindNeighbors(Tcell, dims = pc.num,reduction = "harmony") 
Tcell <- FindClusters(Tcell, resolution = 1.7)
table(Tcell$seurat_clusters)

Tcell <- RunUMAP(Tcell,reduction = "harmony",dims = pc.num)
UMAPPlot(Tcell, reduction = "harmony",label=F,group.by = "seurat_clusters")


cellmarker <- c('CD4','CD8A','CD8B','CD3D','CD3E','CD3G',
                              'TIGIT','LAG3','MIR155HG','LAYN','HAVCR2','TOX','PDCD1','ITGAE',#Inhibitory
                              'CD69','CXCR4','RUNX3','NR4A1',#resident
                              'CCR7','TCF7','LEF1','SELL',#NAIVE 
                              'CD44','IFNG','IL2','LAMTOR3','GNLY',#cytokines
                              'TNFRSF14','CD28','ICOS',#co_stimulatory
                              'FOXP3','CTLA4','IL2RA','IKZF2')
average <- AverageExpression(Tcell, features = cellmarker,group.by = "seurat_clusters")
#转成矩阵
mat <- as.matrix(average[['RNA']])
pheatmap(mat, 
         cluster_cols = T, 
         cluster_rows = F, 
         scale = "row",
         cellwidth = 15, 
         cellheight = 12, 
         border_color = "white",
         color = colorRampPalette(c("#31518a","#fcd6b7", "#820b44"))(50))
DotPlot(Tcell, features =cellmarker, group.by = "seurat_clusters")+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FBE5D6','#F8CBAD','#F4B183','#C55A11'))

(n=length(unique(Tcell@meta.data$seurat_clusters)))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unknown')
celltype[celltype$ClusterID %in% c(4,7),2]='CD4_CD44'
celltype[celltype$ClusterID %in% c(11),2]='CD8_LAG3'
celltype[celltype$ClusterID %in% c(3,14),2]='CD4_CCR7'
celltype[celltype$ClusterID %in% c(0,5),2]='CD4_TIGIT'
celltype[celltype$ClusterID %in% c(1,2),2]='CD8_CXCR4'
celltype[celltype$ClusterID %in% c(6),2]='CD8_LAYN'
celltype[celltype$ClusterID %in% c(8),2]='CD8_LEF1'
celltype[celltype$ClusterID %in% c(12),2]='CD8_GNLY'
celltype[celltype$ClusterID %in% c(13),2]='CD4_Treg'
celltype[celltype$ClusterID %in% c(9),2]='CD8_IFNG'

Tcell@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  Tcell@meta.data[which(Tcell@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

average <- AverageExpression(Tcell, features = cellmarker,group.by = "celltype")
#转成矩阵
mat <- as.matrix(average[['RNA']])
pheatmap(mat, 
         cluster_cols = F, 
         cluster_rows = F, 
         scale = "row",
         cellwidth = 15, 
         cellheight = 12, 
         border_color = "white",
         color = colorRampPalette(c("#31518a","#fcd6b7", "#820b44"))(50))

UMAPPlot(Tcell,group.by = "seurat_clusters")
UMAPPlot(Tcell,group.by = "celltype")

library(monocle)
data <- as(as.matrix(Tcell@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = Tcell@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
#expressionFamily参数用于指定表达矩阵的数据类型，有几个选项可以选择：
#稀疏矩阵用negbinomial.size()，FPKM值用tobit()，logFPKM值用gaussianff()
save(mycds, file = "mycds_raw_cd8.Rdata")

load("mycds_raw.Rdata")
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE) #!
#完成数据导入和预处理后，就可以考虑选择特定基因代表细胞的发育特征
#这里可以选取我们之前挑选的marker gene
save(mycds, file = "mycds_raw_cd8.Rdata")
rm()
load("mycds_raw_cd8.Rdata")
markers.gene <- marker.gene$gene
mycds <- setOrderingFilter(mycds, markers.gene)
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree') #!
#耗时，耗内存
#排序
save(mycds,file = "mycds_reduced_cd8.Rdata")
load("mycds_reduced.Rdata")
cds <- orderCells(mycds)
save(cds,file = "cds_reduced_order_cd8.Rdata")
plot_cell_trajectory(cds,color_by = "Pseudotime")
plot_cell_trajectory(cds,color_by = "State")

Tcell$state <- cds$State
table(Tcell$tissue,Tcell$state)
DF <- data.frame(table(Tcell$tissue,Tcell$state))
library(tidyr)
library(reshape2)
data <- spread(DF,Var1,Freq)
data<- melt(data,value.name =  "Var2" )
colnames(data) <- c("Cluster","Sample","Number")
unique(data$Cluster)
cluster=c("state1","state2","state3","state4","state5")
data$Cluster <- factor(data$Cluster,levels = cluster)
library(RColorBrewer)
sample_color <- c(brewer.pal(12, "Set3"),brewer.pal(9,'BrBG'))
#sample_color <- c("#539ad3","#9b7fb3","#c670a2","#e0fed0")
p <- ggplot(data = data, aes(x = Sample, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+#
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=10, colour = "black"))+
  theme(axis.text.x = element_text(size=10, colour = "black",angle = 90))
p

statemarker <- c('TIGIT','LAG3','MIR155HG','LAYN',#Inhibitory
                 'CD44','IFNG','IL2','LAMTOR3'#cytokines
)

library(ggplot2)
library(pheatmap)

average <- AverageExpression(Tcell, features = statemarker,group.by = "state")
#转成矩阵
mat <- as.matrix(average[['RNA']])
pheatmap(mat, 
         cluster_cols = T, 
         cluster_rows = F, 
         scale = "row",
         cellwidth = 15, 
         cellheight = 12, 
         border_color = "white",
         color = colorRampPalette(c("#31518a","#fcd6b7", "#820b44"))(50))

library(scMetabolism)
Tcell <-sc.metabolism.Seurat(obj = Tcell, method = "Vision", imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway<-c("Pyruvate metabolism", "Oxidative phosphorylation", "Glycolysis / Gluconeogenesis","Glutathione metabolism","Fatty acid degradation","Citrate cycle (TCA cycle)")
DotPlot.metabolism(obj = Tcell, pathway = input.pathway, phenotype = "state", norm = "y")

