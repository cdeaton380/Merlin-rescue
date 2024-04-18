##The code used for analysing the single cell data from CDE241 and 287

setwd("/raleighlab/data1/cde/287_241/")
h_pdx <- readRDS("human/h_pdx.rds")
rm(list=ls())
install.packages("openxlsx", dependencies = TRUE)
install.packages("sigPathway")
install.packages(c("dplyr", "Seurat", "ggplot2", "sctransform", "harmony", "SeuratDisk", "patchwork", "reshape2", "RColorBrewer", "SeuratObject", "plyr", "readr", "sigPathway"), dependencies = TRUE)
?install.packages()
library(openxlsx)
library(HGNChelper)
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(SeuratDisk)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(SeuratObject)
library(plyr)
library(readr)

setwd("/raleighlab/data1/cde/287_241/")

CDE287_2.data <- Read10X(data.dir = "287-2/filtered_feature_bc_matrix")
CDE287_3.data <- Read10X(data.dir = "287-3/filtered_feature_bc_matrix")
CDE287_4.data <- Read10X(data.dir = "287-4/filtered_feature_bc_matrix")
CDE287_5.data <- Read10X(data.dir = "287-5/filtered_feature_bc_matrix")
CDE287_6.data <- Read10X(data.dir = "287-6/filtered_feature_bc_matrix")

CDE241_2.data <- Read10X(data.dir = "241-2/filtered_feature_bc_matrix")
CDE241_3.data <- Read10X(data.dir = "241-3/filtered_feature_bc_matrix")
CDE241_4.data <- Read10X(data.dir = "241-4/filtered_feature_bc_matrix")
CDE241_5.data <- Read10X(data.dir = "241-5/filtered_feature_bc_matrix")
CDE241_6.data <- Read10X(data.dir = "241-6/filtered_feature_bc_matrix")
CDE241_8.data <- Read10X(data.dir = "241-8/filtered_feature_bc_matrix")
CDE241_10.data <- Read10X(data.dir = "241-10/filtered_feature_bc_matrix")

CDE287_2 <- CreateSeuratObject(counts=CDE287_2.data, project = "287-2", min.cells = 3, min.features = 100)
CDE287_3 <- CreateSeuratObject(counts=CDE287_3.data, project = "283-3", min.cells = 3, min.features = 100)
CDE287_4 <- CreateSeuratObject(counts=CDE287_4.data, project = "287-4", min.cells = 3, min.features = 100)
CDE287_5 <- CreateSeuratObject(counts=CDE287_5.data, project = "287-5", min.cells = 3, min.features = 100)
CDE287_6 <- CreateSeuratObject(counts=CDE287_6.data, project = "287-6", min.cells = 3, min.features = 100)

CDE241_2 <- CreateSeuratObject(counts=CDE241_2.data, project = "241-2", min.cells = 3, min.features = 100)
CDE241_3 <- CreateSeuratObject(counts=CDE241_3.data, project = "241-3", min.cells = 3, min.features = 100)
CDE241_4 <- CreateSeuratObject(counts=CDE241_4.data, project = "241-4", min.cells = 3, min.features = 100)
CDE241_5 <- CreateSeuratObject(counts=CDE241_5.data, project = "241-5", min.cells = 3, min.features = 100)
CDE241_6 <- CreateSeuratObject(counts=CDE241_6.data, project = "241-6", min.cells = 3, min.features = 100)
CDE241_8 <- CreateSeuratObject(counts=CDE241_8.data, project = "241-8", min.cells = 3, min.features = 100)
CDE241_10 <- CreateSeuratObject(counts=CDE241_10.data, project = "241-10", min.cells = 3, min.features = 100)

CDE287_2 <- subset(CDE287_2, subset = nFeature_RNA > 100 & nFeature_RNA < 7500)
CDE287_3 <- subset(CDE287_3, subset = nFeature_RNA > 100 & nFeature_RNA < 7500)
CDE287_4 <- subset(CDE287_4, subset = nFeature_RNA > 100 & nFeature_RNA < 7500)
CDE287_5 <- subset(CDE287_5, subset = nFeature_RNA > 100 & nFeature_RNA < 7500)
CDE287_6 <- subset(CDE287_6, subset = nFeature_RNA > 100 & nFeature_RNA < 7500)

CDE241_2 <- subset(CDE241_2, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_3 <- subset(CDE241_3, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_4 <- subset(CDE241_4, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_5 <- subset(CDE241_5, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_6 <- subset(CDE241_6, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_8 <- subset(CDE241_8, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
CDE241_10 <- subset(CDE241_10, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)



rm(CDE287_2.data)
rm(CDE287_3.data)
rm(CDE287_4.data)
rm(CDE287_5.data)
rm(CDE287_6.data)

rm(CDE241_2.data)
rm(CDE241_3.data)
rm(CDE241_4.data)
rm(CDE241_5.data)
rm(CDE241_6.data)
rm(CDE241_8.data)
rm(CDE241_10.data)


##Filter cells by organism
CDE287_2[["percent.mm"]] <- PercentageFeatureSet(CDE287_2, pattern = "^mm10-")
CDE287_3[["percent.mm"]] <- PercentageFeatureSet(CDE287_3, pattern = "^mm10-")
CDE287_4[["percent.mm"]] <- PercentageFeatureSet(CDE287_4, pattern = "^mm10-")
CDE287_5[["percent.mm"]] <- PercentageFeatureSet(CDE287_5, pattern = "^mm10-")
CDE287_6[["percent.mm"]] <- PercentageFeatureSet(CDE287_6, pattern = "^mm10-")

CDE241_2[["percent.mm"]] <- PercentageFeatureSet(CDE241_2, pattern = "^mm10-")
CDE241_3[["percent.mm"]] <- PercentageFeatureSet(CDE241_3, pattern = "^mm10-")
CDE241_4[["percent.mm"]] <- PercentageFeatureSet(CDE241_4, pattern = "^mm10-")
CDE241_5[["percent.mm"]] <- PercentageFeatureSet(CDE241_5, pattern = "^mm10-")
CDE241_6[["percent.mm"]] <- PercentageFeatureSet(CDE241_6, pattern = "^mm10-")
CDE241_8[["percent.mm"]] <- PercentageFeatureSet(CDE241_8, pattern = "^mm10-")
CDE241_10[["percent.mm"]] <- PercentageFeatureSet(CDE241_10, pattern = "^mm10-")


CDE287_2[["percent.hu"]] <- PercentageFeatureSet(CDE287_2, pattern = "^hg19-")
CDE287_3[["percent.hu"]] <- PercentageFeatureSet(CDE287_3, pattern = "^hg19-")
CDE287_4[["percent.hu"]] <- PercentageFeatureSet(CDE287_4, pattern = "^hg19-")
CDE287_5[["percent.hu"]] <- PercentageFeatureSet(CDE287_5, pattern = "^hg19-")
CDE287_6[["percent.hu"]] <- PercentageFeatureSet(CDE287_6, pattern = "^hg19-")

CDE241_2[["percent.hu"]] <- PercentageFeatureSet(CDE241_2, pattern = "^hg19-")
CDE241_3[["percent.hu"]] <- PercentageFeatureSet(CDE241_3, pattern = "^hg19-")
CDE241_4[["percent.hu"]] <- PercentageFeatureSet(CDE241_4, pattern = "^hg19-")
CDE241_5[["percent.hu"]] <- PercentageFeatureSet(CDE241_5, pattern = "^hg19-")
CDE241_6[["percent.hu"]] <- PercentageFeatureSet(CDE241_6, pattern = "^hg19-")
CDE241_8[["percent.hu"]] <- PercentageFeatureSet(CDE241_8, pattern = "^hg19-")
CDE241_10[["percent.hu"]] <- PercentageFeatureSet(CDE241_10, pattern = "^hg19-")



CDE287_2[["human"]] = PercentageFeatureSet(CDE287_2, pattern="^hg19")
CDE287_2[["mouse"]] = PercentageFeatureSet(CDE287_2, pattern="^mm10")
CDE287_3[["human"]] = PercentageFeatureSet(CDE287_3, pattern="^hg19")
CDE287_3[["mouse"]] = PercentageFeatureSet(CDE287_3, pattern="^mm10")
CDE287_4[["human"]] = PercentageFeatureSet(CDE287_4, pattern="^hg19")
CDE287_4[["mouse"]] = PercentageFeatureSet(CDE287_4, pattern="^mm10")
CDE287_5[["human"]] = PercentageFeatureSet(CDE287_5, pattern="^hg19")
CDE287_5[["mouse"]] = PercentageFeatureSet(CDE287_5, pattern="^mm10")
CDE287_6[["human"]] = PercentageFeatureSet(CDE287_6, pattern="^hg19")
CDE287_6[["mouse"]] = PercentageFeatureSet(CDE287_6, pattern="^mm10")

CDE241_2[["human"]] = PercentageFeatureSet(CDE241_2, pattern="^hg19")
CDE241_2[["mouse"]] = PercentageFeatureSet(CDE241_2, pattern="^mm10")
CDE241_3[["human"]] = PercentageFeatureSet(CDE241_3, pattern="^hg19")
CDE241_3[["mouse"]] = PercentageFeatureSet(CDE241_3, pattern="^mm10")
CDE241_4[["human"]] = PercentageFeatureSet(CDE241_4, pattern="^hg19")
CDE241_4[["mouse"]] = PercentageFeatureSet(CDE241_4, pattern="^mm10")
CDE241_5[["human"]] = PercentageFeatureSet(CDE241_5, pattern="^hg19")
CDE241_5[["mouse"]] = PercentageFeatureSet(CDE241_5, pattern="^mm10")
CDE241_6[["human"]] = PercentageFeatureSet(CDE241_6, pattern="^hg19")
CDE241_6[["mouse"]] = PercentageFeatureSet(CDE241_6, pattern="^mm10")
CDE241_8[["human"]] = PercentageFeatureSet(CDE241_8, pattern="^hg19")
CDE241_8[["mouse"]] = PercentageFeatureSet(CDE241_8, pattern="^mm10")
CDE241_10[["human"]] = PercentageFeatureSet(CDE241_10, pattern="^hg19")
CDE241_10[["mouse"]] = PercentageFeatureSet(CDE241_10, pattern="^mm10")

CDE287_2@meta.data$human.umi = (CDE287_2@meta.data$nCount_RNA * CDE287_2@meta.data$human)/100
CDE287_2@meta.data$mouse.umi = (CDE287_2@meta.data$nCount_RNA * CDE287_2@meta.data$mouse)/100
CDE287_3@meta.data$human.umi = (CDE287_3@meta.data$nCount_RNA * CDE287_3@meta.data$human)/100
CDE287_3@meta.data$mouse.umi = (CDE287_3@meta.data$nCount_RNA * CDE287_3@meta.data$mouse)/100
CDE287_4@meta.data$human.umi = (CDE287_4@meta.data$nCount_RNA * CDE287_4@meta.data$human)/100
CDE287_4@meta.data$mouse.umi = (CDE287_4@meta.data$nCount_RNA * CDE287_4@meta.data$mouse)/100
CDE287_5@meta.data$human.umi = (CDE287_5@meta.data$nCount_RNA * CDE287_5@meta.data$human)/100
CDE287_5@meta.data$mouse.umi = (CDE287_5@meta.data$nCount_RNA * CDE287_5@meta.data$mouse)/100
CDE287_6@meta.data$human.umi = (CDE287_6@meta.data$nCount_RNA * CDE287_6@meta.data$human)/100
CDE287_6@meta.data$mouse.umi = (CDE287_6@meta.data$nCount_RNA * CDE287_6@meta.data$mouse)/100

CDE241_2@meta.data$human.umi = (CDE241_2@meta.data$nCount_RNA * CDE241_2@meta.data$human)/100
CDE241_2@meta.data$mouse.umi = (CDE241_2@meta.data$nCount_RNA * CDE241_2@meta.data$mouse)/100
CDE241_3@meta.data$human.umi = (CDE241_3@meta.data$nCount_RNA * CDE241_3@meta.data$human)/100
CDE241_3@meta.data$mouse.umi = (CDE241_3@meta.data$nCount_RNA * CDE241_3@meta.data$mouse)/100
CDE241_4@meta.data$human.umi = (CDE241_4@meta.data$nCount_RNA * CDE241_4@meta.data$human)/100
CDE241_4@meta.data$mouse.umi = (CDE241_4@meta.data$nCount_RNA * CDE241_4@meta.data$mouse)/100
CDE241_5@meta.data$human.umi = (CDE241_5@meta.data$nCount_RNA * CDE241_5@meta.data$human)/100
CDE241_5@meta.data$mouse.umi = (CDE241_5@meta.data$nCount_RNA * CDE241_5@meta.data$mouse)/100
CDE241_6@meta.data$human.umi = (CDE241_6@meta.data$nCount_RNA * CDE241_6@meta.data$human)/100
CDE241_6@meta.data$mouse.umi = (CDE241_6@meta.data$nCount_RNA * CDE241_6@meta.data$mouse)/100
CDE241_8@meta.data$human.umi = (CDE241_8@meta.data$nCount_RNA * CDE241_8@meta.data$human)/100
CDE241_8@meta.data$mouse.umi = (CDE241_8@meta.data$nCount_RNA * CDE241_8@meta.data$mouse)/100
CDE241_10@meta.data$human.umi = (CDE241_10@meta.data$nCount_RNA * CDE241_10@meta.data$human)/100
CDE241_10@meta.data$mouse.umi = (CDE241_10@meta.data$nCount_RNA * CDE241_10@meta.data$mouse)/100


#FeatureScatter(CDE287_2, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE287_3, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE287_4, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE287_5, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE287_6, feature1="human.umi", feature2="mouse.umi")

#FeatureScatter(CDE241_2, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_3, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_4, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_5, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_6, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_8, feature1="human.umi", feature2="mouse.umi")
#FeatureScatter(CDE241_10, feature1="human.umi", feature2="mouse.umi")

CDE287_2@meta.data$cell.type = "doublet"
CDE287_3@meta.data$cell.type = "doublet"
CDE287_4@meta.data$cell.type = "doublet"
CDE287_5@meta.data$cell.type = "doublet"
CDE287_6@meta.data$cell.type = "doublet"

CDE241_2@meta.data$cell.type = "doublet"
CDE241_3@meta.data$cell.type = "doublet"
CDE241_4@meta.data$cell.type = "doublet"
CDE241_5@meta.data$cell.type = "doublet"
CDE241_6@meta.data$cell.type = "doublet"
CDE241_8@meta.data$cell.type = "doublet"
CDE241_10@meta.data$cell.type = "doublet"


CDE287_2@meta.data$cell.type[CDE287_2@meta.data$human > 99] = "human"
CDE287_2@meta.data$cell.type[CDE287_2@meta.data$mouse > 99] = "mouse"
CDE287_3@meta.data$cell.type[CDE287_3@meta.data$human > 99] = "human"
CDE287_3@meta.data$cell.type[CDE287_3@meta.data$mouse > 99] = "mouse"
CDE287_4@meta.data$cell.type[CDE287_4@meta.data$human > 99] = "human"
CDE287_4@meta.data$cell.type[CDE287_4@meta.data$mouse > 99] = "mouse"
CDE287_5@meta.data$cell.type[CDE287_5@meta.data$human > 99] = "human"
CDE287_5@meta.data$cell.type[CDE287_5@meta.data$mouse > 99] = "mouse"
CDE287_6@meta.data$cell.type[CDE287_6@meta.data$human > 99] = "human"
CDE287_6@meta.data$cell.type[CDE287_6@meta.data$mouse > 99] = "mouse"

CDE241_2@meta.data$cell.type[CDE241_2@meta.data$human > 99] = "human"
CDE241_2@meta.data$cell.type[CDE241_2@meta.data$mouse > 99] = "mouse"
CDE241_3@meta.data$cell.type[CDE241_3@meta.data$human > 99] = "human"
CDE241_3@meta.data$cell.type[CDE241_3@meta.data$mouse > 99] = "mouse"
CDE241_4@meta.data$cell.type[CDE241_4@meta.data$human > 99] = "human"
CDE241_4@meta.data$cell.type[CDE241_4@meta.data$mouse > 99] = "mouse"
CDE241_5@meta.data$cell.type[CDE241_5@meta.data$human > 99] = "human"
CDE241_5@meta.data$cell.type[CDE241_5@meta.data$mouse > 99] = "mouse"
CDE241_6@meta.data$cell.type[CDE241_6@meta.data$human > 99] = "human"
CDE241_6@meta.data$cell.type[CDE241_6@meta.data$mouse > 99] = "mouse"
CDE241_8@meta.data$cell.type[CDE241_8@meta.data$human > 99] = "human"
CDE241_8@meta.data$cell.type[CDE241_8@meta.data$mouse > 99] = "mouse"
CDE241_10@meta.data$cell.type[CDE241_10@meta.data$human > 99] = "human"
CDE241_10@meta.data$cell.type[CDE241_10@meta.data$mouse > 99] = "mouse"


hCDE287_2 = subset(CDE287_2, subset = cell.type == "human")
hCDE287_3 = subset(CDE287_3, subset = cell.type == "human")
hCDE287_4 = subset(CDE287_4, subset = cell.type == "human")
hCDE287_5 = subset(CDE287_5, subset = cell.type == "human")
hCDE287_6 = subset(CDE287_6, subset = cell.type == "human")

hCDE241_2 = subset(CDE241_2, subset = cell.type == "human")
hCDE241_3 = subset(CDE241_3, subset = cell.type == "human")
hCDE241_4 = subset(CDE241_4, subset = cell.type == "human")
hCDE241_5 = subset(CDE241_5, subset = cell.type == "human")
hCDE241_6 = subset(CDE241_6, subset = cell.type == "human")
hCDE241_8 = subset(CDE241_8, subset = cell.type == "human")
hCDE241_10 = subset(CDE241_10, subset = cell.type == "human")

hCDE287_2[["percent.mt"]] <- PercentageFeatureSet(hCDE287_2, pattern = "^hg19-MT-")
hCDE287_3[["percent.mt"]] <- PercentageFeatureSet(hCDE287_3, pattern = "^hg19-MT-")
hCDE287_4[["percent.mt"]] <- PercentageFeatureSet(hCDE287_4, pattern = "^hg19-MT-")
hCDE287_5[["percent.mt"]] <- PercentageFeatureSet(hCDE287_5, pattern = "^hg19-MT-")
hCDE287_6[["percent.mt"]] <- PercentageFeatureSet(hCDE287_6, pattern = "^hg19-MT-")

hCDE241_2[["percent.mt"]] <- PercentageFeatureSet(hCDE241_2, pattern = "^hg19-MT-")
hCDE241_3[["percent.mt"]] <- PercentageFeatureSet(hCDE241_3, pattern = "^hg19-MT-")
hCDE241_4[["percent.mt"]] <- PercentageFeatureSet(hCDE241_4, pattern = "^hg19-MT-")
hCDE241_5[["percent.mt"]] <- PercentageFeatureSet(hCDE241_5, pattern = "^hg19-MT-")
hCDE241_6[["percent.mt"]] <- PercentageFeatureSet(hCDE241_6, pattern = "^hg19-MT-")
hCDE241_8[["percent.mt"]] <- PercentageFeatureSet(hCDE241_8, pattern = "^hg19-MT-")
hCDE241_10[["percent.mt"]] <- PercentageFeatureSet(hCDE241_10, pattern = "^hg19-MT-")

hCDE287_2[["percent.nf2"]] <- PercentageFeatureSet(hCDE287_2, pattern = "^hg19-NF2")
hCDE287_3[["percent.nf2"]] <- PercentageFeatureSet(hCDE287_3, pattern = "^hg19-NF2")
hCDE287_4[["percent.nf2"]] <- PercentageFeatureSet(hCDE287_4, pattern = "^hg19-NF2")
hCDE287_5[["percent.nf2"]] <- PercentageFeatureSet(hCDE287_5, pattern = "^hg19-NF2")
hCDE287_6[["percent.nf2"]] <- PercentageFeatureSet(hCDE287_6, pattern = "^hg19-NF2")

hCDE241_2[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_2, pattern = "^hg19-NF2")
hCDE241_3[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_3, pattern = "^hg19-NF2")
hCDE241_4[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_4, pattern = "^hg19-NF2")
hCDE241_5[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_5, pattern = "^hg19-NF2")
hCDE241_6[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_6, pattern = "^hg19-NF2")
hCDE241_8[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_8, pattern = "^hg19-NF2")
hCDE241_10[["percent.nf2"]] <- PercentageFeatureSet(hCDE241_10, pattern = "^hg19-NF2")

# hCDE287_2 = subset(hCDE287_2, subset = percent.mt < 20 & percent.nf2 == 0.0)
# hCDE287_3 = subset(hCDE287_3, subset = percent.mt < 20 & percent.nf2 == 0.0)
# hCDE287_4 = subset(hCDE287_4, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE287_5 = subset(hCDE287_5, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE287_6 = subset(hCDE287_6, subset = percent.mt < 20 & percent.nf2 > 0.0)
# 
# hCDE241_2 = subset(hCDE241_2, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE241_3 = subset(hCDE241_3, subset = percent.mt < 20 & percent.nf2 == 0.0)
# hCDE241_4 = subset(hCDE241_4, subset = percent.mt < 20 & percent.nf2 == 0.0)
# hCDE241_5 = subset(hCDE241_5, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE241_6 = subset(hCDE241_6, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE241_8 = subset(hCDE241_5, subset = percent.mt < 20 & percent.nf2 > 0.0)
# hCDE241_10 = subset(hCDE241_6, subset = percent.mt < 20 & percent.nf2 == 0.)


hCDE287_2[["sample.id"]] <- "hCDE287_2"
hCDE287_3[["sample.id"]] <- "hCDE287_3"
hCDE287_4[["sample.id"]] <- "hCDE287_4"
hCDE287_5[["sample.id"]] <- "hCDE287_5"
hCDE287_6[["sample.id"]] <- "hCDE287_6"

hCDE241_2[["sample.id"]] <- "hCDE241_2"
hCDE241_3[["sample.id"]] <- "hCDE241_3"
hCDE241_4[["sample.id"]] <- "hCDE241_4"
hCDE241_5[["sample.id"]] <- "hCDE241_5"
hCDE241_6[["sample.id"]] <- "hCDE241_6"
hCDE241_8[["sample.id"]] <- "hCDE241_8"
hCDE241_10[["sample.id"]] <- "hCDE241_10"

hCDE287_2[["batch.id"]] <- "first_Veh"
hCDE287_3[["batch.id"]] <- "first_Veh"
hCDE287_4[["batch.id"]] <- "first_Rescue"
hCDE287_5[["batch.id"]] <- "first_Rescue"
hCDE287_6[["batch.id"]] <- "first_Rescue"

hCDE241_2[["batch.id"]] <- "second_Rescue"
hCDE241_3[["batch.id"]] <- "second_Veh"
hCDE241_4[["batch.id"]] <- "second_Veh"
hCDE241_5[["batch.id"]] <- "second_Rescue"
hCDE241_6[["batch.id"]] <- "second_Rescue"
hCDE241_8[["batch.id"]] <- "second_Rescue"
hCDE241_10[["batch.id"]] <- "second_Veh"

hCDE287_2[["condition"]] <- "Veh"
hCDE287_3[["condition"]] <- "Veh"
hCDE287_4[["condition"]] <- "Rescue"
hCDE287_5[["condition"]] <- "Rescue"
hCDE287_6[["condition"]] <- "Rescue"

hCDE241_2[["condition"]] <- "Rescue"
hCDE241_3[["condition"]] <- "Veh"
hCDE241_4[["condition"]] <- "Veh"
hCDE241_5[["condition"]] <- "Rescue"
hCDE241_6[["condition"]] <- "Rescue"
hCDE241_8[["condition"]] <- "Rescue"
hCDE241_10[["condition"]] <- "Veh"

hCDE287_2[["batch"]] <- "first"
hCDE287_3[["batch"]] <- "first"
hCDE287_4[["batch"]] <- "first"
hCDE287_5[["batch"]] <- "first"
hCDE287_6[["batch"]] <- "first"

hCDE241_2[["batch"]] <- "second"
hCDE241_3[["batch"]] <- "second"
hCDE241_4[["batch"]] <- "second"
hCDE241_5[["batch"]] <- "second"
hCDE241_6[["batch"]] <- "second"
hCDE241_8[["batch"]] <- "second"
hCDE241_10[["batch"]] <- "second"

rm(CDE287_2)
rm(CDE287_3)
rm(CDE287_4)
rm(CDE287_5)
rm(CDE287_6)

rm(CDE241_2)
rm(CDE241_3)
rm(CDE241_4)
rm(CDE241_5)
rm(CDE241_6)
rm(CDE241_8)
rm(CDE241_10)

#Merge
h_pdx <- merge(x = hCDE287_2, y = c(hCDE287_3, hCDE287_4, hCDE287_5, hCDE287_6, hCDE241_2, hCDE241_3, hCDE241_4, hCDE241_5, hCDE241_6, hCDE241_8, hCDE241_10), project = "Merlin Rescue")


#SCTransform
h_pdx <- SCTransform(h_pdx, vars.to.regress = 'percent.mm', verbose = T)
h_pdx <- RunPCA(h_pdx, assay = 'SCT', npcs = 30, verbose = T)
h_pdx <- RunUMAP(h_pdx, assay = 'SCT', reduction = "pca", dims = 1:20, verbose = T)
plot1 <- DimPlot(h_pdx, reduction = "umap") + plot_annotation(title = "Merlin Rescue_Human")
plot1

#Normalise
h_pdx <- RunHarmony(h_pdx, group.by.vars = "orig.ident", assay.use = 'SCT')
h_pdx.Embeddings <- Embeddings(h_pdx, 'harmony')
p1 <- DimPlot(object = h_pdx, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = h_pdx, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
p1+p2

## Run Harmony
h_pdx <- h_pdx %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = T) %>% 
  FindNeighbors(reduction = "harmony", k.param = 50, dims = 1:30, verbose = T) %>% 
  FindClusters() %>% 
  identity()
h_pdx <- SetIdent(h_pdx,value = "orig.ident")
DimPlot(h_pdx, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()
h_pdx <- SetIdent(h_pdx, value = "seurat_clusters")
pplot1 <- DimPlot(h_pdx, reduction = 'umap', label.size = 4,repel =,split.by = "condition", label = T)
#ggsave("human/harmony_umap.png", pplot1, width = 16, height = 6, limitsize = F)


#Find markers
hpdx.markers <- FindAllMarkers(h_pdx, only.pos = , min.pct = 0.4, logfc.threshold = 0.25, assay = "SCT")
hpdx_marker_table <- hpdx.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)
head(hpdx.markers, 10)
write.table(hpdx.markers, "human/hpdx_markers.csv", sep = ",")

#n.b. this takes ~5min
saveRDS(h_pdx, file="human/h_pdx.rds")


##Do dot plot
dotplot1 <- DotPlot(h_pdx, features = c("hg19-IER5", "hg19-HSPA6", "hg19-DNAJB1", "hg19-DUSP1", "hg19-KLF6",
                            "hg19-RPL31", "hg19-RPS26", "hg19-RPL21", "hg19-RPS13", "hg19-RPL6", "hg19-RPS27",
                            "hg19-NDRG1", "hg19-SLC2A1", "hg19-ERRFI1", "hg19-BNIP3", "hg19-ANGPTL4",
                            "hg19-RPLP1", "hg19-RPL39", "hg19-RPL32", "hg19-RPS29",
                            "hg19-CTNNB1", "hg19-TCF4", "hg19-PDGFRB", "hg19-WNT5A", "hg19-MMP2", "hg19-LRP1",
                            "hg19-CDK1", "hg19-CCNB2", "hg19-CCNA2", "hg19-FOXM1", "hg19-CCNB1", "hg19-AURKA", "hg19-PLK1",
                            "hg19-HIST1H1B", "hg19-CCNE2", "hg19-BRCA2", "hg19-USP1", "hg19-PCNA", "hg19-MCM4", "hg19-PARP1",
                            "hg19-NPM1", "hg19-HMGB1",
                            "hg19-MT-ATP6", "hg19-MT-CO2", "hg19-MT-ND6", "hg19-MT-ND3",
                            "hg19-IFI6", "hg19-IFIT2", "hg19-IFIT3", "hg19-STAT1", "hg19-MX1","hg19-IFIT1",
                            "hg19-NEAT1", "hg19-NKTR", "hg19-VMP1", "hg19-OGT", "hg19-RMB39",
                            "hg19-IL1B", "hg19-IL8", "hg19-SOD2", "hg19-CSF3", "hg19-NFKBIA",
                            "hg19-RPL10A",
                            "hg19-MT-RNR2", "hg19-MT-RNR1", "hg19-MT-NDL4", 
                            "hg19-HIST1H1E",
                            "hg19-PPP1R15A", "hg19-PIK3R3", "hg19-PPP1R10", "hg19-RHOB",
                            "hg19-KRT81", "hg19-TGM2", "hg19-KRT18", "hg19-THBS1"), 
                    group.by = , cols = ) +
                    ylab("Cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                    theme(panel.background = element_rect(fill = ))
ggsave("human/dotplot1.png", dotplot1, width = 20, height = 6, limitsize = F)

dotplot1 <- DotPlot(h_pdx, features = c("hg19-IER5", "hg19-DUSP1",
                                        "hg19-RPL31", "hg19-RPS26",
                                        "hg19-NDRG1", "hg19-ANGPTL4",
                                        "hg19-RPLP1", "hg19-RPS29",
                                        "hg19-CTNNB1", "hg19-LRP1",
                                        "hg19-FOXM1", "hg19-CCNB1",
                                        "hg19-HIST1H1B","hg19-PCNA",
                                        "hg19-NPM1", "hg19-HMGB1",
                                        "hg19-MT-ATP6", "hg19-MT-CO2", 
                                        "hg19-IFI6","hg19-STAT1", 
                                        "hg19-NKTR", "hg19-RBM39",
                                        "hg19-IL8", "hg19-CSF3", 
                                        "hg19-RPL10A",
                                        "hg19-MT-RNR2", "hg19-MT-RNR1", "hg19-MT-NDL4", 
                                        "hg19-HIST1H1E",
                                        "hg19-PIK3R3", "hg19-PPP1R10",
                                        "hg19-KRT81", "hg19-THBS1"), 
                    group.by = , cols = ) +
  ylab("Cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(panel.background = element_rect(fill = ))
ggsave("human/dotplot.pdf", dotplot1, width = 15, height = 6, limitsize = F)

## Do Heatmap
DimHeatmap(h_pdx, dims = 1, cells = 500, balanced = T, nfeatures = 200)
ggsave("")
top10 <- head(h_pdx@assays[["SCT"]]@var.features)


hpdx.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
hpdx_heatmap <- DoHeatmap(h_pdx, features = top10$gene, assay = , slot = "counts") 
DoHeatmap(h_pdx, features = top10$gene, assay = , slot = ) 

ggplot2::ggsave(filename = "heatmap.tiff", plot = hpdx_heatmap)
ggsave("human/heatmap.png", hpdx_heatmap, width = 8, height = 6, limitsize = F)
top10
h_pdx@SCT@scale.data

##Do cluster distribution
n_cells <- FetchData(h_pdx, assay = "SCT",
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
setwd("/raleighlab/data1/cde/287_241/")
write.table(n_cells, "human/cluster_distribution.csv", sep = ",")

## Remove "hg-19" from beginning of gene name for further analysis
h_pdx@assays[["RNA"]]@data@Dimnames[[1]] <- gsub("hg19-", "", h_pdx@assays[["RNA"]]@data@Dimnames[[1]])

h_pdx@assays[["SCT"]]@var.features <- gsub("hg19-", "", h_pdx@assays[["SCT"]]@var.features)
h_pdx@assays[["SCT"]]@var.features <- gsub("mm10-", "", h_pdx@assays[["SCT"]]@var.features)

#Plotting Wnt graphs
#Violin plot of CTNNB1 expression
CTNNB1vln <- VlnPlot(h_pdx, features = "hg19-CTNNB1", split.by = "condition")
ggsave("CTNNB1-vln.png", CTNNB1vln, width = 20, height = 6, limitsize = F)
CTNNB1vln

##Feature plots
CTNNB1 <- FeaturePlot(h_pdx, features = c("hg19-CTNNB1"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-CTNNB1.png", CTNNB1, width = 8, height = 8, limitsize = F)
SSTR2 <- FeaturePlot(h_pdx, features = c("hg19-SSTR2"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-SSTR2.png", SSTR2, width = 8, height = 8, limitsize = F)
PTGES2 <- FeaturePlot(h_pdx, features = c("hg19-PTGES2"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-PTGES2.png", PTGES2, width = 8, height = 8, limitsize = F)
CDK6 <- FeaturePlot(h_pdx, features = c("hg19-CDK6"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-CDK6.png", CTNNB1, width = 8, height = 8, limitsize = F)
FOXM1 <- FeaturePlot(h_pdx, features = c("hg19-FOXM1"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-FOXM1.png", FOXM1, width = 8, height = 8, limitsize = F)
TOP2A <- FeaturePlot(h_pdx, features = c("hg19-TOP2A"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-TOP2A.png", TOP2A, width = 8, height = 8, limitsize = F)
NF2 <- FeaturePlot(h_pdx, features = c("hg19-NF2"), dims = c(1,2), min.cutoff = , order = T)  
ggsave("human/hg19-NF2.png", NF2, width = 8, height = 8, limitsize = F)
SOD2 <- FeaturePlot(h_pdx, features = c("hg19-SOD2"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-SOD2.png", SOD2, width = 8, height = 8, limitsize = F)
IL8 <- FeaturePlot(h_pdx, features = c("hg19-IL8", ) , dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-IL8.png", IL8, width = 8, height = 8, limitsize = F)
CXCL2 <- FeaturePlot(h_pdx, features = c("hg19-CXCL2"), dims = c(1,2), min.cutoff = , order = T)  
ggsave("human/hg19-CXCL2.png", CXCL2, width = 8, height = 8, limitsize = F)
CXCL5 <- FeaturePlot(h_pdx, features = c("hg19-CXCL5"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-CXCL5.png", CXCL5, width = 8, height = 8, limitsize = F)
CXCL3 <- FeaturePlot(h_pdx, features = c("hg19-CXCL3"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-CXCL3.png", CXCL3, width = 8, height = 8, limitsize = F)
TNFAIP3 <- FeaturePlot(h_pdx, features = c("hg19-TNFAIP3"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-TNFAIP3.png", TNFAIP3, width = 8, height = 8, limitsize = F)
NFKBIZ <- FeaturePlot(h_pdx, features = c("hg19-NFKBIZ"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-NFKBIZ.png", NFKBIZ, width = 8, height = 8, limitsize = F)
IFIT3 <- FeaturePlot(h_pdx, features = c("hg19-IFIT1-3"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IFIT3.png", IFIT3, width = 8, height = 8, limitsize = F)
ISG15 <- FeaturePlot(h_pdx, features = c("hg19-ISG15"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-ISG15.png", ISG15, width = 8, height = 8, limitsize = F)
IFI6 <- FeaturePlot(h_pdx, features = c("hg19-IFI6"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IFI6.png", IFI6, width = 8, height = 8, limitsize = F)
HMGA2 <- FeaturePlot(h_pdx, features = c("hg19-HMGA2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-HMGA2.png", HMGA2, width = 8, height = 8, limitsize = F)
SQLE <- FeaturePlot(h_pdx, features = c("hg19-SQLE"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-SQLE.png", SQLE, width = 8, height = 8, limitsize = F)
S100A4 <- FeaturePlot(h_pdx, features = c("hg19-S100A4"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-S100A4.png", S100A4, width = 8, height = 8, limitsize = F)
TGFB1 <- FeaturePlot(h_pdx, features = c("hg19-TGFB1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-TGFB1.png", TGFB1, width = 8, height = 8, limitsize = F)
FN1 <- FeaturePlot(h_pdx, features = c("hg19-FN1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-FN1.png", FN1, width = 8, height = 8, limitsize = F)
ATF3 <- FeaturePlot(h_pdx, features = c("hg19-ATF3"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-ATF3.png", ATF3, width = 8, height = 8, limitsize = F)
FOS <- FeaturePlot(h_pdx, features = c("hg19-FOS"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-FOS.png", FOS, width = 8, height = 8, limitsize = F)
JUN <- FeaturePlot(h_pdx, features = c("hg19-JUN"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-JUN.png", JUN, width = 8, height = 8, limitsize = F)
IE3A <- FeaturePlot(h_pdx, features = c("hg19-IE3A"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IE3A.png", TNFAIP3, width = 8, height = 8, limitsize = F)
DUSP1 <- FeaturePlot(h_pdx, features = c("hg19-DUSP1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-DUSP1.png", DUSP1, width = 8, height = 8, limitsize = F)
EMP1 <- FeaturePlot(h_pdx, features = c("hg19-EMP1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-EMP1.png", EMP1, width = 8, height = 8, limitsize = F)
MCL1 <- FeaturePlot(h_pdx, features = c("hg19-MCL1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MCL1.png", MCL1, width = 8, height = 8, limitsize = F)
IER3 <- FeaturePlot(h_pdx, features = c("hg19-IER3"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-IER3.png", IER3, width = 8, height = 8, limitsize = F)
CYR61 <- FeaturePlot(h_pdx, features = c("hg19-CYR61"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-CYR61.png", CYR61, width = 8, height = 8, limitsize = F)
GADD45A <- FeaturePlot(h_pdx, features = c("hg19-GADD45A"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-GADD45A.png", GADD45A, width = 8, height = 8, limitsize = F)
MYC <- FeaturePlot(h_pdx, features = c("hg19-MYC"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MYC.png", MYC, width = 8, height = 8, limitsize = F)

## FEATURE PLOTS FOR CLUSTER VALIDATION
DUSP1 <- FeaturePlot(h_pdx, features = c("hg19-DUSP1"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-DUSP1.png", DUSP1, width = 8, height = 8, limitsize = F)
JUN <- FeaturePlot(h_pdx, features = c("hg19-JUN"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-JUN.png", JUN, width = 8, height = 8, limitsize = F)
MYC <- FeaturePlot(h_pdx, features = c("hg19-MYC"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-MYC.png", MYC, width = 8, height = 8, limitsize = F)
FOS <- FeaturePlot(h_pdx, features = c("hg19-FOS"), dims = c(1,2), min.cutoff = , order = T)
ggsave("human/hg19-FOS.png", FOS, width = 8, height = 8, limitsize = F)
RPSL31 <- FeaturePlot(h_pdx, features = c("hg19-RPSL31"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-RPSL31.png", RPSL31, width = 8, height = 8, limitsize = F)
RPS26 <- FeaturePlot(h_pdx, features = c("hg19-RPS26"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-RPS26.png", RPS26, width = 8, height = 8, limitsize = F)
RPL21 <- FeaturePlot(h_pdx, features = c("hg19-RPL21"), dims = c(1,2), min.cutoff = , order = T)  
ggsave("human/hg19-RPL21.png", RPL21, width = 8, height = 8, limitsize = F)
RPS13 <- FeaturePlot(h_pdx, features = c("hg19-RPS13"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-RPS13.png", RPS13, width = 8, height = 8, limitsize = F)
RPL6 <- FeaturePlot(h_pdx, features = c("hg19-RPL6", ) , dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-RPL6.png", RPL6, width = 8, height = 8, limitsize = F)
RPS27 <- FeaturePlot(h_pdx, features = c("hg19-RPS27"), dims = c(1,2), min.cutoff = , order = T)  
ggsave("human/hg19-RPS27.png", RPS27, width = 8, height = 8, limitsize = F)
VEGFA <- FeaturePlot(h_pdx, features = c("hg19-VEGFA"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-VEGFA.png", VEGFA, width = 8, height = 8, limitsize = F)
ADAMTS1 <- FeaturePlot(h_pdx, features = c("hg19-ADAMTS1"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-ADAMTS1.png", ADAMTS1, width = 8, height = 8, limitsize = F)
ADAMTS1.1 <- FeaturePlot(h_pdx, features = c("hg19-ADAMTS1"), dims = c(1,2), min.cutoff = 2, order = T) 
ggsave("human/hg19-ADAMTS1-1.png", ADAMTS1.1, width = 8, height = 8, limitsize = F)
SERPINE1 <- FeaturePlot(h_pdx, features = c("hg19-SERPINE1"), dims = c(1,2), min.cutoff = , order = T) 
ggsave("human/hg19-SERPINE1.png", SERPINE1, width = 8, height = 8, limitsize = F)
PGK1 <- FeaturePlot(h_pdx, features = c("hg19-PGK1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-PGK1.png", PGK1, width = 8, height = 8, limitsize = F)
ANGPTL4 <- FeaturePlot(h_pdx, features = c("hg19-ANGPTL4"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-ANGPTL4.png", ANGPTL4, width = 8, height = 8, limitsize = F)
RPLP1 <- FeaturePlot(h_pdx, features = c("hg19-RPLP1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-RPLP1.png", RPLP1, width = 8, height = 8, limitsize = F)
RPL39 <- FeaturePlot(h_pdx, features = c("hg19-RPL39"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-RPL39.png", RPL39, width = 8, height = 8, limitsize = F)
RPL32 <- FeaturePlot(h_pdx, features = c("hg19-RPL32"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-RPL32.png", RPL32, width = 8, height = 8, limitsize = F)
RPS29 <- FeaturePlot(h_pdx, features = c("hg19-RPS29"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-RPS29.png", RPS29, width = 8, height = 8, limitsize = F)
TCF21 <- FeaturePlot(h_pdx, features = c("hg19-TCF21"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-TCF21.png", TCF21, width = 8, height = 8, limitsize = F)
PDGFRB <- FeaturePlot(h_pdx, features = c("hg19-PDGFRB"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-PDGFRB.png", PDGFRB, width = 8, height = 8, limitsize = F)
WNT5A <- FeaturePlot(h_pdx, features = c("hg19-WNT5A"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-WNT5A.png", WNT5A, width = 8, height = 8, limitsize = F)
MIK67 <- FeaturePlot(h_pdx, features = c("hg19-MIK67"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MIK67.png", MIK67, width = 8, height = 8, limitsize = F)
TOP2A <- FeaturePlot(h_pdx, features = c("hg19-TOP2A"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-TOP2A.png", TOP2A, width = 8, height = 8, limitsize = F)
CKS2 <- FeaturePlot(h_pdx, features = c("hg19-CKS2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-CKS2.png", CKS2, width = 8, height = 8, limitsize = F)
CDK1 <- FeaturePlot(h_pdx, features = c("hg19-CDK1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-CDK1.png", CDK1, width = 8, height = 8, limitsize = F)
CCNB2 <- FeaturePlot(h_pdx, features = c("hg19-CCNB2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-CCNB2.png", CCNB2, width = 8, height = 8, limitsize = F)
CCNA2 <- FeaturePlot(h_pdx, features = c("hg19-CCNA2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-CCNA2.png", CCNA2, width = 8, height = 8, limitsize = F)
FOXM1 <- FeaturePlot(h_pdx, features = c("hg19-FOXM1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-FOXM1.png", FOXM1, width = 8, height = 8, limitsize = F)
HIST1H4C <- FeaturePlot(h_pdx, features = c("hg19-HIST1H4C"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-HIST1H4C.png", HIST1H4C, width = 8, height = 8, limitsize = F)
HIST1H1B <- FeaturePlot(h_pdx, features = c("hg19-HIST1H1B"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-HIST1H1B.png", HIST1H1B, width = 8, height = 8, limitsize = F)
PCNA <- FeaturePlot(h_pdx, features = c("hg19-PCNA"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-PCNA.png", PCNA, width = 8, height = 8, limitsize = F)
CCNE2 <- FeaturePlot(h_pdx, features = c("hg19-CCNE2"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-CCNE2.png", CCNE2, width = 8, height = 8, limitsize = F)

BRCA2 <- FeaturePlot(h_pdx, features = c("hg19-BRCA2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-BRCA2.png", BRCA2, width = 8, height = 8, limitsize = F)
NPM1 <- FeaturePlot(h_pdx, features = c("hg19-NPM1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-NPM1.png", NPM1, width = 8, height = 8, limitsize = F)
HMGB1 <- FeaturePlot(h_pdx, features = c("hg19-HMGB1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-HMGB1.png", HMGB1, width = 8, height = 8, limitsize = F)
MT-ATP6 <- FeaturePlot(h_pdx, features = c("hg19-MT-ATP6"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MT-ATP6.png", MT-ATP6, width = 8, height = 8, limitsize = F)
MT-CO2 <- FeaturePlot(h_pdx, features = c("hg19-MT-CO2"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MT-CO2.png", MT-CO2, width = 8, height = 8, limitsize = F)
MT-ND6 <- FeaturePlot(h_pdx, features = c("hg19-MT-ND6"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MT-ND6.png", MT-ND6, width = 8, height = 8, limitsize = F)
MT-ND3 <- FeaturePlot(h_pdx, features = c("hg19-MT-ND3"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MT-ND3.png", MT-ND3, width = 8, height = 8, limitsize = F)
IFI6 <- FeaturePlot(h_pdx, features = c("hg19-IFI6"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-IFI6.png", IFI6, width = 8, height = 8, limitsize = F)

IFIT2 <- FeaturePlot(h_pdx, features = c("hg19-IFIT2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IFIT2.png", IFIT2, width = 8, height = 8, limitsize = F)
IFI27 <- FeaturePlot(h_pdx, features = c("hg19-IFI27"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IFI27.png", IFI27, width = 8, height = 8, limitsize = F)
IFIT3 <- FeaturePlot(h_pdx, features = c("hg19-IFIT3"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IFIT3.png", IFIT3, width = 8, height = 8, limitsize = F)
STAT1 <- FeaturePlot(h_pdx, features = c("hg19-STAT1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-STAT1.png", STAT1, width = 8, height = 8, limitsize = F)
LY6E <- FeaturePlot(h_pdx, features = c("hg19-LY6E"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-LY6E.png", LY6E, width = 8, height = 8, limitsize = F)
MALAT1 <- FeaturePlot(h_pdx, features = c("hg19-MALAT1"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MALAT1.png", MALAT1, width = 8, height = 8, limitsize = F)
MEG3 <- FeaturePlot(h_pdx, features = c("hg19-MEG3"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-MEG3.png", MEG3, width = 8, height = 8, limitsize = F)
RBM39 <- FeaturePlot(h_pdx, features = c("hg19-RBM39"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-RBM39.png", RBM39, width = 8, height = 8, limitsize = F)

RGBM5 <- FeaturePlot(h_pdx, features = c("hg19-RGBM5"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-RGBM5.png", RGBM5, width = 8, height = 8, limitsize = F)
IL1 <- FeaturePlot(h_pdx, features = c("hg19-IL1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IL1.png", IL1, width = 8, height = 8, limitsize = F)
IL8 <- FeaturePlot(h_pdx, features = c("hg19-IL8"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-IL8.png", IL8, width = 8, height = 8, limitsize = F)
CXCL1 <- FeaturePlot(h_pdx, features = c("hg19-CXCL1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-CXCL1.png", CXCL1, width = 8, height = 8, limitsize = F)
CXCL2 <- FeaturePlot(h_pdx, features = c("hg19-CXCL2"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-CXCL2.png", CXCL2, width = 8, height = 8, limitsize = F)
CXCL3 <- FeaturePlot(h_pdx, features = c("hg19-CXCL3"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-CXCL3.png", CXCL3, width = 8, height = 8, limitsize = F)
RPL6 <- FeaturePlot(h_pdx, features = c("hg19-RPL6"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-RPL6.png", RPL6, width = 8, height = 8, limitsize = F)
RPL10A <- FeaturePlot(h_pdx, features = c("hg19-RPL10A"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-RPL10A.png", RPL10A, width = 8, height = 8, limitsize = F)

MT-RNR2 <- FeaturePlot(h_pdx, features = c("hg19-MT-RNR2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MT-RNR2.png", MT-RNR2, width = 8, height = 8, limitsize = F)
MT-RNR1 <- FeaturePlot(h_pdx, features = c("hg19-MT-RNR1"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MT-RNR1.png", MT-RNR1, width = 8, height = 8, limitsize = F)
MT-ND4L <- FeaturePlot(h_pdx, features = c("hg19-MT-ND4L"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MT-ND4L.png", MT-ND4L, width = 8, height = 8, limitsize = F)
MT-ND2 <- FeaturePlot(h_pdx, features = c("hg19-MT-ND2"), dims = c(1,2), min.cutoff = 0, order = T) 
ggsave("human/hg19-MT-ND2.png", MT-ND2, width = 8, height = 8, limitsize = F)
HIST1H1E <- FeaturePlot(h_pdx, features = c("hg19-HIST1H1E"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-HIST1H1E.png", HIST1H1E, width = 8, height = 8, limitsize = F)
JUND <- FeaturePlot(h_pdx, features = c("hg19-JUND"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-JUND.png", JUND, width = 8, height = 8, limitsize = F)
PIK3R3 <- FeaturePlot(h_pdx, features = c("hg19-PIK3R3"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-PIK3R3.png", PIK3R3, width = 8, height = 8, limitsize = F)
FOSB <- FeaturePlot(h_pdx, features = c("hg19-FOSB"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-FOSB.png", FOSB, width = 8, height = 8, limitsize = F)

JUNB <- FeaturePlot(h_pdx, features = c("hg19-JUNB"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-JUNB.png", JUNB, width = 8, height = 8, limitsize = F)
COL4A1 <- FeaturePlot(h_pdx, features = c("hg19-COL4A1"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-COL4A1.png", COL4A1, width = 8, height = 8, limitsize = F)
COL4A2 <- FeaturePlot(h_pdx, features = c("hg19-COL4A2"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-COL4A2.png", COL4A2, width = 8, height = 8, limitsize = F)
THBS1 <- FeaturePlot(h_pdx, features = c("hg19-THBS1"), dims = c(1,2), min.cutoff = 0, order = T)
ggsave("human/hg19-THBS1.png", THBS1, width = 8, height = 8, limitsize = F)

FeaturePlot(h_pdx, features = c("hg19-OCI5"), dims = c(1,2), min.cutoff = 0, order = T)

plot101 <- FeaturePlot(h_pdx, features = c("hg19-CTNNB1", "hg19-NF2"), order = T, blend = T, cols = c("orange", "brown"))
ggsave("human/NF2vsCTNNB1.png", plot101, width = 8, height = 8, limitsize = F)

plot102 <- FeatureScatter(h_pdx, feature1 = "hg19-CTNNB1", feature2 = "hg19-NF2", slot = "data")
ggsave("human/NF2vsCTNNB1_2.png", plot102, width = 8, height = 8, limitsize = F)
plot103 <- FeatureScatter(h_pdx, feature1 = "hg19-TP53", feature2 = "hg19-NF2", slot = "scale.data")
ggsave("human/NF2vsCTNNB1_2.png", plot102, width = 8, height = 8, limitsize = F)

plot104 <- VlnPlot(h_pdx, features = c("hg19-CTNNB1", "hg19-NF2"))
plot105 <- DotPlot(h_pdx, features = "hg19-NF2", split.by = "condition")
plot106 <- DotPlot(h_pdx, features = "hg19-CTNNB1", split.by = "condition")
ggsave("human/NF2_dot.pdf", plot105, width = 8, height = 8, limitsize = F)
ggsave("human/CTNNB1_dot.pdf", plot106, width = 8, height = 8, limitsize = F)


plot#Geneset plots
setwd("/raleighlab/data1/cde/genelists")
Wntgenes <- read.delim(file = "GOBP_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY.v7.5.1.grp", sep = "")
h_pdx <- AddModuleScore(h_pdx, features = Wntgenes, assay = 'SCT', ctrl = 3, name = "Wntscore", search = F)
Wntplot2 <- FeaturePlot(h_pdx, features = 'Wntscore1', order = T, split.by = "condition", min.cutoff = 0, pt.size = .3)
setwd("/raleighlab/data1/cde/287_241/")
ggsave("human/Wntplot.png", Wntplot2, width = 16, height = 8, limitsize = F)
Wntplot2

## More Wnt plots
WntNF2plot <- FeaturePlot(h_pdx, feature = c('hg19-NF2', "Wntscore1"), order = T, split.by = 'condition', pt.size = 0.1, min.cutoff = 0)
WntNF2plot
WntBplot <- FeaturePlot(h_pdx, feature = c('hg19-CTNNB1', "Wntscore1"), order = T, split.by = 'condition', pt.size = 0.1, min.cutoff = 0)
WntBplot
pplot1

## Find out total cells in each cluster per condition
n_cells <- FetchData(h_pdx, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
setwd("/raleighlab/data1/cde/287_241/")
write.table(n_cells, "human/cluster_distribution.csv", sep = ",")


## Apical gene set enrichment
setwd("/raleighlab/data1/cde/genelists")
Apicalgenes <- read.delim(file = 'HALLMARK_APICAL_JUNCTION.v7.5.1.grp', sep = '')
h_pdx <- AddModuleScore(h_pdx, features = Apicalgenes, assay = 'RNA', ctrl = 3, slot = data, name = "Apicalscore")
Apicalplot <- FeaturePlot(h_pdx, features = 'Apicalscore1', order = T, split.by = 'orig.ident', min.cutoff = 0)
Apicalplot

##Cell cycle analysis
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hpdx_cc <- CellCycleScoring(h_pdx, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = 'RNA')
head(hpdx_cc[[]])
hpdx_cc <- RunPCA(hpdx_cc, features = c(s.genes, g2m.genes))
plot40 <- DimPlot(hpdx_cc)
plot40
setwd("/raleighlab/data1/cde/287_241/")
ggsave("human/cell_cycle_clusters.png", plot40, width = 8, height = 6, limitsize = F)

##UMAP coordinates
head(h_pdx[["umap"]]@cell.embeddings)
hpdx_umap_coords_1 <- h_pdx[["umap"]]@cell.embeddings
hpdx_umap_coords_2 <- h_pdx@meta.data$seurat_clusters
hpdx_umap_coords <- cbind(hpdx_umap_coords_1, hpdx_umap_coords_2)
hpdx_umap_coords_3 <- hpdx_cc@meta.data$Phase
hpdx_umap_coords <- cbind(hpdx_umap_coords, hpdx_umap_coords_3)
setwd("/raleighlab/data1/cde/287_241/")
write.table(hpdx_umap_coords, 'human/hpdx_umap_coords.csv', sep=',')
head(hpdx_umap_coords)


##scType analysis
install.packages("devtools", dependencies = TRUE)
install.packages("openxlsx", dependencies = TRUE)
install.packages("HGNChelper", dependencies = TRUE)
library(HGNChelper)
library(openxlsx)
library(devtools)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = hpdx[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(hpdx@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(hpdx@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(hpdx@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

hpdx@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  hpdx@meta.data$customclassif[hpdx@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(hpdx, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  


