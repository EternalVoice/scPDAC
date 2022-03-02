# CRA001160
library(Seurat)
# Normal
N1 <- Read10X("data/PDAC/N1/")
N1 <- CreateSeuratObject(counts = N1, min.cells = 3)
N1$SampleName <- "N1"
N1$Group <- "Healthy"
N1$Patient <- "Normal"
N1$Age <- 64
N1$Gender <- "Female"
N1$GEO <- "CRA001160"

N2 <- Read10X("data/PDAC/N2/")
N2 <- CreateSeuratObject(counts = N2, min.cells = 3)
N2$SampleName <- "N2"
N2$Group <- "Healthy"
N2$Patient <- "Normal"
N2$Age <- 55
N2$Gender <- "Male"
N2$GEO <- "CRA001160"

N3 <- Read10X("data/PDAC/N3/")
N3 <- CreateSeuratObject(counts = N3, min.cells = 3)
N3$SampleName <- "N3"
N3$Group <- "Healthy"
N3$Patient <- "Normal"
N3$Age <- 50
N3$Gender <- "Male"
N3$GEO <- "CRA001160"

N4 <- Read10X("data/PDAC/N4/")
N4 <- CreateSeuratObject(counts = N4, min.cells = 3)
N4$SampleName <- "N4"
N4$Group <- "Healthy"
N4$Patient <- "Normal"
N4$Age <- 53
N4$Gender <- "Male"
N4$GEO <- "CRA001160"

N5 <- Read10X("data/PDAC/N5/")
N5 <- CreateSeuratObject(counts = N5, min.cells = 3)
N5$SampleName <- "N5"
N5$Group <- "Healthy"
N5$Patient <- "Normal"
N5$Age <- 52
N5$Gender <- "Female"
N5$GEO <- "CRA001160"

N6 <- Read10X("data/PDAC/N6/")
N6 <- CreateSeuratObject(counts = N6, min.cells = 3)
N6$SampleName <- "N6"
N6$Group <- "Healthy"
N6$Patient <- "Normal"
N6$Age <- 31
N6$Gender <- "Female"
N6$GEO <- "CRA001160"

N7 <- Read10X("data/PDAC/N7/")
N7 <- CreateSeuratObject(counts = N7, min.cells = 3)
N7$SampleName <- "N7"
N7$Group <- "Healthy"
N7$Patient <- "Normal"
N7$Age <- 42
N7$Gender <- "Female"
N7$GEO <- "CRA001160"

N8 <- Read10X("data/PDAC/N8/")
N8 <- CreateSeuratObject(counts = N8, min.cells = 3)
N8$SampleName <- "N8"
N8$Group <- "Healthy"
N8$Patient <- "Normal"
N8$Age <- 41
N8$Gender <- "Male"
N8$GEO <- "CRA001160"

N9 <- Read10X("data/PDAC/N9/")
N9 <- CreateSeuratObject(counts = N9, min.cells = 3)
N9$SampleName <- "N9"
N9$Group <- "Healthy"
N9$Patient <- "Normal"
N9$Age <- 34
N9$Gender <- "Male"
N9$GEO <- "CRA001160"

N10 <- Read10X("data/PDAC/N10/")
N10 <- CreateSeuratObject(counts = N10, min.cells = 3)
N10$SampleName <- "N10"
N10$Group <- "Healthy"
N10$Patient <- "Normal"
N10$Age <- 65
N10$Gender <- "Female"
N10$GEO <- "CRA001160"

N11 <- Read10X("data/PDAC/N11/")
N11 <- CreateSeuratObject(counts = N11, min.cells = 3)
N11$SampleName <- "N11"
N11$Group <- "Healthy"
N11$Patient <- "Normal"
N11$Age <- 30
N11$Gender <- "Female"
N11$GEO <- "CRA001160"

# Tumor
T1 <- Read10X("data/PDAC/T1/")
T1 <- CreateSeuratObject(counts = T1, min.cells = 3)
T1$SampleName <- "T1"
T1$Group <- "PDACIII-Primary"
T1$Patient <- "PDAC"
T1$Age <- 64
T1$Gender <- "Male"
T1$GEO <- "CRA001160"

T2 <- Read10X("data/PDAC/T2/")
T2 <- CreateSeuratObject(counts = T2, min.cells = 3)
T2$SampleName <- "T2"
T2$Group <- "PDACIIB-Primary"
T2$Patient <- "PDAC"
T2$Age <- 52
T2$Gender <- "Male"
T2$GEO <- "CRA001160"

T3 <- Read10X("data/PDAC/T3/")
T3 <- CreateSeuratObject(counts = T3, min.cells = 3)
T3$SampleName <- "T3"
T3$Group <- "PDACIB-Primary"
T3$Patient <- "PDAC"
T3$Age <- 58
T3$Gender <- "Female"
T3$GEO <- "CRA001160"

T4 <- Read10X("data/PDAC/T4/")
T4 <- CreateSeuratObject(counts = T4, min.cells = 3)
T4$SampleName <- "T4"
T4$Group <- "PDACIIB-Primary"
T4$Patient <- "PDAC"
T4$Age <- 72
T4$Gender <- "Female"
T4$GEO <- "CRA001160"

T5 <- Read10X("data/PDAC/T5/")
T5 <- CreateSeuratObject(counts = T5, min.cells = 3)
T5$SampleName <- "T5"
T5$Group <- "PDACIB-Primary"
T5$Patient <- "PDAC"
T5$Age <- 65
T5$Gender <- "Female"
T5$GEO <- "CRA001160"

T6 <- Read10X("data/PDAC/T6/")
T6 <- CreateSeuratObject(counts = T6, min.cells = 3)
T6$SampleName <- "T6"
T6$Group <- "PDACIIA-Primary"
T6$Patient <- "PDAC"
T6$Age <- 64
T6$Gender <- "Male"
T6$GEO <- "CRA001160"

T7 <- Read10X("data/PDAC/T7/")
T7 <- CreateSeuratObject(counts = T7, min.cells = 3)
T7$SampleName <- "T7"
T7$Group <- "PDACIIB-Primary"
T7$Patient <- "PDAC"
T7$Age <- 70
T7$Gender <- "Male"
T7$GEO <- "CRA001160"

T8 <- Read10X("data/PDAC/T8/")
T8 <- CreateSeuratObject(counts = T8, min.cells = 3)
T8$SampleName <- "T8"
T8$Group <- "PDACIII-Primary"
T8$Patient <- "PDAC"
T8$Age <- 66
T8$Gender <- "Female"
T8$GEO <- "CRA001160"

T9 <- Read10X("data/PDAC/T9/")
T9 <- CreateSeuratObject(counts = T9, min.cells = 3)
T9$SampleName <- "T9"
T9$Group <- "PDACIIA-Primary"
T9$Patient <- "PDAC"
T9$Age <- 36
T9$Gender <- "Male"
T9$GEO <- "CRA001160"

T10 <- Read10X("data/PDAC/T10/")
T10 <- CreateSeuratObject(counts = T10, min.cells = 3)
T10$SampleName <- "T10"
T10$Group <- "PDACIB-Primary"
T10$Patient <- "PDAC"
T10$Age <- 61
T10$Gender <- "Male"
T10$GEO <- "CRA001160"

T11 <- Read10X("data/PDAC/T11/")
T11 <- CreateSeuratObject(counts = T11, min.cells = 3)
T11$SampleName <- "T11"
T11$Group <- "PDACIIB-Primary"
T11$Patient <- "PDAC"
T11$Age <- 51
T11$Gender <- "Male"
T11$GEO <- "CRA001160"

T12 <- Read10X("data/PDAC/T12/")
T12 <- CreateSeuratObject(counts = T12, min.cells = 3)
T12$SampleName <- "T12"
T12$Group <- "PDACIII-Primary"
T12$Patient <- "PDAC"
T12$Age <- 54
T12$Gender <- "Male"
T12$GEO <- "CRA001160"

T13 <- Read10X("data/PDAC/T13/")
T13 <- CreateSeuratObject(counts = T13, min.cells = 3)
T13$SampleName <- "T13"
T13$Group <- "PDACIII-Primary"
T13$Patient <- "PDAC"
T13$Age <- 58
T13$Gender <- "Female"
T13$GEO <- "CRA001160"

T14 <- Read10X("data/PDAC/T14/")
T14 <- CreateSeuratObject(counts = T14, min.cells = 3)
T14$SampleName <- "T14"
T14$Group <- "PDACIIB-Primary"
T14$Patient <- "PDAC"
T14$Age <- 67
T14$Gender <- "Female"
T14$GEO <- "CRA001160"

T15 <- Read10X("data/PDAC/T15/")
T15 <- CreateSeuratObject(counts = T15, min.cells = 3)
T15$SampleName <- "T15"
T15$Group <- "PDACIIB-Primary"
T15$Patient <- "PDAC"
T15$Age <- 54
T15$Gender <- "Female"
T15$GEO <- "CRA001160"

T16 <- Read10X("data/PDAC/T16/")
T16 <- CreateSeuratObject(counts = T16, min.cells = 3)
T16$SampleName <- "T16"
T16$Group <- "PDACIIB-Primary"
T16$Patient <- "PDAC"
T16$Age <- 56
T16$Gender <- "Female"
T16$GEO <- "CRA001160"

T17 <- Read10X("data/PDAC/T17/")
T17 <- CreateSeuratObject(counts = T17, min.cells = 3)
T17$SampleName <- "T17"
T17$Group <- "PDACIB-Primary"
T17$Patient <- "PDAC"
T17$Age <- 71
T17$Gender <- "Female"
T17$GEO <- "CRA001160"

T18 <- Read10X("data/PDAC/T18/")
T18 <- CreateSeuratObject(counts = T18, min.cells = 3)
T18$SampleName <- "T18"
T18$Group <- "PDACIB-Primary"
T18$Patient <- "PDAC"
T18$Age <- 68
T18$Gender <- "Female"
T18$GEO <- "CRA001160"

T19 <- Read10X("data/PDAC/T19/")
T19 <- CreateSeuratObject(counts = T19, min.cells = 3)
T19$SampleName <- "T19"
T19$Group <- "PDACIB-Primary"
T19$Patient <- "PDAC"
T19$Age <- 59
T19$Gender <- "Female"
T19$GEO <- "CRA001160"

T20 <- Read10X("data/PDAC/T20/")
T20 <- CreateSeuratObject(counts = T20, min.cells = 3)
T20$SampleName <- "T20"
T20$Group <- "PDACIIB-Primary"
T20$Patient <- "PDAC"
T20$Age <- 59
T20$Gender <- "Male"
T20$GEO <- "CRA001160"

T21 <- Read10X("data/PDAC/T21/")
T21 <- CreateSeuratObject(counts = T21, min.cells = 3)
T21$SampleName <- "T21"
T21$Group <- "PDACIB-Primary"
T21$Patient <- "PDAC"
T21$Age <- 59
T21$Gender <- "Male"
T21$GEO <- "CRA001160"

T22 <- Read10X("data/PDAC/T22/")
T22 <- CreateSeuratObject(counts = T22, min.cells = 3)
T22$SampleName <- "T22"
T22$Group <- "PDACIB-Primary"
T22$Patient <- "PDAC"
T22$Age <- 67
T22$Gender <- "Female"
T22$GEO <- "CRA001160"

T23 <- Read10X("data/PDAC/T23/")
T23 <- CreateSeuratObject(counts = T23, min.cells = 3)
T23$SampleName <- "T23"
T23$Group <- "PDACIIB-Primary"
T23$Patient <- "PDAC"
T23$Age <- 54
T23$Gender <- "Male"
T23$GEO <- "CRA001160"

T24 <- Read10X("data/PDAC/T24/")
T24 <- CreateSeuratObject(counts = T24, min.cells = 3)
T24$SampleName <- "T24"
T24$Group <- "PDACIB-Primary"
T24$Patient <- "PDAC"
T24$Age <- 44
T24$Gender <- "Female"
T24$GEO <- "CRA001160"

seurat.list <- list(
  N1 = N1,N2 = N2,N3 = N3,N4 = N4,N5 = N5,N6 = N6,N7 = N7,N8 = N8,
  N9 = N9,N10 = N10,N11 = N11,T1 = T1,T2 = T2,T3 = T3,T4 = T4,T5 = T5,
  T6 = T6,T7 = T7,T8 = T8,T9 = T9,T10 = T10,T11 = T11,T12 = T12,T13 = T13,
  T14 = T14,T15 = T15,T16 = T16,T17 = T17,T18 = T18,T19 = T19,T20 = T20,
  T21 = T21,T22 = T22,T23 = T23,T24 = T24
)
saveRDS(seurat.list,file = "raw.seurat.list.rds")

seurat.comb <- merge(seurat.list[[1]],seurat.list[2:length(seurat.list)])
saveRDS(seurat.comb, file = "seurat.comb.rds")

# seurat.list <- lapply(seurat.list, function(x){
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x,selection.method="vst",nfeatures=2000)
# })
# saveRDS(seurat.list, file = "rds/seurat.list.rds")
# 
# scRNA.anchors <- FindIntegrationAnchors(object.list = seurat.list, verbose = F)
# scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors,verbose = F)
# saveRDS(scRNA.anchors, "CRA001160.anchors.rds")
# saveRDS(scRNA.integrated, "CRA001160.integrated.rds")

# --------------------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(rlang)

if(!file.exists("validation")) dir.create("validation")

# preprocessing
seurat.comb <- readRDS("validation/seurat.comb.rds")
Idents(seurat.comb) <- seurat.comb$SampleName
seurat.comb <- NormalizeData(seurat.comb, normalization.method = "LogNormalize",scale.factor = 10000)
seurat.comb <- FindVariableFeatures(seurat.comb,selection.method = "vst")
#Percent Mitochondrial Genes
#QC Metric used to remove cells with overabundant Mitochondrial genes, typically associated with nuclear wash out during sequencing
seurat.comb[["percent.mt"]] <- PercentageFeatureSet(seurat.comb, pattern = "^MT-")
p1 <- VlnPlot(seurat.comb, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 1e-5)
ggsave(paste0("validation/","metrics.beforeQC.pdf"),p1,width = 20,height = 5)

## Quality Control
# QC criteria: 
#     i). remove cells had either lower than 200 or higher than 5000 expressed genes
#    ii). discarded cells with more than 30,000 UMIs
#   iii). mitochondria content higher than 30%
seurat.comb <- subset(seurat.comb, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 30)
p2 <- VlnPlot(seurat.comb, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 0)
ggsave(paste0("validation/","metrics.afterQC.pdf"),p1,width = 20,height = 5)

# Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat.comb <- CellCycleScoring(seurat.comb,s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
seurat.comb$CC.Difference <- seurat.comb$S.Score - seurat.comb$G2M.Score
seurat.comb <- ScaleData(object = seurat.comb, vars.to.regress = "CC.Difference", features = VariableFeatures(seurat.comb))
## Dimensional reduction
# Run PCA and Determine Dimensions for 90% Variance
seurat.comb <- RunPCA(object = seurat.comb, features = VariableFeatures(seurat.comb))

PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}

PC.num <- PCDeterminators(seurat.comb)
# Find Neighbors + Find Clusters
seurat.comb <- FindNeighbors(seurat.comb, dims = 1:PC.num)
seurat.comb <- FindClusters(seurat.comb, resolution = 0.2)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
seurat.comb <- RunUMAP(seurat.comb, dims = 1:PC.num)
seurat.comb <- RunTSNE(seurat.comb,dims = 1:PC.num)

mypal <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
           "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
           "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
           "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
           "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
           "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
           "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
           "#e598d7", "#343b20", "#bbb5d9", "#975223")
mypal2 <- c("#c14089","#6f5553","#E5C494","#738f4c","#bb6240","#66C2A5","#2dfd29","#0c0fdc")
mycol <- c("grey88","DarkCyan")


# Find neighbors and clusters WITH harmony batch correction
Idents(seurat.comb) <- seurat.comb$SampleName
seurat.comb <- RunHarmony(seurat.comb, group.by.vars = "SampleName")
seurat.comb.harmony <- FindNeighbors(seurat.comb,dims = 1:PC.num, reduction = "harmony") %>%
  FindClusters(resolution = 0.2, reduction = "harmony")

saveRDS(seurat.comb.harmony, file = "validation/seurat.comb.harmony.rds")

p.clust <- DimPlot(seurat.comb.harmony,reduction = "tsne",group.by = "seurat_clusters")
ggsave("validation/seurat.clusters.pdf",p.clust,width = 6,height = 5)
p.sn <- DimPlot(seurat.comb.harmony,reduction = "tsne",group.by = "SampleName")

# T_cells <- c("CD3D","CD3E","CD4","CD8A","CD3G","IL7R","LEF1")
# B_cell <- c("MS4A1","CD79A","CD79B","CD52","CD19","SDC1","IGJ","IGLL5","CXCR4","KIT","CD27","HLA-DRA")
# Epithelial_cell <- c("EPCAM","ACTA2","KRT7","KRT8","KRT18","KRT19","CDH1","PRSS1","CTRB2",
#                      "REG1A","CLU","MKI67","SPINK1","TFF1","MUC1")
# ductal_cell <- c("AMBP","CFTR","MMP7","KRT19","KRT7","TSPAN8","SLPI")
# Acinar_cell <- c("PRSS1","CTRB1","CTRB2","REG1B","SPINK1","AMY2A")
# Stellate_cell <- c("RGS5","ACTA2","PDGFRB","ADIRF")
# Endothelial_cell <- c("CDH5","PLVAP","VWF","CLDN5","KDR","PECAM1")
# Fibroblast_cell <- c("LUM","DCN","COL1A1","ACTA2","SPARC","CDH11","PDGFRA","PDGFRB","COL3A1",
#                      "RGS5","IGFBP7","PDPN","MCAM","IL6","APOE","GLI1","GLI2","GLI3","PDGFA")
# NK <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38")

# m1.fea <- c("HLA-DRA","HLA-DRB5","HLA-DRB1","ITGAX","CD86","NOS2","STAT1")
# m2.fea <- c("CD40","SPP1","MARCO","APOE","CHIT1","FABP5","CCL18","HLA-DQA1","COX17","LY6E","LAMP1","HAVCR2","SERPINB6","IL18","CCL2","ATF5","CXCL3","VEGFB","SLC2A1","CD163","MSR1","VEGFA","MAF")
# Mono_CD14: FCN1, S100A8, S100A9
# Mono_CD16: FCGR3A, LST1, LILRB2
# Macrophage:
#     Macro_C1QC: C1QC, C1QA, APOE
#     Macro_INHBA: INHBA, IL1RN, CCL4
#     Macro_NLRP3: NLRP3, EREG, IL1B
#     Macro_LYVE1: LYVE1, PLTP, SEPP1

library(MySeuratWrappers)
Myeloid_cell <- c("APOE","C1QA","LYZ","HLA-DRA")
p.mye <- VlnPlot(seurat.comb.harmony,features = Myeloid_cell,stacked = T,group.by = "seurat_clusters",
                 pt.size = 0,combine = T,direction = "horizontal",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("validation/mye.fea.pdf",p.mye,width = 4,height = 6)

curr.ids <- c(0:15)
new.ids <- c("Unknown","Unknown","Unknown","Unknown",
             "myeloid","Unknown","Unknown","Unknown",
             "Unknown","Unknown","myeloid","myeloid",
             "myeloid","myeloid","myeloid","myeloid")
seurat.comb.harmony$annot <- plyr::mapvalues(seurat.comb.harmony$seurat_clusters,from = curr.ids,to = new.ids)

mye <- subset(seurat.comb.harmony,cells = rownames(subset(seurat.comb.harmony@meta.data,annot == "myeloid")))

saveRDS(mye,file = "validation/mye.rds")

rm(list = ls());gc()

mye <- readRDS("validation/mye.rds")

mye <- NormalizeData(mye, normalization.method = "LogNormalize",scale.factor = 10000)
mye <- FindVariableFeatures(mye,selection.method = "vst")
mye <- ScaleData(object = mye, features = VariableFeatures(mye))
mye <- RunPCA(object = mye, features = VariableFeatures(mye))
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(mye)
# Find Neighbors + Find Clusters
mye <- FindNeighbors(mye, dims = 1:PC.num)
mye <- FindClusters(mye, resolution = 0.2)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
mye <- RunUMAP(mye, dims = 1:PC.num)
mye <- RunTSNE(mye,dims = 1:PC.num)

p.mye <- DimPlot(mye,reduction = "umap",group.by = "seurat_clusters") + 
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/mye.pdf",p.mye,width = 4.3,height = 4)

Macrophage_cell <- c("AIF1","FCGR1A","CD14","CD68")

p.macro <- VlnPlot(mye,features = Macrophage_cell,stacked = T,
                   pt.size = 0,combine = T,direction = "horizontal",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("validation/macro.fea.pdf",p.macro,width = 4,height = 6)


CD16mono <- c("FCGR3A","LST1","LILRB2")
p.CD16mono <- VlnPlot(mye,features = CD16mono,stacked = T,
                      pt.size = 0,combine = T,direction = "horizontal",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("validation/CD16mono.fea.pdf",p.CD16mono,width = 3,height = 6)

mast <- c("CLU")
p.mast <- VlnPlot(mye,features = mast,stacked = T,pt.size = 0,combine = T,direction = "horizontal",
                  x.lab = '',y.lab = '') + NoLegend() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(size = 1,colour = "black"))
ggsave("validation/mast.fea.pdf",p.mast,width = 2,height = 6)


curr.ids <- c(0:14)
new.ids <- c("macro","CD16mono","macro","macro","mast","mast","mast","mast",
             "mast","macro","unknown","mast","CD16mono","unknown","mast")
mye$ct <- plyr::mapvalues(mye$seurat_clusters,from = curr.ids,to = new.ids)

p.mye.ct <- DimPlot(mye,group.by = "ct") +
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/p.mye.ct.pdf",p.mye.ct,width = 5,height = 4)

saveRDS(mye,file = "validation/mye.annot.rds")

# macro
macro <- subset(mye,cells = rownames(subset(macro@meta.data, ct=="macro")))
macro <- NormalizeData(macro, normalization.method = "LogNormalize",scale.factor = 10000)
macro <- FindVariableFeatures(macro,selection.method = "vst")
macro <- ScaleData(object = macro, features = VariableFeatures(macro))
macro <- RunPCA(object = macro, features = VariableFeatures(macro))
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(macro)
# Find Neighbors + Find Clusters
macro <- FindNeighbors(macro, dims = 1:PC.num)
macro <- FindClusters(macro, resolution = 1.2)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
macro <- RunUMAP(macro, dims = 1:PC.num)
macro <- RunTSNE(macro,dims = 1:PC.num)

p.macro <- DimPlot(macro,reduction = "umap") + NoLegend() +
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/macro.pdf",p.macro,width = 5,height = 5)

rtm.fea <- c("RGS1","HLA-DPA1","CD74","C1QC","FCGBP","SEPP1","CSF1R","TREM2","S100A4")
p.rtm <- VlnPlot(macro,features = rtm.fea,stacked = T,
                 pt.size = 0,combine = T,direction = "horizontal",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave("validation/rtm.fea.pdf",p.rtm,width = 12,height = 6)

rtm <- subset(macro,cells = rownames(subset(macro@meta.data,seurat_clusters == 0 | seurat_clusters == 1)))
rtm <- NormalizeData(rtm, normalization.method = "LogNormalize",scale.factor = 10000)
rtm <- FindVariableFeatures(rtm,selection.method = "vst")
rtm <- ScaleData(object = rtm, features = VariableFeatures(rtm))
rtm <- RunPCA(object = rtm, features = VariableFeatures(rtm))
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(rtm)
# Find Neighbors + Find Clusters
rtm <- FindNeighbors(rtm, dims = 1:PC.num)
rtm <- FindClusters(rtm, resolution = 1.2)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
rtm <- RunUMAP(rtm, dims = 1:PC.num)
rtm <- RunTSNE(rtm,dims = 1:PC.num)
p.rtm <- DimPlot(rtm,reduction = "tsne")
ggsave("validation/rtm.pdf",p.rtm,width = 4.2,height = 4)

p1 <- FeaturePlot(rtm,features = "HSP90AA1",reduction = "tsne",cols = c("grey88","DarkCyan"),pt.size = 0.01)

mytheme <- theme(panel.border = element_rect(size = 1.5,colour = "black"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "bottom")

cell.high <- rownames(subset(rtm@meta.data,seurat_clusters != 7))
rtm$Celltypes <- ifelse(rtm$seurat_clusters == 7, "GLUL-SQSTM1-","GLUL+SQSTM1+")
p.rtm <- DimPlot(rtm,reduction = "tsne") + NoLegend() +
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/rtm.pdf",p.rtm,width = 5,height = 5)

highlight <- DimPlot(rtm,reduction = "tsne",cells.highlight = cell.high,cols.highlight = "#DF78C1") + mytheme

ggsave("validation/GLUL.SQSTM1.rtm.pdf",highlight,width = 5,height = 5.2)

saveRDS(rtm, file = "validation/rtm.rds")

# mast
mast <- subset(mye,cells = rownames(subset(mye@meta.data, ct=="mast")))
mast <- NormalizeData(mast, normalization.method = "LogNormalize",scale.factor = 10000)
mast <- FindVariableFeatures(mast,selection.method = "vst")
mast <- ScaleData(object = mast, features = VariableFeatures(mast))
mast <- RunPCA(object = mast, features = VariableFeatures(mast))
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(mast)
# Find Neighbors + Find Clusters
mast <- FindNeighbors(mast, dims = 1:PC.num)
mast <- FindClusters(mast, resolution = 0.4)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
mast <- RunUMAP(mast, dims = 1:PC.num)
mast <- RunTSNE(mast,dims = 1:PC.num)

p.mast <- DimPlot(mast,reduction = "tsne") + NoLegend() +
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/mast.pdf",p.mast,width = 5,height = 5)

cell.high <- rownames(subset(mast@meta.data,seurat_clusters == 0 | seurat_clusters == 1 |
                               seurat_clusters == 2 | seurat_clusters == 3 |
                               seurat_clusters == 5 | seurat_clusters == 7 |
                               seurat_clusters == 9 | seurat_clusters == 10))
curr.ids <- c(0:17)
new.ids <- c("HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+",
             "HSP90AA1-HSP90AB1-","HSP90AA1+HSP90AB1+","HSP90AA1-HSP90AB1-","HSP90AA1+HSP90AB1+",
             "HSP90AA1-HSP90AB1-","HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+","HSP90AA1-HSP90AB1-",
             "HSP90AA1-HSP90AB1-","HSP90AA1-HSP90AB1-","HSP90AA1-HSP90AB1-","HSP90AA1-HSP90AB1-",
             "HSP90AA1-HSP90AB1-","HSP90AA1-HSP90AB1-")
mast$Celltype <- plyr::mapvalues(mast$seurat_clusters,from = curr.ids,to = new.ids)

highlight <- DimPlot(mast,reduction = "tsne",cells.highlight = cell.high,cols.highlight = "#00B050") + mytheme
ggsave("validation/HSP90AA1.HSP90AB1.Mast.pdf",highlight,width = 5,height = 5.2)

saveRDS(mast,file = "validation/mast.rds")

# CD16mono
CD16mono <- subset(mye,cells = rownames(subset(mye@meta.data, ct=="CD16mono")))
CD16mono <- NormalizeData(CD16mono, normalization.method = "LogNormalize",scale.factor = 10000)
CD16mono <- FindVariableFeatures(CD16mono,selection.method = "vst")
CD16mono <- ScaleData(object = CD16mono, features = VariableFeatures(CD16mono))
CD16mono <- RunPCA(object = CD16mono, features = VariableFeatures(CD16mono))
PCDeterminators <- function(object){
  stdev <- object@reductions$pca@stdev
  var <- stdev^2
  endVar <- 0
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    exp.var <- numerator/total
    if(endVar == 0){
      if(exp.var > 0.9){
        endVar <- endVar + 1
        PC.num <- i
      }
    }
  }
  sum(var[1:PC.num])/sum(var)
  return(PC.num)
}
PC.num <- PCDeterminators(CD16mono)
# Find Neighbors + Find Clusters
CD16mono <- FindNeighbors(CD16mono, dims = 1:PC.num)
CD16mono <- FindClusters(CD16mono, resolution = 1.5)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
CD16mono <- RunUMAP(CD16mono, dims = 1:PC.num)
CD16mono <- RunTSNE(CD16mono,dims = 1:PC.num)

p.CD16mono <- DimPlot(CD16mono,reduction = "tsne") + NoLegend() +
  theme(panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
ggsave("validation/CD16mono.pdf",p.CD16mono,width = 5,height = 5)

CD16mono$Celltype <- ifelse(CD16mono$seurat_clusters == 1, "JAK3-TLR4-","JAK3+TLR4+")
cell.high <- rownames(subset(CD16mono@meta.data,seurat_clusters != 1))
highlight <- DimPlot(CD16mono,reduction = "tsne",cells.highlight = cell.high,cols.highlight = "#FF7E07") + mytheme
ggsave("validation/JAK3.TLR4.CD16mono.pdf",highlight,width = 5,height = 5.2)
saveRDS(CD16mono,file = "validation/CD16mono.rds")
