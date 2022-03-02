library(Seurat)
library(MySeuratWrappers)
library(ggplot2)

total.integrated <- readRDS("input/total.integrated.harmony.rds")

DimRedu <- function(object){
  object$CellTypes <- as.character(object$CellTypes)
  object <- NormalizeData(object = object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object,selection.method="vst",nfeatures=2000)
  DefaultAssay(object) <- "integrated"
  object <- ScaleData(object, features = rownames(object))
  object <- RunPCA(object = object, features = VariableFeatures(object))
  # Run PCA and Determine Dimensions for 90% Variance
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
  PC.num <- PCDeterminators(object)
  # Find neighbors and clusters with harmony batch correction
  object <- FindNeighbors(object, dims = 1:PC.num, reduction = "pca")
  object <- FindClusters(object, resolution = 0.2, reduction = "pca")
  # myeloid <- myeloid %>% RunHarmony("SampleName",plot_convergence = F)
  object <- RunUMAP(object, dims = 1:PC.num, reduction = "pca")
  object <- RunTSNE(object,dims = 1:PC.num)
  return(object)
}


## ----- 1. T cells ---------

Tsub <- subset(total.integrated,cells = rownames(subset(total.integrated@meta.data, CellTypes=='T cells')))
saveRDS(Tsub,file = "input/cellSub/Tsub.rds")
Tsub <- DimRedu(Tsub)

mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "right",
                 plot.title = element_blank(),
                 legend.text = element_text(size = 14))


p1 <- DimPlot(Tsub,reduction = "tsne") + mytheme
ggsave(paste0("figures/celltypes/Tcell/","clusters.pdf"),p1,width = 4.5,height = 4)

p2 <- VlnPlot(Tsub,features = c("CD8A","CD8B","CD4"),stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/Tcell/","features.pdf"),p2,width = 3,height = 4)

curr.ids <- c(0:13)
new.ids <- c("CD4T","CD8T","Unknown","Unknown","CD4T","CD4T","CD4T","CD8T","CD8T","CD8T","CD8T","CD8T","Unknown","CD4T")
Tsub$CellTypes <- plyr::mapvalues(Tsub$seurat_clusters, from = curr.ids, to = new.ids)

p3 <- DimPlot(Tsub,reduction = "tsne",group.by = "CellTypes") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/Tcell/","celltypes.pdf"),p3,width = 3,height = 3.2)

saveRDS(Tsub, file = "input/cellSub/Tsub.cell.annot.rds")

# Epithelial subset
Epi <- subset(total.integrated,cells = rownames(subset(total.integrated@meta.data, CellTypes=='Epithelial')))
rm(total.integrated);gc()
saveRDS(Epi,file = "input/cellSub/Epithelial.rds")


## ---- 2.Mast cells subset ----

myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
mye.mast <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, CellTypes == "Mast cells")))
mye.mast <- DimRedu(mye.mast)
p1 <- DimPlot(mye.mast,reduction = "tsne") + mytheme 
ggsave(paste0("figures/celltypes/mast/","dimRedu.pdf"),p1,width = 4.4,height = 4)
# HSP90AB1
mypal <- c("grey88","DarkCyan")
DefaultAssay(mye.mast) <- "RNA"
f.HSP90AB1 <- FeaturePlot(mye.mast,features = "HSP90AB1",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/mast/","HSP90AB1.pdf"),f.HSP90AB1,width = 4.5,height = 4)
# HSP90AA1
f.HSP90AA1 <- FeaturePlot(mye.mast,features = "HSP90AA1",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/mast/","HSP90AA1.pdf"),f.HSP90AA1,width = 4.5,height = 4)

p2 <- VlnPlot(mye.mast,features = c("HSP90AB1","HSP90AA1"),stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/mast/","features.pdf"),p2,width = 2,height = 4)

curr.ids <- c(0:5)
new.ids <- c("HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+","HSP90AA1-HSP90AB1-","HSP90AA1+HSP90AB1+","HSP90AA1+HSP90AB1+")

mye.mast$CellTypes <- plyr::mapvalues(mye.mast$seurat_clusters,from = curr.ids,to = new.ids)
saveRDS(mye.mast,file = "input/cellSub/mye.pdac.mast.sub.cell.annot.rds")
p3 <- DimPlot(mye.mast,group.by = "CellTypes",reduction = "tsne") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/mast/","HSP90.pdf"),p3,width = 4,height = 4.2)
all.markers <- FindMarkers(mye.mast,ident.1 = "HSP90AA1+HSP90AB1+",ident.2 = "HSP90AA1-HSP90AB1-",group.by = "CellTypes")

all.markers <- all.markers[!grepl("^RP[SL]",rownames(all.markers)),]
all.markers <- all.markers[!grepl("^MT-",rownames(all.markers)),]
all.markers <- all.markers[order(all.markers[,2],decreasing = T),]
top20 <- rbind(head(all.markers,20),tail(all.markers,20))

up.genes <- c("TPSAB1","RGS1","AREG","CPA3","JUN","LTC4S","CLU","KIT","FAU","HSPA1B",
              "CD63","HSPA1A","NACA","HSP90AB1","HPGDS","VWA5A","DDIT4","CTSG")
p4 <- VlnPlot(mye.mast,features = up.genes,stacked = T,direction = "horizontal",pt.size = 0,
        group.by = "CellTypes.2",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/mast/","up18.pdf"),p4,width = 12,height = 3)

down.genes <- c("CLC","RUNX1","FAM101B","SORL1","PIM1","CSF3R","ATP10D","MAF","MYO1F",
                "IL17RA","MXD1","MS4A3","IL3RA","TP53INP1","FAM65B","S100A9")
p5 <- VlnPlot(mye.mast,features = down.genes,stacked = T,direction = "horizontal",pt.size = 0,
              group.by = "CellTypes.2",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/mast/","down16.pdf"),p5,width = 12,height = 3)

## correlation with CD8A & epithelial cells

# this cohort

## correlation with CD8A
CD8T <- readRDS("input/CD8T.rds")
gene.1 <- c("CD8A","CD8B")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
# up genes
Idents(mye.mast) <- mye.mast$SampleName
avgExp.T2 <- AverageExpression(mye.mast,assays = "RNA",slot = "data",features = up.genes)
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)
dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
# down genes
avgExp.T3 <- AverageExpression(mye.mast,assays = "RNA",slot = "data",features = down.genes)
avgExp.T3 <- data.frame(t(avgExp.T3$RNA),stringsAsFactors = F)
avgExp.T3$Sample <- rownames(avgExp.T3)
avgExp <- dplyr::left_join(avgExp.T3,avgExp,by = "Sample")

rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

## correlation with epithelial cells
Epigene <- c("EPCAM","ACTA2","KRT7","KRT8","KRT18","KRT19","CDH1","PRSS1","CTRB2",
             "REG1A","CLU","MKI67","SPINK1","TFF1","MUC1")
Epi <- readRDS("input/cellSub/Epithelial.rds")
Idents(Epi) <- Epi$SampleName
avgExp.E <- AverageExpression(Epi,assays = "RNA",slot = "data",features = Epigene)
avgExp.E <- data.frame(t(avgExp.E$RNA),stringsAsFactors = F)
avgExp.E$Sample <- rownames(avgExp.E)

avgExp <- dplyr::left_join(avgExp,avgExp.E,by = "Sample")
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

# up
for(i in 17:34){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[35])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/thisCorhort/CD8T/up/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

# down
for(i in 1:16){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[35])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/thisCorhort/CD8T/down/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

## epi
# up
for(i in 17:34){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[37])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/thisCorhort/Epi/up/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}
# down
for(i in 1:16){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[37])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/thisCorhort/Epi/down/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}
save(avgExp,file = "input/mye.mast.corr.with.CD8T.Epi.rda")

# correlation analysis by TCGA
load("input/tcga/TCGA_36cancers.rnaseq.RData")
PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
# gene <- c("HSP90AA1","HSP90AB1")
# up
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% up.genes)])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
for(i in 2:19){
  p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/TCGA/CD8T/up/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}
# down
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% down.genes)])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
for(i in 2:17){
  p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/mast/corr/TCGA/CD8T/down/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

## visualize mast subset splitBy samples
mast <- readRDS("input/cellSub/mye.pdac.mast.sub.cell.annot.rds")
p.mast <- DimPlot(mast, reduction = "tsne", label = F, pt.size = 0.2, split.by = "SampleName",
                  group.by = "CellTypes",ncol = 5)
ggsave("figures/celltypes/mast/splitBySampleName.pdf",p.mast,width = 15,height = 15)


## ---- 3.CD14 mono ----

myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
mye.CD14mono <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, CellTypes == "CD14 Mono.")))
mye.CD14mono <- DimRedu(mye.CD14mono)
p1 <- DimPlot(mye.CD14mono,reduction = "tsne") + mytheme
ggsave(paste0("figures/celltypes/CD14mono/","clusters.pdf"),p1,width = 4.5,height = 4)
# BIRC3
mypal <- c("grey88","DarkCyan")
DefaultAssay(mye.CD14mono) <- "RNA"
f.BIRC3 <- FeaturePlot(mye.CD14mono,features = "BIRC3",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/CD14mono/","BIRC3.pdf"),f.BIRC3,width = 4.5,height = 4)
p.BIRC3 <- VlnPlot(mye.CD14mono,features = "BIRC3",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/CD14mono/","features.pdf"),p.BIRC3,width = 2,height = 5)

curr.ids <- c(0:7)
new.ids <- c("BIRC3+","BIRC3+","BIRC3-","BIRC3-","BIRC3+","BIRC3-","BIRC3+","BIRC3+")
mye.CD14mono$CellTypes <- plyr::mapvalues(mye.CD14mono$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.CD14mono,reduction = "tsne",group.by = "CellTypes") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/CD14mono/","celltypes.pdf"),p.ct,width = 3,height = 3.2)
saveRDS(mye.CD14mono,file = "input/cellSub/mye.pdac.CD14mono.sub.cell.annot.rds")

## BIRC3 correlation with CD8A
CD8T <- readRDS("input/CD8T.rds")
gene.1 <- c("CD8A","CD8B")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
Idents(mye.CD14mono) <- mye.CD14mono$SampleName
avgExp.T2 <- AverageExpression(mye.CD14mono,assays = "RNA",slot = "data",features = "BIRC3")
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)
dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL
for(i in 2:3){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD14mono/",colnames(avgExp)[i],"_vs_","BIRC3.pdf"),p,width = 5,height = 5)
}


## ---- 4.CD16 mono ----
myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
mye.CD16mono <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, CellTypes == "CD16 Mono.")))
mye.CD16mono <- DimRedu(mye.CD16mono)
p1 <- DimPlot(mye.CD16mono,reduction = "tsne") + mytheme
ggsave(paste0("figures/celltypes/CD16mono/","clusters.pdf"),p1,width = 4.5,height = 4)
mypal <- c("grey88","DarkCyan")
DefaultAssay(mye.CD16mono) <- "RNA"
# BIRC3
f.BIRC3 <- FeaturePlot(mye.CD16mono,features = "BIRC3",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/CD16mono/","BIRC3.pdf"),f.BIRC3,width = 4.5,height = 4)
p.BIRC3 <- VlnPlot(mye.CD16mono,features = "BIRC3",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/CD16mono/","fea.BIRC3.pdf"),p.BIRC3,width = 2,height = 5)

curr.ids <- c(0:6)
new.ids <- c("BIRC3+","BIRC3-","BIRC3+","BIRC3+","BIRC3-","BIRC3+","BIRC3+")
mye.CD16mono$CellTypes.1 <- plyr::mapvalues(mye.CD16mono$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.CD16mono,reduction = "tsne",group.by = "CellTypes.1") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/CD16mono/","celltypes.pdf"),p.ct,width = 3,height = 3.2)
# JAK3
f.JAK3 <- FeaturePlot(mye.CD16mono,features = "JAK3",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/CD16mono/","JAK3.pdf"),f.JAK3,width = 4.5,height = 4)
p.JAK3 <- VlnPlot(mye.CD16mono,features = "JAK3",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/CD16mono/","fea.JAK3.pdf"),p.JAK3,width = 2,height = 5)

curr.ids <- c(0:6)
new.ids <- c("JAK3+","JAK3-","JAK3+","JAK3+","JAK3-","JAK3-","JAK3+")
mye.CD16mono$CellTypes.2 <- plyr::mapvalues(mye.CD16mono$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.CD16mono,reduction = "tsne",group.by = "CellTypes.2") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/CD16mono/","celltypes.2.pdf"),p.ct,width = 3,height = 3.2)

# PPIA
Idents(mye.CD16mono) <- mye.CD16mono$seurat_clusters
f.PPIA <- FeaturePlot(mye.CD16mono,features = "PPIA",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/CD16mono/","PPIA.pdf"),f.PPIA,width = 4.5,height = 4)
p.PPIA <- VlnPlot(mye.CD16mono,features = "PPIA",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/CD16mono/","fea.PPIA.pdf"),p.PPIA,width = 2,height = 5)

# TLR4
f.TLR4 <- FeaturePlot(mye.CD16mono,features = "TLR4",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/CD16mono/","TLR4.pdf"),f.TLR4,width = 4.5,height = 4)
p.TLR4 <- VlnPlot(mye.CD16mono,features = "TLR4",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/CD16mono/","fea.TLR4.pdf"),p.TLR4,width = 2,height = 5)

curr.ids <- c(0:6)
new.ids <- c("TLR4+","TLR4-","TLR4+","TLR4+","TLR4-","TLR4-","TLR4+")
mye.CD16mono$CellTypes.3 <- plyr::mapvalues(mye.CD16mono$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.CD16mono,reduction = "tsne",group.by = "CellTypes.3") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/CD16mono/","celltypes.3.pdf"),p.ct,width = 3,height = 3.2)

curr.ids <- c(0:6)
new.ids <- c("JAK3+TLR4+","JAK3-TLR4-","JAK3+TLR4+","JAK3+TLR4+","JAK3-TLR4-","JAK3-TLR4-","JAK3+TLR4+")
mye.CD16mono$CellTypes.4 <- plyr::mapvalues(mye.CD16mono$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.CD16mono,reduction = "tsne",group.by = "CellTypes.4") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/CD16mono/","celltypes.4.pdf"),p.ct,width = 3,height = 3.2)

saveRDS(mye.CD16mono,file = "input/cellSub/mye.pdac.CD16mono.sub.cell.annot.rds")

## BIRC3 correlation with CD8A
# this cohort
CD8T <- readRDS("input/CD8T.rds")
gene.1 <- c("CD8A","CD8B")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
Idents(mye.CD16mono) <- mye.CD16mono$SampleName
avgExp.T2 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = "BIRC3")
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)
dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[1], y = colnames(avgExp)[2])) +
  geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
  ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
ggsave(paste0("figures/celltypes/CD16mono/",colnames(avgExp)[2],"_vs_","BIRC3.pdf"),p,width = 5,height = 5)

# TCGA
load("input/tcga/TCGA_36cancers.rnaseq.RData")
PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% "BIRC3")])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
p <- ggplot(dt,aes_string(x = colnames(dt)[2], y = colnames(dt)[1])) +
  geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
  ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
ggsave(paste0("figures/celltypes/CD16mono/",colnames(dt)[2],"_vs_","CD8A.tcga.pdf"),p,width = 5,height = 5)

all.markers <- FindMarkers(mye.CD16mono,ident.1 = "BIRC3+",ident.2 = "BIRC3-",group.by = "CellTypes.1")
all.markers <- all.markers[!grepl("^RP[SL]",rownames(all.markers)),]
all.markers <- all.markers[!grepl("^MT-",rownames(all.markers)),]
all.markers <- all.markers[order(all.markers[,2],decreasing = T),]
top20 <- rbind(head(all.markers,20),tail(all.markers,20))

up.genes <- c(head(rownames(all.markers),20))
p.up <- VlnPlot(mye.CD16mono,features = up.genes,stacked = T,group.by = "CellTypes.1",y.max = 5,pt.size = 0,
        direction = "horizontal",x.lab = "",y.lab = "") +  NoLegend() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/CD16mono/","up20.pdf"),p.up,width = 12,height = 3)
dn.genes <- tail(rownames(all.markers),20)
p.dn <- VlnPlot(mye.CD16mono,features = dn.genes,stacked = T,group.by = "CellTypes.1",y.max = 5,
        direction = "horizontal",x.lab = "",y.lab = "",pt.size = 0.001) +  NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("BIRC3+","BIRC3-")), label = "p.signif") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/CD16mono/","down20.pdf"),p.dn,width = 20,height = 3)


# correlation with CD8A & epithelial cells
# this cohort
## correlation with CD8A
# CD8T <- readRDS("input/CD8T.rds")
# gene.1 <- c("CD8A","CD8B")
# avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
# avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
# avgExp.T$Sample <- rownames(avgExp.T)
# # up genes
# Idents(mye.CD16mono) <- mye.CD16mono$SampleName
# avgExp.T2 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = up.genes)
# avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
# avgExp.T2$Sample <- rownames(avgExp.T2)
# dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
# # down genes
# avgExp.T3 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = dn.genes)
# avgExp.T3 <- data.frame(t(avgExp.T3$RNA),stringsAsFactors = F)
# avgExp.T3$Sample <- rownames(avgExp.T3)
# avgExp <- dplyr::left_join(avgExp.T3,avgExp,by = "Sample")
# 
# ## correlation with epithelial cells
# Epigene <- c("EPCAM")
# Epi <- readRDS("input/cellSub/Epithelial.rds")
# Idents(Epi) <- Epi$SampleName
# avgExp.E <- AverageExpression(Epi,assays = "RNA",slot = "data",features = Epigene)
# avgExp.E <- data.frame(t(avgExp.E$RNA),stringsAsFactors = F)
# avgExp.E$Sample <- rownames(avgExp.E)
# 
# avgExp <- dplyr::left_join(avgExp,avgExp.E,by = "Sample")
# rownames(avgExp) <- avgExp$Sample
# avgExp$Sample <- NULL
# 
# if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/"))) dir.create("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/",recursive = T)
# # up
# for(i in 21:40){
#   p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[41])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
# }
# 
# if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/"))) dir.create("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/",recursive = T)
# # down
# for(i in 1:20){
#   p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[41])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
# }
# 
# if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/")))
#   dir.create("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/",recursive = T)
# ## epi
# # up
# for(i in 21:40){
#   p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[43])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
# }
# 
# if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/")))
#   dir.create("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/",recursive = T)
# down
# for(i in 1:20){
#   p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[43])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
# }
# save(avgExp,file = "input/mye.CD16mono.corr.with.CD8T.Epi.rda")


# correlation analysis by TCGA
# load("input/tcga/TCGA_36cancers.rnaseq.RData")
# PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
# # gene <- c("HSP90AA1","HSP90AB1")
# # up
# in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% up.genes)])
# dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
# rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
# for(i in 2:19){
#   p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/CD16mono/corr/TCGA/CD8T/up/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
# }
# # down
# in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% dn.genes)])
# dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
# rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
# for(i in 2:17){
#   p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
#     geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
#     theme_bw() +
#     theme(axis.title = element_text(size = 16),
#           axis.text = element_text(size = 14),
#           axis.ticks.length = unit(0.25, 'cm'),
#           axis.ticks = element_line(size = 1),
#           panel.border = element_rect(size = 1.5),
#           panel.grid = element_blank()) +
#     geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
#     ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
#   ggsave(paste0("figures/celltypes/mast/corr/TCGA/CD8T/down/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
# }

## JAK3&TLR4 correlation with CD8A
# this cohort
CD8T <- readRDS("input/CD8T.rds")
gene.1 <- c("CD8A")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
Idents(mye.CD16mono) <- mye.CD16mono$SampleName
avgExp.T2 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = c("JAK3","TLR4"))
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)
dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL
for(i in 1:2){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[3])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

# in TCGA
load("input/tcga/TCGA_36cancers.rnaseq.RData")
PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% c("JAK3","TLR4"))])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
for(i in 2:3){
  p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/",colnames(dt)[i],"_vs_","CD8A.tcga.pdf"),p,width = 5,height = 5)
}

## epithelial cells
Epigene <- c("EPCAM")
Epi <- readRDS("input/cellSub/Epithelial.rds")
Idents(Epi) <- Epi$SampleName
avgExp.E <- AverageExpression(Epi,assays = "RNA",slot = "data",features = Epigene)
avgExp.E <- data.frame(t(avgExp.E$RNA),stringsAsFactors = F)
avgExp.E$Sample <- rownames(avgExp.E)

Idents(mye.CD16mono) <- mye.CD16mono$SampleName
avgExp.T2 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = c("JAK3","TLR4"))
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)
dplyr::left_join(avgExp.T2,avgExp.E,by = "Sample") -> avgExp
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL
for(i in 1:2){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[3])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}

## Find Differentially expressed genes
all.markers <- FindMarkers(mye.CD16mono,ident.1 = "JAK3+TLR4+",ident.2 = "JAK3-TLR4-",group.by = "CellTypes.4")
all.markers <- all.markers[!grepl("^RP[SL]",rownames(all.markers)),]
all.markers <- all.markers[!grepl("^MT-",rownames(all.markers)),]
all.markers <- all.markers[order(all.markers[,2],decreasing = T),]
top20 <- rbind(head(all.markers,20),tail(all.markers,20))

up.genes <- c("CRIP1","GABARAP","IFITM1","ZBTB7A","ITGB1","KRTCAP2","ZYX","MIF","WARS","MX1","IFI6")
p.up <- VlnPlot(mye.CD16mono,features = up.genes,stacked = T,group.by = "CellTypes.4",y.max = 5,pt.size = 0,
                direction = "horizontal",x.lab = "",y.lab = "") +  NoLegend() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  ggpubr::stat_compare_means(comparisons = list(c("JAK3+TLR4+","JAK3-TLR4-")), label = "p.signif")
ggsave(paste0("figures/celltypes/CD16mono/","up20.2.pdf"),p.up,width = 12,height = 3)
dn.genes <- c("TNF","TMEM66","NPM1","RGS2","MS4A4A","OCIAD1","BTF3","EEF1D")
p.dn <- VlnPlot(mye.CD16mono,features = dn.genes,stacked = T,group.by = "CellTypes.4",y.max = 6,
                direction = "horizontal",x.lab = "",y.lab = "",pt.size = 0.001) +  NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("JAK3+TLR4+","JAK3-TLR4-")), label = "p.signif") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/CD16mono/","down20.2.pdf"),p.dn,width = 16,height = 3)

# correlation with CD8A & epithelial cells
# this cohort
## correlation with CD8A
CD8T <- readRDS("input/cellSub/CD8T.rds")
gene.1 <- c("CD8A")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
# up genes
Idents(mye.CD16mono) <- mye.CD16mono$SampleName

DataFrame <- function(object,feature){
  tmp <- AverageExpression(object,assays = "RNA",slot = "data",features = feature)
  tmp <- data.frame(t(tmp$RNA),stringsAsFactors = F)
  tmp$Sample <- rownames(tmp)
  return(tmp)
}
avgExp.T2 <- DataFrame(mye.CD16mono,feature=up.genes)

dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
# down genes
avgExp.T3 <- AverageExpression(mye.CD16mono,assays = "RNA",slot = "data",features = dn.genes)
avgExp.T3 <- data.frame(t(avgExp.T3$RNA),stringsAsFactors = F)
avgExp.T3$Sample <- rownames(avgExp.T3)
avgExp <- dplyr::left_join(avgExp.T3,avgExp,by = "Sample")

## correlation with epithelial cells
Epigene <- c("EPCAM")
Epi <- readRDS("input/cellSub/Epithelial.rds")
Idents(Epi) <- Epi$SampleName
avgExp.E <- AverageExpression(Epi,assays = "RNA",slot = "data",features = Epigene)
avgExp.E <- data.frame(t(avgExp.E$RNA),stringsAsFactors = F)
avgExp.E$Sample <- rownames(avgExp.E)

avgExp <- dplyr::left_join(avgExp,avgExp.E,by = "Sample")
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/"))) dir.create("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/",recursive = T)
# up
for(i in 9:23){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[24])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/up/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/"))) dir.create("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/",recursive = T)
# down
for(i in 1:8){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[24])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/CD8T/down/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

## epi
if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/")))
  dir.create("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/",recursive = T)
# up
for(i in 9:23){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[25])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/up/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}

if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/")))
  dir.create("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/",recursive = T)
# down
for(i in 1:8){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[25])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/thisCorhort/Epi/down/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}

# correlation analysis by TCGA
if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/TCGA/CD8T/up/")))
  dir.create("figures/celltypes/CD16mono/corr/TCGA/CD8T/up/",recursive = T)
load("input/tcga/TCGA_36cancers.rnaseq.RData")
PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
# gene <- c("HSP90AA1","HSP90AB1")
# up
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% up.genes)])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
for(i in 2:12){
  p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/TCGA/CD8T/up/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}
# down
if(!dir.exists(paste0("figures/celltypes/CD16mono/corr/TCGA/CD8T/down/")))
  dir.create("figures/celltypes/CD16mono/corr/TCGA/CD8T/down/",recursive = T)
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% dn.genes)])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
for(i in 2:8){
  p <- ggplot(dt,aes_string(x = colnames(dt)[i], y = colnames(dt)[1])) +
    geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25, 'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
    ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
  ggsave(paste0("figures/celltypes/CD16mono/corr/TCGA/CD8T/down/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

## visualize CD16mono subset splitBy samples
CD16mono <- readRDS("input/cellSub/mye.pdac.CD16mono.sub.cell.annot.rds")
p.CD16mono <- DimPlot(CD16mono, reduction = "tsne", label = F, pt.size = 0.2, split.by = "SampleName",
                      group.by = "CellTypes.4",ncol = 5)
ggsave("figures/celltypes/CD16mono/splitBySampleName.pdf",p.CD16mono,width = 15,height = 15)


## ---- 5.DCs ----
myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
mye.DC <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, CellTypes == "DCs")))
mye.DC <- DimRedu(mye.DC)
p1 <- DimPlot(mye.DC,reduction = "tsne") + mytheme
ggsave(paste0("figures/celltypes/DCs/","clusters.pdf"),p1,width = 4.5,height = 4)
mypal <- c("grey88","DarkCyan")
DefaultAssay(mye.DC) <- "RNA"
# CHMP1B
f.CHMP1B <- FeaturePlot(mye.DC,features = "CHMP1B",reduction = "tsne",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/DCs/","CHMP1B.pdf"),f.CHMP1B,width = 4.5,height = 4)
p.CHMP1B <- VlnPlot(mye.DC,features = "CHMP1B",stacked = T,direction = "horizontal",pt.size = 0,x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.border = element_rect(size = 1,colour = "black")) +
  NoLegend()
ggsave(paste0("figures/celltypes/DCs/","fea.CHMP1B.pdf"),p.CHMP1B,width = 2,height = 5)

curr.ids <- c(0:7)
new.ids <- c("CHMP1B-","CHMP1B-","CHMP1B-","CHMP1B-","CHMP1B-","CHMP1B-","CHMP1B+","CHMP1B-")
mye.DC$CellTypes.1 <- plyr::mapvalues(mye.DC$seurat_clusters,from = curr.ids,to = new.ids)
p.ct <- DimPlot(mye.DC,reduction = "tsne",group.by = "CellTypes.1") + mytheme + theme(legend.position = "bottom")
ggsave(paste0("figures/celltypes/DCs/","celltypes.pdf"),p.ct,width = 3,height = 3.2)

saveRDS(mye.DC,file = "input/cellSub/mye.pdac.DC.sub.cell.annot.rds")

# correlation with CD8A & epithelial cells
# this cohort
## correlation with CD8A
CD8T <- readRDS("input/cellSub/CD8T.rds")
gene.1 <- c("CD8A")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
# up genes
Idents(mye.DC) <- mye.DC$SampleName
DataFrame <- function(object,feature){
  tmp <- AverageExpression(object,assays = "RNA",slot = "data",features = feature)
  tmp <- data.frame(t(tmp$RNA),stringsAsFactors = F)
  tmp$Sample <- rownames(tmp)
  return(tmp)
}
avgExp.T2 <- DataFrame(mye.DC,feature="CHMP1B")

dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp

rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[1], y = colnames(avgExp)[2])) +
  geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
  ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
ggsave(paste0("figures/celltypes/DCs/",colnames(avgExp)[1],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
