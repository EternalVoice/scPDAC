## ------------------ Myeloid subset ---------------------------

## ---- 1.Dimensional Reduction & Cell Type Mapping ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(rlang)

total.integrated.harmony <- readRDS("input/total.integrated.harmony.rds")
myeloid.cells <- rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Myeloid" | CellTypes == "Mast cells"))
myeloid <- subset(total.integrated.harmony, cells = myeloid.cells)
myeloid <- NormalizeData(object = myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
myeloid <- FindVariableFeatures(myeloid,selection.method="vst",nfeatures=2000)
DefaultAssay(myeloid) <- "integrated"
myeloid <- ScaleData(myeloid, features = rownames(myeloid))
myeloid <- RunPCA(object = myeloid, features = VariableFeatures(myeloid))
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
PC.num <- PCDeterminators(myeloid)
# Find neighbors and clusters with harmony batch correction
myeloid <- FindNeighbors(myeloid, dims = 1:PC.num, reduction = "pca")
myeloid <- FindClusters(myeloid, resolution = 0.2, reduction = "pca")
# myeloid <- myeloid %>% RunHarmony("SampleName",plot_convergence = F)
myeloid <- RunUMAP(myeloid, dims = 1:PC.num, reduction = "pca")
myeloid <- RunTSNE(myeloid,dims = 1:PC.num)

saveRDS(myeloid, file = "input/mye.with.mast.dr.rds")

mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "bottom",
                 plot.title = element_blank(),
                 legend.text = element_text(size = 14))
p.myeloid <- DimPlot(myeloid,reduction = "tsne",label = T) + mytheme
ggsave(paste0("figures/celltypes/myeloid/","myeloid.pdf"),p.myeloid,width = 5,height = 5.8)
# Mast cells: KIT, CPA3, TPSAB1
# Granulocyte: CXCR2, FCGR3B, IFITM2, SLC25A37, IL1R2, CXCR1, SIRPA, S100A8
# Mono_CD14: FCN1, S100A8, S100A9
# Mono_CD16: FCGR3A, LST1, LILRB2
# Macrophage:
#     Macro_C1QC: C1QC, C1QA, APOE
#     Macro_INHBA: INHBA, IL1RN, CCL4
#     Macro_NLRP3: NLRP3, EREG, IL1B
#     Macro_LYVE1: LYVE1, PLTP, SEPP1
# 
# DCs:
#     pDC_LILRA4: LILRA4, GZMB, IL3RA, CLIC3, IRF7
#     cDC1_CLEC9A: CLEC9A, FLT3, IDO1, CADM1, CLNK
#     cDC2_CD1C: CD1C, FCER1A, HLA-DQA1,LYZ, CD1A, IL22RA2, LTB, CD52, CASP1
#     cDC3_LAMP3: LAMP3, CCR7, FSCN1, CCL22, IL7R, MMP7, IFIT2

DefaultAssay(myeloid) <- "RNA"
p.mast <- DotPlot(myeloid,features = c("KIT","CPA3","TPSAB1")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","mast.pdf"),p.mast,width = 4,height = 5)
p.gra <- DotPlot(myeloid,features = c("CXCR2","FCGR3B","IFITM2","SLC25A37","IL1R2","CXCR1","SIRPA","S100A8")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","granulocyte.pdf"),p.gra,width = 8,height = 4)
# Monocyte
p.Mono.CD14 <- DotPlot(myeloid, features = c("CD14","FCN1","S100A8","S100A9")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","Mono.CD14.pdf"),p.Mono.CD14,width = 5,height = 5)
p.Mono.CD16 <- DotPlot(myeloid,features = c("FCGR3A","LST1","LILRB2")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","Mono.CD16.pdf"),p.Mono.CD16,width = 4,height = 5)
# Macrophage
p.macro.C1QC <- DotPlot(myeloid, features = c("C1QC","C1QA","APOE")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","macro.C1QC.pdf"),p.macro.C1QC,width = 4,height = 5)
p.macro.INHBA <- DotPlot(myeloid, features = c("INHBA","IL1RN","CCL4")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","macro.INHBA.pdf"),p.macro.INHBA,width = 4,height = 5)
p.macro.NLRP3 <- DotPlot(myeloid, features = c("NLRP3","EREG","IL1B")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","macro.NLRP3.pdf"),p.macro.NLRP3,width = 4,height = 5)
p.macro.LYVE1 <- DotPlot(myeloid, features = c("LYVE1","PLTP","SEPP1")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","macro.LYVE1.pdf"),p.macro.LYVE1,width = 4,height = 5)
# Dentritic cells
p.pDC.LILRA4 <- DotPlot(myeloid, features = c("LILRA4","GZMB","IL3RA","CLIC3","IRF7")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","pDC.LILRA4.pdf"),p.pDC.LILRA4,width = 6,height = 5)
p.cDC1.CLEC9A <- DotPlot(myeloid, features = c("CLEC9A","FLT3","IDO1","CADM1","CLNK")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","cDC1.CLEC9A.pdf"),p.cDC1.CLEC9A,width = 6,height = 5)
p.cDC2.CD1C <- DotPlot(myeloid, features = c("CD1C","FCER1A","HLA-DQA1","LYZ","CD1A","IL22RA2","LTB","CD52","CASP1")) + theme_bw()
ggsave(paste0("figures/celltypes/myeloid/","cDC2.CD1C.pdf"),p.cDC2.CD1C,width = 9,height = 5)
p.cDC3.LAMP3 <- DotPlot(myeloid, features = c("LAMP3","CCR7","FSCN1","CCL22","IL7R","MMP7","IFIT2"))
ggsave(paste0("figures/celltypes/myeloid/","cDC3.LAMP3.pdf"),p.cDC3.LAMP3,width = 8,height = 5)

## Find markers
Idents(myeloid) <- myeloid$seurat_clusters
mye.all.markers <- FindAllMarkers(myeloid,only.pos = T,min.pct = 0.25)
save(mye.all.markers, file = "output/mye.all.markers.rda")

# C10: IL7R, IL32, CCR7, GZMB
# C11: LYZ, IL3RA
curr.ids <- c(0:13)
new.ids <- c("Granulocyte","CD14 Mono.","Granulocyte","Macrophage","CD14 Mono.","Granulocyte","Macrophage",
             "Mast cells", "CD16 Mono.", "Macrophage", "DCs", "DCs", "Granulocyte", "Granulocyte")
myeloid$CellTypes <- plyr::mapvalues(myeloid$seurat_clusters, from = curr.ids, to = new.ids)

saveRDS(myeloid, file = "input/mye.with.mast.annot.rds")

# Macrophage  #DF78C1
# Mast cells  #00B050
# CD14 Mono. #DB4840
# CD16 Mono. #FF7E07
# DCs        #7030A0
mycols <- c("#4E9BD3","#DB4840","#DF78C1","#00B050","#FF7E07","#7030A0")
# mycols <- c('#1f77b4', '#39d486', '#43315c', '#00c8ff', '#de4343', '#ff7f0e', '#e377c2')

# total (PDAC + Normal)
p.celltype <- DimPlot(myeloid, reduction = "tsne", group.by = "CellTypes", cols = mycols) + 
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(paste0("figures/celltypes/myeloid/","celltypes.pdf"),p.celltype,width = 5.7,height = 6)

## ------------ 2.Cell component comparison ----------------------
# patient vs normal
p.patient_vs_normal <- DimPlot(myeloid, reduction = "tsne", split.by = "Patient",cols = mycols) + NoLegend()
ggsave(paste0("figures/celltypes/myeloid/","p.patient_vs_normal.pdf"),p.patient_vs_normal,width = 8,height = 4.6)
# class
p.class <- DimPlot(myeloid, reduction = "tsne", split.by = "Class",cols = mycols,pt.size = 0.8) + NoLegend()
ggsave(paste0("figures/celltypes/myeloid/","p.class.pdf"),p.class,width = 16,height = 4.6)

## PDAC subset: primary, metastatic, paracancerous
myeloid.pdac <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, Patient == "PDAC")))
myeloid.pdac$State <- ifelse(myeloid.pdac$Class=="Paracancerous", "Para-tumor", "tumor")
# visualize celltypes
p.celltype.pdac <- DimPlot(myeloid.pdac, reduction = "tsne", group.by = "CellTypes", cols = mycols) + 
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(paste0("figures/celltypes/myeloid/","celltypes.pdac.pdf"),p.celltype.pdac,width = 5.7,height = 6)
saveRDS(myeloid.pdac, file = "input/mye.with.mast.pdac.sub.rds")
# PBMC: PDAC & Normal
myeloid.pbmc <- subset(myeloid, cells = rownames(subset(myeloid@meta.data, Tissue == "PBMC")))
# visualize celltypes
p.celltype.pbmc <- DimPlot(myeloid.pbmc, reduction = "tsne", group.by = "CellTypes", cols = mycols) + 
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        plot.title = element_blank(),
        legend.text = element_text(size = 14))
ggsave(paste0("figures/celltypes/myeloid/","celltypes.pbmc.pdf"),p.celltype.pbmc,width = 5.7,height = 6)
saveRDS(myeloid.pbmc, file = "input/mye.with.mast.pbmc.sub.rds")

## ----------- Myeloid: cell heterogeneity ---------------

mytheme.1 <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                   panel.grid = element_blank(),
                   axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
                   axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
                   legend.position = "top", legend.direction = "horizontal",
                   legend.text = element_text(size = 14),legend.title = element_text(size = 15))
# PDAC: (primary, metatastic, paracanerous)
p.class <- ggplot(myeloid.pdac@meta.data, aes(x = Class, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(paste0("figures/celltypes/myeloid/","mye.pdac.cell.prop_vs_class.pdf"),p.class, width = 11,height = 2)
# PDAC: (Tissue, PBMC)
p.ts <- ggplot(myeloid.pdac@meta.data, aes(x = Tissue, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(paste0("figures/celltypes/myeloid/","mye.pdac.cell.prop_vs_tissue.pdf"),p.ts, width = 11,height = 2)
# PBMC: (PDAC, Healthy)
p.pbmc.tumor_vs_normal <- ggplot(myeloid.pbmc@meta.data, aes(x = Patient, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(paste0("figures/celltypes/myeloid/","mye.pbmc.cell.prop_vs_patient.pdf"),p.pbmc.tumor_vs_normal, width = 11,height = 2)

# Patient vs. Celltypes
# p.pt <- ggplot(myeloid@meta.data, aes(x = Patient, fill = CellTypes)) +
#   geom_bar(position = "fill") + theme_bw() + coord_flip() +
#   scale_fill_manual(values = mycols) +
#   mytheme.1 +
#   guides(fill = guide_legend(nrow = 1)) +
#   labs(x = "Patient", y = "Cell Fraction")
# ggsave(paste0("figures/celltypes/myeloid/","cell.prop_vs_Patient.pdf"),p.pt, width = 11,height = 2)

# Tissue vs. Celltypes
# p.ts <- ggplot(myeloid@meta.data, aes(x = Tissue, fill = CellTypes)) +
#   geom_bar(position = "fill") + theme_bw() + coord_flip() +
#   scale_fill_manual(values = mycols) +
#   mytheme.1 +
#   guides(fill = guide_legend(nrow = 1)) +
#   labs(x = "Tissue", y = "Cell Fraction")
# ggsave(paste0("figures/celltypes/myeloid/","cell.prop_vs_Tissue.pdf"),p.ts, width = 11,height = 2)

## ----------- immune-suppressive markers check ---------------
rm(list = ls());gc()
myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")

mytheme.2 <- theme(panel.border = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   legend.position = "none")
mycols.2 <- c("grey90","#39d486")
# SPP1
p.SPP1 <- FeaturePlot(myeloid, features = "SPP1", cols = mycols.2, reduction = "tsne") + mytheme.2
ggsave(paste0("figures/celltypes/myeloid/featurePlot/","SPP1.pdf"),p.SPP1,width = 4,height = 4.2)
# MARCO
p.MARCO <- FeaturePlot(myeloid, features = "MARCO", cols = mycols.2, reduction = "tsne") + mytheme.2
ggsave(paste0("figures/celltypes/myeloid/featurePlot/","MARCO.pdf"),p.MARCO,width = 4,height = 4.2)
# APOE
p.APOE <- FeaturePlot(myeloid, features = "APOE", cols = mycols.2, reduction = "tsne") + mytheme.2
ggsave(paste0("figures/celltypes/myeloid/featurePlot/","APOE.pdf"),p.APOE,width = 4,height = 4.2)
# CD68
p.CD68 <- FeaturePlot(myeloid, features = "CD68", cols = mycols.2, reduction = "tsne") + mytheme.2
ggsave(paste0("figures/celltypes/myeloid/featurePlot/","CD68.pdf"),p.CD68,width = 4,height = 4.2)
# SIRPA
p.SIRPA <- FeaturePlot(myeloid, features = "SIRPA", cols = mycols.2, reduction = "tsne") + mytheme.2
ggsave(paste0("figures/celltypes/myeloid/featurePlot/","SIRPA.pdf"),p.SIRPA,width = 4,height = 4.2)


## ----------- Visualize marker genes ---------------
rm(list = ls());gc()
myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
# Mast: KIT, CPA3, TPSAB1
# Granulocyte: FCGR3B, IFITM2, CXCR2, S100A8, SLC25A37
# CD14 Mono.: CD14, FCN1, S100A8, S100A9
# CD16 Mono.: FCGR3A, LST1, LILRB2
# Macrophage: CCL4, C1QC, C1QA, APOE
# DCs: IL7R, IL32, CCR7, GZMB, IL8
features <- c("KIT","CPA3","TPSAB1",
              "FCGR3B","IFITM2","CXCR2","SLC25A37",
              "CD14","FCN1","S100A8","S100A9",
              "FCGR3A","LST1","LILRB2",
              "CCL4","C1QC","C1QA","APOE",
              "IL7R","IL32","CCR7","GZMB","IL8")

p.features <- DotPlot(myeloid, features = features,group.by = "CellTypes") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 14,face = "bold"),
        axis.text.y = element_text(size = 14,face = "bold"),
        panel.border = element_rect(size = 1.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlBu"))) +
  guides(size = guide_legend(title = "pct.exp"), color = guide_colorbar(title = "ave.exp"))
ggsave(paste0("figures/celltypes/myeloid/","features.pdf"),p.features,width = 10,height = 4)

# p.features.2 <- p.features + theme(legend.position = "none")
# ggsave(paste0("figures/celltypes/myeloid/","features.2.pdf"),p.features.2,width = 12,height = 3)

## ----------- 3.Survival validation ---------------
rm(list = ls());gc()
#' @param TCGA.exp.mtx TCGA expression matrix
#' @param clinical.data Paired clinical data
#' @param signatures Signatures to validate
#' @param title legend title for plot
#' 
dataProcess <- function(TCGA.exp.mtx, clinical.data, signatures){
  expr <- read.table(TCGA.exp.mtx, sep = "\t", header = T)
  expr <- expr[!duplicated(expr$Hugo_Symbol),]
  rownames(expr) <- expr$Hugo_Symbol
  expr <- expr[, -c(1,2)]
  colnames(expr) <- gsub("\\.", "-", colnames(expr))
  expr <- data.frame(t(expr))
  # just grab genes in signatures
  expr.2 <- data.frame(expr[,signatures])
  rownames(expr.2) <- rownames(expr)
  colnames(expr.2) <- gsub("expr.","",colnames(expr.2))
  expr.3 <- data.frame(scale(expr.2))
  expr.3$PATIENT_ID <- gsub("-01","",rownames(expr.3))
  # average over those with multiple measurements
  ave.exp <- aggregate(.~PATIENT_ID, FUN = mean, data = expr.3)
  
  # Prepare clinical data
  clin <- read.table(clinical.data, skip = 4, header = T, sep = "\t")
  small_clin <- data.frame(clin$PATIENT_ID, clin$OS_MONTHS, clin$OS_STATUS)
  small_clin <- small_clin[small_clin$clin.OS_STATUS != "[Not Available]",]
  # create binary variable if dead or alive
  Dead <- rep(NA, dim(small_clin)[1])
  for(i in 1:length(small_clin$clin.OS_STATUS)){
    if(stringr::str_split(small_clin$clin.OS_STATUS,":",simplify = T)[,2][i] == "DECEASED"){
      Dead[i] <- 1
    }else if(stringr::str_split(small_clin$clin.OS_STATUS,":",simplify = T)[,2][i] == "LIVING"){
      Dead[i] <- 0
    }
  }
  small_clin$Dead <- Dead
  small_clin$clin.OS_MONTHS <- as.integer(small_clin$clin.OS_MONTHS)
  # merge data
  merged.data <- merge(small_clin, ave.exp, by.x = "clin.PATIENT_ID", by.y = "PATIENT_ID")
  
  # survival analysis
  coxfit <- survival::coxph(survival::Surv(clin.OS_MONTHS, Dead) ~., data = merged.data[,c(2,4:ncol(merged.data))])
  coefs <- as.matrix(coxfit$coefficients)
  risk.score <- as.matrix((merged.data[,5:ncol(merged.data)] - colMeans(merged.data[,5:ncol(merged.data)]))) %*% coefs
  merged.data$risk_score <- risk.score
  group <- rep(NA, dim(merged.data)[1])
  for(i in 1:dim(merged.data)[1]){
    if(risk.score[i]>0){
      group[i] <- "High"
    }else{
      group[i] <- "Low"
    }
  }
  merged.data$group <- group
  
  return(merged.data)
}

mytheme.3 <- theme(legend.title = element_text(size = 14), 
                   legend.text = element_text(size = 14),
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14),
                   axis.line.x = element_line(size = 1),
                   axis.line.y = element_line(size = 1),
                   axis.title.x = element_text(size = 16),
                   axis.title.y = element_text(size = 16))

## --- macrophage ---
signatures <- c("C1QC","C1QA","APOE","MARCO","INHBA","IL1RN","CCL4",
                "NLRP3","EREG","IL1B","LYVE1","PLTP","SEPP1")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "Macrophage\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","macro.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","macro.surv.2.pdf"),p.surv.2,width = 4,height = 4)

## --- granulocyte ---
signatures <- c("FCGR3B", "IFITM2", "CXCR2", "S100A8", "SLC25A37","CXCR1","IL1R2")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 4, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "Granulocyte\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","gra.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","gra.surv.2.pdf"),p.surv.2,width = 4,height = 4)

## --- Mast cells ---
signatures <- c("KIT","CPA3","TPSAB1","HDC","GATA2","HPGDS","TPSD1","SLC18A2","MS4A2","IL1RL1","VWA5A")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "Mast\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","mast.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","mast.surv.2.pdf"),p.surv.2,width = 4,height = 4)

## --- CD14 Mono. ---
signatures <- c("CD14","FCN1","S100A8","S100A9","S100A12","VCAN","CD36")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "CD14_Mono.\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","CD14_Mono.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","CD14_Mono.surv.2.pdf"),p.surv.2,width = 4,height = 4)

## --- CD16 Mono. ---
signatures <- c("FCGR3A","LST1","LILRB2","IFITM2","SIGLEC10","CX3CR1","LILRB1","LILRA1","TCF7L2","MTSS1","RHOC")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "CD16_Mono.\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","CD16_Mono.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","CD16_Mono.surv.2.pdf"),p.surv.2,width = 4,height = 4)

## --- DCs ---
signatures <- c("GZMB","JCHAIN","MZB1","CLIC3","CXCL8","IL7R","CCR7","MMP7","IL32")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = signatures)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "DCs\nsignatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","DCs.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","DCs.surv.2.pdf"),p.surv.2,width = 4,height = 4)

rm(list = ls()[!grepl("myeloid",ls())]);gc()


## ----------- 4.Tissue Preference Estimation ---------------

rm(list = ls());gc()
# myeloid all
myeloid <- readRDS("input/mye.with.mast.cell.annot.rds")
chiq.N <- chisq.test(table(myeloid$CellTypes, myeloid$Patient)[,1])
chiq.T <- chisq.test(table(myeloid$CellTypes, myeloid$Patient)[,2])
Normal <- chiq.N$observed/chiq.N$expected
Tumor <- chiq.T$observed/chiq.T$expected
result <- rbind(Tumor, Normal)
result <- t(result)
wilcox.test(result)$p.value

dat <- reshape2::melt(result)
colnames(dat) <- c("CellTypes","Group","Roe")
p.Roe.NT <- ggplot(dat, aes(x=CellTypes,y = Roe,fill=Group)) + 
  geom_bar(stat = "identity",position = "dodge",color="black",size=1.3) +
  geom_text(aes(label=signif(Roe,3),size=6), position = position_dodge2(width = 0.9,preserve = 'single'),
            vjust = -0.20, hjust = 0.5) +
  scale_fill_manual(values = ggsci::pal_jama(alpha = 1)(2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,hjust = 1,vjust = 1,face = "bold"),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.3),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave(paste0("figures/celltypes/myeloid/TissuePreference/mye.Roe.pdf"),p.Roe.NT ,width = 8,height = 4)

# myeloid pdac (class: anatomic sites)
myeloid.pdac <- readRDS("input/mye.with.mast.pdac.sub.rds")
chiq.Met <- chisq.test(table(myeloid.pdac$CellTypes, myeloid.pdac$Class)[,1])
chiq.Para <- chisq.test(table(myeloid.pdac$CellTypes, myeloid.pdac$Class)[,2])
chiq.Pri <- chisq.test(table(myeloid.pdac$CellTypes, myeloid.pdac$Class)[,3])
Paracancerous <- chiq.Para$observed/chiq.Para$expected
Primary <- chiq.Pri$observed/chiq.Pri$expected
Metastatic <- chiq.Met$observed/chiq.Met$expected
result <- rbind(Paracancerous, Primary, Metastatic)
result <- t(result)
wilcox.test(result)$p.value

dat <- reshape2::melt(result)
colnames(dat) <- c("CellTypes","Group","Roe")
p.pdac.as <- ggplot(dat, aes(x=CellTypes,y = Roe,fill=Group)) + 
  geom_bar(stat = "identity",position = "dodge",color="black",size=1.3) +
  geom_text(aes(label=signif(Roe,3),size=6), position = position_dodge2(width = 0.9,preserve = 'single')) +
  scale_fill_manual(values = ggsci::pal_jama(alpha = 1)(3)) +
  theme_bw() +coord_flip() +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14,hjust = 1,vjust = 1,face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.3),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) 
ggsave(paste0("figures/celltypes/myeloid/TissuePreference/mye.pdac.anatomicSite.pdf"),p.pdac.as ,width = 6,height = 7)

# myeloid pdac (tissue: tissue, pbmc)
chiq.pbmc <- chisq.test(table(myeloid.pdac$CellTypes, myeloid.pdac$Tissue)[,1])
chiq.tiss <- chisq.test(table(myeloid.pdac$CellTypes, myeloid.pdac$Tissue)[,2])
Tissue <- chiq.tiss$observed/chiq.tiss$expected
PBMC <- chiq.pbmc$observed/chiq.pbmc$expected
result <- rbind(Tissue, PBMC)
result <- t(result)
wilcox.test(result)$p.value

dat <- reshape2::melt(result)
colnames(dat) <- c("CellTypes","Group","Roe")
p.pdac.ts <- ggplot(dat, aes(x=CellTypes,y = Roe,fill=Group)) + 
  geom_bar(stat = "identity",position = "dodge",color="black",size=1.3) +
  geom_text(aes(label=signif(Roe,3),size=6), position = position_dodge2(width = 0.9,preserve = 'single')) +
  scale_fill_manual(values = ggsci::pal_jama(alpha = 1)(2)) +
  theme_bw() + coord_flip() +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14,hjust = 1,vjust = 1,face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.3),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave(paste0("figures/celltypes/myeloid/TissuePreference/mye.pdac.tissue.pdf"),p.pdac.ts,width = 6,height = 6)

# myeloid pbmc (PDAC vs Healthy)
myeloid.pbmc <- readRDS("input/mye.with.mast.pbmc.sub.rds")
chiq.Normal <- chisq.test(table(myeloid.pbmc$CellTypes, myeloid.pbmc$Patient)[,1])
chiq.Tumor <- chisq.test(table(myeloid.pbmc$CellTypes, myeloid.pbmc$Patient)[,2])
PDAC <- chiq.Tumor$observed/chiq.Tumor$expected
Healthy <- chiq.Normal$observed/chiq.Normal$expected
result <- t(rbind(PDAC, Healthy))
wilcox.test(result)$p.value

dat <- reshape2::melt(result)
colnames(dat) <- c("CellTypes","Group","Roe")
p.pbmc.nt <- ggplot(dat, aes(x=CellTypes,y = Roe,fill=Group)) + 
  geom_bar(stat = "identity",position = "dodge",color="black",size=1.3) +
  geom_text(aes(label=signif(Roe,3),size=6), position = position_dodge2(width = 0.9,preserve = 'single')) +
  scale_fill_manual(values = ggsci::pal_jama(alpha = 1)(2)) +
  theme_bw() + coord_flip() +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14,hjust = 1,vjust = 1,face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.3),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
ggsave(paste0("figures/celltypes/myeloid/TissuePreference/mye.pbmc.tumor_vs_normal.pdf"),p.pbmc.nt,width = 6,height = 6)

rm(list = ls()[!grepl("myeloid",ls())]);gc()


## ----------- 5.Differential Gene Analysis ---------------

library(ggplot2)
library(ggrepel)

mytheme.4 <- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                   axis.text = element_text(size = 12),
                   legend.title = element_text(size = 12),
                   legend.text = element_text(size = 12)) 

#### ---->>> PDAC (tumor vs. para-tumor && tumor.tissue vs. paratumor.tissue && Tissue vs. PBMC) <<<----

myeloid.pdac.tissue <- subset(myeloid.pdac, cells = rownames(subset(myeloid.pdac@meta.data, Tissue == "Tissue")))
saveRDS(myeloid.pdac.tissue, file = "input/mye.with.mast.pdac.tissue.sub.rds")
mye.pdac.tumor <- subset(myeloid.pdac,cells=rownames(subset(myeloid.pdac@meta.data, State == "tumor")))

## --- macrophage ---
## 1.tumor vs para-tumor
Idents(myeloid.pdac) <- myeloid.pdac$CellTypes
mye.pdac.macro.tumor_vs_paratumor <- FindMarkers(myeloid.pdac, ident.1 = "tumor",group.by = "State",subset.ident = "Macrophage",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.macro.tumor_vs_paratumor <- mye.pdac.macro.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.macro.tumor_vs_paratumor)),]
mye.pdac.macro.tumor_vs_paratumor <- mye.pdac.macro.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.macro.tumor_vs_paratumor)),]
mye.pdac.macro.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.macro.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.macro.tumor_vs_paratumor <- mye.pdac.macro.tumor_vs_paratumor[order(mye.pdac.macro.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.macro.tumor_vs_paratumor,10),tail(mye.pdac.macro.tumor_vs_paratumor,10))

p.mye.pdac.macro.tumor_vs_paratumor <- ggplot(mye.pdac.macro.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.macro.tumor_vs_paratumor.pdf"),p.mye.pdac.macro.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tissue.macro.tumor_vs_paratumor <- FindMarkers(myeloid.pdac.tissue,ident.1 = "tumor",group.by = "State",subset.ident = "Macrophage")
# remove ribosomal and mitochondrial genes
mye.pdac.tissue.macro.tumor_vs_paratumor <- mye.pdac.tissue.macro.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tissue.macro.tumor_vs_paratumor)),]
mye.pdac.tissue.macro.tumor_vs_paratumor <- mye.pdac.tissue.macro.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tissue.macro.tumor_vs_paratumor)),]
mye.pdac.tissue.macro.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tissue.macro.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tissue.macro.tumor_vs_paratumor <- mye.pdac.tissue.macro.tumor_vs_paratumor[order(mye.pdac.tissue.macro.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tissue.macro.tumor_vs_paratumor,10), tail(mye.pdac.tissue.macro.tumor_vs_paratumor,10))
p.mye.pdac.tissue.macro.tumor_vs_paratumor <- ggplot(mye.pdac.tissue.macro.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.macro.tumor_vs_paratumor.pdf"),p.mye.pdac.tissue.macro.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.macro.tissue_vs_pbmc <- FindMarkers(mye.pdac.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "Macrophage",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.macro.tissue_vs_pbmc <- mye.pdac.macro.tissue_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.macro.tissue_vs_pbmc)),]
mye.pdac.macro.tissue_vs_pbmc <- mye.pdac.macro.tissue_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.macro.tissue_vs_pbmc)),]
mye.pdac.macro.tissue_vs_pbmc$Significance <- ifelse(mye.pdac.macro.tissue_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.macro.tissue_vs_pbmc <- mye.pdac.macro.tissue_vs_pbmc[order(mye.pdac.macro.tissue_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.macro.tissue_vs_pbmc,10), tail(mye.pdac.macro.tissue_vs_pbmc,10))
p.mye.pdac.macro.tissue_vs_pbmc <- ggplot(mye.pdac.macro.tissue_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.macro.tissue_vs_pbmc.pdf"),p.mye.pdac.macro.tissue_vs_pbmc,width = 7,height = 5)

save(mye.pdac.macro.tumor_vs_paratumor, 
     mye.pdac.tissue.macro.tumor_vs_paratumor,
     mye.pdac.macro.tissue_vs_pbmc, 
     file = "output/myeloid/DEgenes/mye.pdac.marco.DEgenes.rda")

rm(list = ls()[!grepl("myeloid",ls())]);gc()

## --- Mast cells ---
# 1.tumor vs para-tumor
Idents(myeloid.pdac) <- myeloid.pdac$CellTypes
mye.pdac.mast.tumor_vs_paratumor <- FindMarkers(myeloid.pdac, ident.1 = "tumor",group.by = "State",subset.ident = "Mast cells",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.mast.tumor_vs_paratumor <- mye.pdac.mast.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.mast.tumor_vs_paratumor)),]
mye.pdac.mast.tumor_vs_paratumor <- mye.pdac.mast.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.mast.tumor_vs_paratumor)),]
mye.pdac.mast.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.mast.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.mast.tumor_vs_paratumor <- mye.pdac.mast.tumor_vs_paratumor[order(mye.pdac.mast.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.mast.tumor_vs_paratumor,10),tail(mye.pdac.mast.tumor_vs_paratumor,10))
p.mye.pdac.mast.tumor_vs_paratumor <- ggplot(mye.pdac.mast.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.mast.tumor_vs_paratumor.pdf"),p.mye.pdac.mast.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tissue.mast.tumor_vs_paratumor <- FindMarkers(myeloid.pdac.tissue,ident.1 = "tumor",group.by = "State",subset.ident = "Mast cells")
# remove ribosomal and mitochondrial genes
mye.pdac.tissue.mast.tumor_vs_paratumor <- mye.pdac.tissue.mast.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tissue.mast.tumor_vs_paratumor)),]
mye.pdac.tissue.mast.tumor_vs_paratumor <- mye.pdac.tissue.mast.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tissue.mast.tumor_vs_paratumor)),]
mye.pdac.tissue.mast.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tissue.mast.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tissue.mast.tumor_vs_paratumor <- mye.pdac.tissue.mast.tumor_vs_paratumor[order(mye.pdac.tissue.mast.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tissue.mast.tumor_vs_paratumor,10),tail(mye.pdac.tissue.mast.tumor_vs_paratumor,10))
p.mye.pdac.tissue.mast.tumor_vs_paratumor <- ggplot(mye.pdac.tissue.mast.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.mast.tumor_vs_paratumor.pdf"),p.mye.pdac.tissue.mast.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.mast.tissue_vs_pbmc <- FindMarkers(mye.pdac.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "Mast cells",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.mast.tissue_vs_pbmc <- mye.pdac.mast.tissue_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.mast.tissue_vs_pbmc)),]
mye.pdac.mast.tissue_vs_pbmc <- mye.pdac.mast.tissue_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.mast.tissue_vs_pbmc)),]
mye.pdac.mast.tissue_vs_pbmc$Significance <- ifelse(mye.pdac.mast.tissue_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.mast.tissue_vs_pbmc <- mye.pdac.mast.tissue_vs_pbmc[order(mye.pdac.mast.tissue_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.mast.tissue_vs_pbmc,10), tail(mye.pdac.mast.tissue_vs_pbmc,10))
p.mye.pdac.mast.tissue_vs_pbmc <- ggplot(mye.pdac.mast.tissue_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.mast.tissue_vs_pbmc.pdf"),p.mye.pdac.mast.tissue_vs_pbmc,width = 7,height = 5)

save(mye.pdac.mast.tumor_vs_paratumor, 
     mye.pdac.tissue.mast.tumor_vs_paratumor,
     mye.pdac.mast.tissue_vs_pbmc, 
     file = "output/myeloid/DEgenes/mye.pdac.mast.DEgenes.rda")

rm(list = ls()[!grepl("myeloid",ls())]);gc()

## --- CD14 Mono. ---
## 1.tumor vs para-tumor
mye.pdac.CD14.mono.tumor_vs_paratumor <- FindMarkers(myeloid.pdac, ident.1 = "tumor",group.by = "State",subset.ident = "CD14 Mono.",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.CD14.mono.tumor_vs_paratumor <- mye.pdac.CD14.mono.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.CD14.mono.tumor_vs_paratumor)),]
mye.pdac.CD14.mono.tumor_vs_paratumor <- mye.pdac.CD14.mono.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.CD14.mono.tumor_vs_paratumor)),]
mye.pdac.CD14.mono.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.CD14.mono.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.CD14.mono.tumor_vs_paratumor <- mye.pdac.CD14.mono.tumor_vs_paratumor[order(mye.pdac.CD14.mono.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.CD14.mono.tumor_vs_paratumor,10),tail(mye.pdac.CD14.mono.tumor_vs_paratumor,10))

p.mye.pdac.CD14.mono.tumor_vs_paratumor <- ggplot(mye.pdac.CD14.mono.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.CD14.mono.tumor_vs_paratumor.pdf"),p.mye.pdac.CD14.mono.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor <- FindMarkers(myeloid.pdac.tissue,ident.1 = "tumor",group.by = "State",subset.ident = "CD14 Mono.")
# remove ribosomal and mitochondrial genes
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor)),]
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor)),]
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[order(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,10),tail(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,10))

p.mye.pdac.tissue.CD14.mono.tumor_vs_paratumor <- ggplot(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.CD14.mono.tumor_vs_paratumor.pdf"),p.mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.CD14.mono.tissue_vs_pbmc <- FindMarkers(mye.pdac.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "CD14 Mono.",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.CD14.mono.tissue_vs_pbmc <- mye.pdac.CD14.mono.tissue_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.CD14.mono.tissue_vs_pbmc)),]
mye.pdac.CD14.mono.tissue_vs_pbmc <- mye.pdac.CD14.mono.tissue_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.CD14.mono.tissue_vs_pbmc)),]
mye.pdac.CD14.mono.tissue_vs_pbmc$Significance <- ifelse(mye.pdac.CD14.mono.tissue_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.CD14.mono.tissue_vs_pbmc <- mye.pdac.CD14.mono.tissue_vs_pbmc[order(mye.pdac.CD14.mono.tissue_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.CD14.mono.tissue_vs_pbmc,10), tail(mye.pdac.CD14.mono.tissue_vs_pbmc,10))
p.mye.pdac.CD14.mono.tissue_vs_pbmc <- ggplot(mye.pdac.CD14.mono.tissue_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.CD14.mono.tissue_vs_pbmc.pdf"),p.mye.pdac.CD14.mono.tissue_vs_pbmc,width = 6,height = 5)

save(mye.pdac.CD14.mono.tumor_vs_paratumor, 
     mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,
     mye.pdac.CD14.mono.tissue_vs_pbmc, 
     file = "output/myeloid/DEgenes/mye.pdac.CD14.Mono.DEgenes.rda")

rm(list = ls()[!grepl("myeloid",ls())]);gc()

## --- CD16 Mono. ---
## 1.tumor vs para-tumor
mye.pdac.CD16.mono.tumor_vs_paratumor <- FindMarkers(myeloid.pdac, ident.1 = "tumor",group.by = "State",subset.ident = "CD16 Mono.",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.CD16.mono.tumor_vs_paratumor <- mye.pdac.CD16.mono.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.CD16.mono.tumor_vs_paratumor)),]
mye.pdac.CD16.mono.tumor_vs_paratumor <- mye.pdac.CD16.mono.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.CD16.mono.tumor_vs_paratumor)),]
mye.pdac.CD16.mono.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.CD16.mono.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.CD16.mono.tumor_vs_paratumor <- mye.pdac.CD16.mono.tumor_vs_paratumor[order(mye.pdac.CD16.mono.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.CD16.mono.tumor_vs_paratumor,10),tail(mye.pdac.CD16.mono.tumor_vs_paratumor,10))

p.mye.pdac.CD16.mono.tumor_vs_paratumor <- ggplot(mye.pdac.CD16.mono.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.CD16.mono.tumor_vs_paratumor.pdf"),p.mye.pdac.CD16.mono.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor <- FindMarkers(myeloid.pdac.tissue,ident.1 = "tumor",group.by = "State",subset.ident = "CD16 Mono.")
# remove ribosomal and mitochondrial genes
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor)),]
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor)),]
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor <- mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[order(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,10),tail(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,10))

p.mye.pdac.tissue.CD16.mono.tumor_vs_paratumor <- ggplot(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.CD16.mono.tumor_vs_paratumor.pdf"),p.mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC (PDAC)
mye.pdac.CD16.mono.tissue_vs_pbmc <- FindMarkers(mye.pdac.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "CD16 Mono.",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.CD16.mono.tissue_vs_pbmc <- mye.pdac.CD16.mono.tissue_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.CD16.mono.tissue_vs_pbmc)),]
mye.pdac.CD16.mono.tissue_vs_pbmc <- mye.pdac.CD16.mono.tissue_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.CD16.mono.tissue_vs_pbmc)),]
mye.pdac.CD16.mono.tissue_vs_pbmc$Significance <- ifelse(mye.pdac.CD16.mono.tissue_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.CD16.mono.tissue_vs_pbmc <- mye.pdac.CD16.mono.tissue_vs_pbmc[order(mye.pdac.CD16.mono.tissue_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.CD16.mono.tissue_vs_pbmc,10), tail(mye.pdac.CD16.mono.tissue_vs_pbmc,10))
p.mye.pdac.CD16.mono.tissue_vs_pbmc <- ggplot(mye.pdac.CD16.mono.tissue_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.CD16.mono.tissue_vs_pbmc.pdf"),p.mye.pdac.CD16.mono.tissue_vs_pbmc,width = 6,height = 5)

save(mye.pdac.CD16.mono.tumor_vs_paratumor, 
     mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,
     mye.pdac.CD16.mono.tissue_vs_pbmc, 
     file = "output/myeloid/DEgenes/mye.pdac.CD16.Mono.DEgenes.rda")

rm(list = ls()[!grepl("myeloid",ls())]);gc()

## --- DCs ---
## 1.tumor vs para-tumor
mye.pdac.DC.tumor_vs_paratumor <- FindMarkers(myeloid.pdac, ident.1 = "tumor",group.by = "State",subset.ident = "DCs",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.DC.tumor_vs_paratumor <- mye.pdac.DC.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.DC.tumor_vs_paratumor)),]
mye.pdac.DC.tumor_vs_paratumor <- mye.pdac.DC.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.DC.tumor_vs_paratumor)),]
mye.pdac.DC.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.DC.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.DC.tumor_vs_paratumor <- mye.pdac.DC.tumor_vs_paratumor[order(mye.pdac.DC.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.DC.tumor_vs_paratumor,10),tail(mye.pdac.DC.tumor_vs_paratumor,10))
p.mye.pdac.DC.tumor_vs_paratumor <- ggplot(mye.pdac.DC.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.DC.tumor_vs_paratumor.pdf"),p.mye.pdac.DC.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tissue.DC.tumor_vs_paratumor <- FindMarkers(myeloid.pdac.tissue,ident.1 = "tumor",group.by = "State",subset.ident = "DCs")
# remove ribosomal and mitochondrial genes
mye.pdac.tissue.DC.tumor_vs_paratumor <- mye.pdac.tissue.DC.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tissue.DC.tumor_vs_paratumor)),]
mye.pdac.tissue.DC.tumor_vs_paratumor <- mye.pdac.tissue.DC.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tissue.DC.tumor_vs_paratumor)),]
mye.pdac.tissue.DC.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tissue.DC.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tissue.DC.tumor_vs_paratumor <- mye.pdac.tissue.DC.tumor_vs_paratumor[order(mye.pdac.tissue.DC.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tissue.DC.tumor_vs_paratumor,10),tail(mye.pdac.tissue.DC.tumor_vs_paratumor,10))
p.mye.pdac.tissue.DC.tumor_vs_paratumor <- ggplot(mye.pdac.tissue.DC.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.DC.tumor_vs_paratumor.pdf"),p.mye.pdac.tissue.DC.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC (PDAC)
mye.pdac.DC.tissue_vs_pbmc <- FindMarkers(mye.pdac.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "DCs",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.DC.tissue_vs_pbmc <- mye.pdac.DC.tissue_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.DC.tissue_vs_pbmc)),]
mye.pdac.DC.tissue_vs_pbmc <- mye.pdac.DC.tissue_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.DC.tissue_vs_pbmc)),]
mye.pdac.DC.tissue_vs_pbmc$Significance <- ifelse(mye.pdac.DC.tissue_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.DC.tissue_vs_pbmc <- mye.pdac.DC.tissue_vs_pbmc[order(mye.pdac.DC.tissue_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.DC.tissue_vs_pbmc,10), tail(mye.pdac.DC.tissue_vs_pbmc,10))
p.mye.pdac.DC.tissue_vs_pbmc <- ggplot(mye.pdac.DC.tissue_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.DC.tissue_vs_pbmc.pdf"),p.mye.pdac.DC.tissue_vs_pbmc,width = 6,height = 5)

save(mye.pdac.DC.tumor_vs_paratumor, 
     mye.pdac.tissue.DC.tumor_vs_paratumor,
     mye.pdac.DC.tissue_vs_pbmc, 
     file = "output/myeloid/DEgenes/mye.pdac.DC.DEgenes.rda")

rm(list = ls());gc()

### ---------------->>>>>>>>>>>>> Joint Plot <<<<<<<<<<<<<<<----------------------
load("output/myeloid/DEgenes/mye.pdac.marco.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pdac.mast.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pdac.CD14.Mono.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pdac.CD16.Mono.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pdac.DC.DEgenes.rda")

## 1.tumor vs para-tumor
mye.pdac.macro.tumor_vs_paratumor$celltypes <- "Macrophage"
mye.pdac.macro.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.macro.tumor_vs_paratumor,10)),
                                            rep("",nrow(mye.pdac.macro.tumor_vs_paratumor)-20),
                                            rownames(tail(mye.pdac.macro.tumor_vs_paratumor,10)))
mye.pdac.mast.tumor_vs_paratumor$celltypes <- "Mast cells"
mye.pdac.mast.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.mast.tumor_vs_paratumor,10)),
                                           rep("",nrow(mye.pdac.mast.tumor_vs_paratumor)-20),
                                           rownames(tail(mye.pdac.mast.tumor_vs_paratumor,10)))
mye.pdac.CD14.mono.tumor_vs_paratumor$celltypes <- "CD14 Mono."
mye.pdac.CD14.mono.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.CD14.mono.tumor_vs_paratumor,10)),
                                                rep("",nrow(mye.pdac.CD14.mono.tumor_vs_paratumor)-20),
                                                rownames(tail(mye.pdac.CD14.mono.tumor_vs_paratumor,10)))
mye.pdac.CD16.mono.tumor_vs_paratumor$celltypes <- "CD16 Mono."
mye.pdac.CD16.mono.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.CD16.mono.tumor_vs_paratumor,10)),
                                                rep("",nrow(mye.pdac.CD16.mono.tumor_vs_paratumor)-20),
                                                rownames(tail(mye.pdac.CD16.mono.tumor_vs_paratumor,10)))
mye.pdac.DC.tumor_vs_paratumor$celltypes <- "DCs"
mye.pdac.DC.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.DC.tumor_vs_paratumor,10)),
                                         rep("",nrow(mye.pdac.DC.tumor_vs_paratumor)-20),
                                         rownames(tail(mye.pdac.DC.tumor_vs_paratumor,10)))

mye.pdac.tumor_vs_paratumor.DEgenes <- rbind(mye.pdac.macro.tumor_vs_paratumor,
                                             mye.pdac.mast.tumor_vs_paratumor,
                                             mye.pdac.CD14.mono.tumor_vs_paratumor,
                                             mye.pdac.CD16.mono.tumor_vs_paratumor,
                                             mye.pdac.DC.tumor_vs_paratumor)
save(mye.pdac.tumor_vs_paratumor.DEgenes,file = "output/myeloid/DEgenes/mye.pdac.tumor_vs_paratumor.DEgenes.rda")

mye.pdac.tumor_vs_paratumor.DEgenes$bar.1 <- 0
mye.pdac.tumor_vs_paratumor.DEgenes$bar.2 <- 0
for(i in 1:length(unique(mye.pdac.tumor_vs_paratumor.DEgenes$celltypes))){
  tmp <- which(mye.pdac.tumor_vs_paratumor.DEgenes$celltypes==unique(mye.pdac.tumor_vs_paratumor.DEgenes$celltypes)[i])[1]
  mye.pdac.tumor_vs_paratumor.DEgenes$bar.1[tmp] <- 0.2
  mye.pdac.tumor_vs_paratumor.DEgenes$bar.2[tmp] <- -0.2
}
mye.pdac.tumor_vs_paratumor.DEgenes$celltypes <- factor(mye.pdac.tumor_vs_paratumor.DEgenes$celltypes,levels = c("Macrophage","Mast cells","CD14 Mono.","CD16 Mono.","DCs"))
p <- ggplot(mye.pdac.tumor_vs_paratumor.DEgenes, aes(celltypes,avg_logFC, color=Significance)) +
  geom_point(aes(col = Significance)) +
  geom_jitter() + scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(aes(celltypes, avg_logFC, label = gene),color="black") +
  scale_y_continuous(limits = c(-3,3)) +
  theme_classic() +
  ylab("average log2FC") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size=5))) 

p2 <- p + geom_bar(aes(celltypes,bar.1,fill=celltypes),stat = 'identity',color=NA) +
  geom_bar(aes(celltypes,bar.2,fill=celltypes),stat = 'identity',color=NA) +
  scale_fill_manual(values = c("#DF78C1","#00B050","#DB4840","#FF7E07","#7030A0")) +
  theme(legend.position = "none")
  

ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tumor_vs_paratumor.DEgenes.pdf"),p2,width = 18,height = 12)

## 2.pdac.tissue_vs_pbmc
mye.pdac.macro.tissue_vs_pbmc$celltypes <- "Macrophage"
mye.pdac.macro.tissue_vs_pbmc$gene <- c(rownames(head(mye.pdac.macro.tissue_vs_pbmc,10)),
                                        rep("",nrow(mye.pdac.macro.tissue_vs_pbmc)-20),
                                        rownames(tail(mye.pdac.macro.tissue_vs_pbmc,10)))
mye.pdac.mast.tissue_vs_pbmc$celltypes <- "Mast cells"
mye.pdac.mast.tissue_vs_pbmc$gene <- c(rownames(head(mye.pdac.mast.tissue_vs_pbmc,10)),
                                       rep("",nrow(mye.pdac.mast.tissue_vs_pbmc)-20),
                                       rownames(tail(mye.pdac.mast.tissue_vs_pbmc,10)))
mye.pdac.CD14.mono.tissue_vs_pbmc$celltypes <- "CD14 Mono."
mye.pdac.CD14.mono.tissue_vs_pbmc$gene <- c(rownames(head(mye.pdac.CD14.mono.tissue_vs_pbmc,10)),
                                            rep("",nrow(mye.pdac.CD14.mono.tissue_vs_pbmc)-20),
                                            rownames(tail(mye.pdac.CD14.mono.tissue_vs_pbmc,10)))
mye.pdac.CD16.mono.tissue_vs_pbmc$celltypes <- "CD16 Mono."
mye.pdac.CD16.mono.tissue_vs_pbmc$gene <- c(rownames(head(mye.pdac.CD16.mono.tissue_vs_pbmc,10)),
                                            rep("",nrow(mye.pdac.CD16.mono.tissue_vs_pbmc)-20),
                                            rownames(tail(mye.pdac.CD16.mono.tissue_vs_pbmc,10)))
mye.pdac.DC.tissue_vs_pbmc$celltypes <- "DCs"
mye.pdac.DC.tissue_vs_pbmc$gene <- c(rownames(head(mye.pdac.DC.tissue_vs_pbmc,10)),
                                     rep("",nrow(mye.pdac.DC.tissue_vs_pbmc)-20),
                                     rownames(tail(mye.pdac.DC.tissue_vs_pbmc,10)))
mye.pdac.tissue_vs_pbmc.DEgenes <- rbind(mye.pdac.macro.tissue_vs_pbmc,
                                         mye.pdac.mast.tissue_vs_pbmc,
                                         mye.pdac.CD14.mono.tissue_vs_pbmc,
                                         mye.pdac.CD16.mono.tissue_vs_pbmc,
                                         mye.pdac.DC.tissue_vs_pbmc)
save(mye.pdac.tissue_vs_pbmc.DEgenes,file = "output/myeloid/DEgenes/mye.pdac.tissue_vs_pbmc.DEgenes.rda")

mye.pdac.tissue_vs_pbmc.DEgenes$bar.1 <- 0
mye.pdac.tissue_vs_pbmc.DEgenes$bar.2 <- 0
for(i in 1:length(unique(mye.pdac.tissue_vs_pbmc.DEgenes$celltypes))){
  tmp <- which(mye.pdac.tissue_vs_pbmc.DEgenes$celltypes==unique(mye.pdac.tissue_vs_pbmc.DEgenes$celltypes)[i])[1]
  mye.pdac.tissue_vs_pbmc.DEgenes$bar.1[tmp] <- 0.2
  mye.pdac.tissue_vs_pbmc.DEgenes$bar.2[tmp] <- -0.2
}
mye.pdac.tissue_vs_pbmc.DEgenes$celltypes <- factor(mye.pdac.tissue_vs_pbmc.DEgenes$celltypes,levels = c("Macrophage","Mast cells","CD14 Mono.","CD16 Mono.","DCs"))
p <- ggplot(mye.pdac.tissue_vs_pbmc.DEgenes, aes(celltypes,avg_logFC, color=Significance)) +
  geom_point(aes(col = Significance)) +
  geom_jitter() + scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(aes(celltypes, avg_logFC, label = gene),color="black") +
  scale_y_continuous(limits = c(-3,3)) +
  theme_classic() +
  ylab("average log2FC") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size=5))) 

p2 <- p + geom_bar(aes(celltypes,bar.1,fill=celltypes),stat = 'identity',color=NA) +
  geom_bar(aes(celltypes,bar.2,fill=celltypes),stat = 'identity',color=NA) +
  scale_fill_manual(values = c("#DF78C1","#00B050","#DB4840","#FF7E07","#7030A0")) +
  theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue_vs_pbmc.DEgenes.pdf"),p2,width = 18,height = 12)

## 3.pdac.tissue.tumor_vs_paratumor
mye.pdac.tissue.macro.tumor_vs_paratumor$celltypes <- "Macrophage"
mye.pdac.tissue.macro.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.tissue.macro.tumor_vs_paratumor,10)),
                                                   rep("",nrow(mye.pdac.tissue.macro.tumor_vs_paratumor)-20),
                                                   rownames(tail(mye.pdac.tissue.macro.tumor_vs_paratumor,10)))
mye.pdac.tissue.mast.tumor_vs_paratumor$celltypes <- "Mast cells"
mye.pdac.tissue.mast.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.tissue.mast.tumor_vs_paratumor,10)),
                                                  rep("",nrow(mye.pdac.tissue.mast.tumor_vs_paratumor)-20),
                                                  rownames(tail(mye.pdac.tissue.mast.tumor_vs_paratumor,10)))
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$celltypes <- "CD14 Mono."
mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,10)),
                                                       rep("",nrow(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor)-20),
                                                       rownames(tail(mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,10)))
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$celltypes <- "CD16 Mono."
mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,10)),
                                                       rep("",nrow(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor)-20),
                                                       rownames(tail(mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,10)))
mye.pdac.tissue.DC.tumor_vs_paratumor$celltypes <- "DCs"
mye.pdac.tissue.DC.tumor_vs_paratumor$gene <- c(rownames(head(mye.pdac.tissue.DC.tumor_vs_paratumor,10)),
                                                rep("",nrow(mye.pdac.tissue.DC.tumor_vs_paratumor)-20),
                                                rownames(tail(mye.pdac.tissue.DC.tumor_vs_paratumor,10)))
mye.pdac.tissue.tumor_vs_paratumor.DEgenes <- rbind(mye.pdac.tissue.macro.tumor_vs_paratumor,
                                                    mye.pdac.tissue.mast.tumor_vs_paratumor,
                                                    mye.pdac.tissue.CD14.mono.tumor_vs_paratumor,
                                                    mye.pdac.tissue.CD16.mono.tumor_vs_paratumor,
                                                    mye.pdac.tissue.DC.tumor_vs_paratumor)
save(mye.pdac.tissue.tumor_vs_paratumor.DEgenes,file = "output/myeloid/DEgenes/mye.pdac.tissue.tumor_vs_paratumor.DEgenes.rda")

mye.pdac.tissue.tumor_vs_paratumor.DEgenes$bar.1 <- 0
mye.pdac.tissue.tumor_vs_paratumor.DEgenes$bar.2 <- 0
for(i in 1:length(unique(mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes))){
  tmp <- which(mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes==unique(mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes)[i])[1]
  mye.pdac.tissue.tumor_vs_paratumor.DEgenes$bar.1[tmp] <- 0.2
  mye.pdac.tissue.tumor_vs_paratumor.DEgenes$bar.2[tmp] <- -0.2
}
mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes <- factor(mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes,levels = c("Macrophage","Mast cells","CD14 Mono.","CD16 Mono.","DCs"))
p <- ggplot(mye.pdac.tissue.tumor_vs_paratumor.DEgenes, aes(celltypes,avg_logFC, color=Significance)) +
  geom_point(aes(col = Significance)) +
  geom_jitter() + scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(aes(celltypes, avg_logFC, label = gene),color="black") +
  scale_y_continuous(limits = c(-3,3)) +
  theme_classic() +
  ylab("average log2FC") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size=5))) 

p2 <- p + geom_bar(aes(celltypes,bar.1,fill=celltypes),stat = 'identity',color=NA) +
  geom_bar(aes(celltypes,bar.2,fill=celltypes),stat = 'identity',color=NA) +
  scale_fill_manual(values = c("#DF78C1","#00B050","#DB4840","#FF7E07","#7030A0")) +
  theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pdac.tissue.tumor_vs_paratumor.DEgenes.pdf"),p2,width = 18,height = 12)

#### ---------------->>>>>>>>>>>>> PBMC (PDAC vs Normal) <<<<<<<<<<<<<<<-----------------------------

Idents(myeloid.pbmc) <- myeloid.pbmc$CellTypes
# macrophage
mye.pbmc.macro.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "Macrophage",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.macro.pdac_vs_normal)),]
mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.macro.pdac_vs_normal)),]
mye.pbmc.macro.pdac_vs_normal$Significance <- ifelse(mye.pbmc.macro.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[order(mye.pbmc.macro.pdac_vs_normal[,2],decreasing = T),]
top10 <- rbind(head(mye.pbmc.macro.pdac_vs_normal,10), tail(mye.pbmc.macro.pdac_vs_normal,10))
p.mye.pbmc.macro.pdac_vs_normal <- ggplot(mye.pbmc.macro.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.macro.pdac_vs_normal.pdf"),p.mye.pbmc.macro.pdac_vs_normal,width = 7,height = 5)
save(mye.pbmc.macro.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.macro.DEgenes.rda")

# Mast cells
mye.pbmc.mast.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "Mast cells", logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.mast.pdac_vs_normal)),]
mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.mast.pdac_vs_normal)),]
mye.pbmc.mast.pdac_vs_normal$Significance <- ifelse(mye.pbmc.mast.pdac_vs_normal$p_val < 0.05, TRUE, FALSE)
mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[order(mye.pbmc.mast.pdac_vs_normal[,2],decreasing = T),]
top10 <- rbind(head(mye.pbmc.mast.pdac_vs_normal,10), tail(mye.pbmc.mast.pdac_vs_normal,10))
p.mye.pbmc.mast.pdac_vs_normal <- ggplot(mye.pbmc.mast.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.mast.pdac_vs_normal.pdf"),p.mye.pbmc.mast.pdac_vs_normal,width = 7,height = 5)
save(mye.pbmc.mast.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.mast.DEgenes.rda")

# CD14 Mono.
mye.pbmc.CD14.mono.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "CD14 Mono.", logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.CD14.mono.pdac_vs_normal)),]
mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.CD14.mono.pdac_vs_normal)),]
mye.pbmc.CD14.mono.pdac_vs_normal$Significance <- ifelse(mye.pbmc.CD14.mono.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[order(mye.pbmc.CD14.mono.pdac_vs_normal[,2],decreasing = T),]
top10 <- rbind(head(mye.pbmc.CD14.mono.pdac_vs_normal,10), tail(mye.pbmc.CD14.mono.pdac_vs_normal,10))
p.mye.pbmc.CD14.mono.pdac_vs_normal <- ggplot(mye.pbmc.CD14.mono.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  scale_y_continuous(limits = c(0,50)) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.CD14.mono.pdac_vs_normal.pdf"),p.mye.pbmc.CD14.mono.pdac_vs_normal,width = 7,height = 5)
save(mye.pbmc.CD14.mono.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.CD14.mono.DEgenes.rda")

# CD16 Mono.
mye.pbmc.CD16.mono.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "CD16 Mono.", logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pbmc.CD16.mono.pdac_vs_normal <- mye.pbmc.CD16.mono.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.CD16.mono.pdac_vs_normal)),]
mye.pbmc.CD16.mono.pdac_vs_normal <- mye.pbmc.CD16.mono.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.CD16.mono.pdac_vs_normal)),]
mye.pbmc.CD16.mono.pdac_vs_normal$Significance <- ifelse(mye.pbmc.CD16.mono.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
mye.pbmc.CD16.mono.pdac_vs_normal <- mye.pbmc.CD16.mono.pdac_vs_normal[order(mye.pbmc.CD16.mono.pdac_vs_normal[,2],decreasing = T),]
top10 <- rbind(head(mye.pbmc.CD16.mono.pdac_vs_normal,10), tail(mye.pbmc.CD16.mono.pdac_vs_normal,10))
p.mye.pbmc.CD16.mono.pdac_vs_normal <- ggplot(mye.pbmc.CD16.mono.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.CD16.mono.pdac_vs_normal.pdf"),p.mye.pbmc.CD16.mono.pdac_vs_normal,width = 7,height = 5)
save(mye.pbmc.CD16.mono.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.CD16.mono.DEgenes.rda")

# DCs
mye.pbmc.DC.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "DCs", logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pbmc.DC.pdac_vs_normal <- mye.pbmc.DC.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.DC.pdac_vs_normal)),]
mye.pbmc.DC.pdac_vs_normal <- mye.pbmc.DC.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.DC.pdac_vs_normal)),]
mye.pbmc.DC.pdac_vs_normal$Significance <- ifelse(mye.pbmc.DC.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
mye.pbmc.DC.pdac_vs_normal <- mye.pbmc.DC.pdac_vs_normal[order(mye.pbmc.DC.pdac_vs_normal[,2],decreasing = T),]
top10 <- rbind(head(mye.pbmc.DC.pdac_vs_normal,10), tail(mye.pbmc.DC.pdac_vs_normal,10))
p.mye.pbmc.DC.pdac_vs_normal <- ggplot(mye.pbmc.DC.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme.4
ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.DC.pdac_vs_normal.pdf"),p.mye.pbmc.DC.pdac_vs_normal,width = 7,height = 5)
save(mye.pbmc.DC.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.DC.DEgenes.rda")

# ---------- Output DE genes ------------------

# 1.mye.pdac.tumor_vs_paratumor
load("output/myeloid/DEgenes/mye.pdac.tumor_vs_paratumor.DEgenes.rda")
# macro
macro <- mye.pdac.tumor_vs_paratumor.DEgenes[mye.pdac.tumor_vs_paratumor.DEgenes$celltypes=="Macrophage",]
macro.up <- macro[macro$avg_logFC > 0,]
write(rownames(macro.up),file = "output/myeloid/DEgenes/txt/mye.pdac.macro.tumor_vs_paratumor.up.txt")
macro.down <- macro[macro$avg_logFC < 0,]
write(rownames(macro.down),file = "output/myeloid/DEgenes/txt/mye.pdac.macro.tumor_vs_paratumor.down.txt")
# mast
mast <- mye.pdac.tumor_vs_paratumor.DEgenes[mye.pdac.tumor_vs_paratumor.DEgenes$celltypes=="Mast cells",]
mast.up <- mast[mast$avg_logFC > 0,]
write(rownames(mast.up),file = "output/myeloid/DEgenes/txt/mye.pdac.mast.tumor_vs_paratumor.up.txt")
mast.down <- mast[mast$avg_logFC < 0,]
write(rownames(mast.down),file = "output/myeloid/DEgenes/txt/mye.pdac.mast.tumor_vs_paratumor.down.txt")
# CD14 Mono.
CD14Mono <- mye.pdac.tumor_vs_paratumor.DEgenes[mye.pdac.tumor_vs_paratumor.DEgenes$celltypes=="CD14 Mono.",]
CD14Mono.up <- CD14Mono[CD14Mono$avg_logFC > 0,]
write(rownames(CD14Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.CD14Mono.tumor_vs_paratumor.up.txt")
CD14Mono.down <- CD14Mono[CD14Mono$avg_logFC < 0,]
write(rownames(CD14Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.CD14Mono.tumor_vs_paratumor.down.txt")
# CD16 Mono.
CD16Mono <- mye.pdac.tumor_vs_paratumor.DEgenes[mye.pdac.tumor_vs_paratumor.DEgenes$celltypes=="CD16 Mono.",]
CD16Mono.up <- CD16Mono[CD16Mono$avg_logFC > 0,]
write(rownames(CD16Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.CD16Mono.tumor_vs_paratumor.up.txt")
CD16Mono.down <- CD16Mono[CD16Mono$avg_logFC < 0,]
write(rownames(CD16Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.CD16Mono.tumor_vs_paratumor.down.txt")
# DCs
DC <- mye.pdac.tumor_vs_paratumor.DEgenes[mye.pdac.tumor_vs_paratumor.DEgenes$celltypes=="DCs",]
DC.up <- DC[DC$avg_logFC > 0,]
write(rownames(DC.up),file = "output/myeloid/DEgenes/txt/mye.pdac.DC.tumor_vs_paratumor.up.txt")
DC.down <- DC[DC$avg_logFC < 0,]
write(rownames(DC.down),file = "output/myeloid/DEgenes/txt/mye.pdac.DC.tumor_vs_paratumor.down.txt")

rm(list = ls());gc()

# 2.mye.pdac.tumor_vs_paratumor
load("output/myeloid/DEgenes/mye.pdac.tissue.tumor_vs_paratumor.DEgenes.rda")
# macro
macro <- mye.pdac.tissue.tumor_vs_paratumor.DEgenes[mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes=="Macrophage",]
macro.up <- macro[macro$avg_logFC > 0,]
write(rownames(macro.up),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.macro.tumor_vs_paratumor.up.txt")
macro.down <- macro[macro$avg_logFC < 0,]
write(rownames(macro.down),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.macro.tumor_vs_paratumor.down.txt")
# mast
mast <- mye.pdac.tissue.tumor_vs_paratumor.DEgenes[mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes=="Mast cells",]
mast.up <- mast[mast$avg_logFC > 0,]
write(rownames(mast.up),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.mast.tumor_vs_paratumor.up.txt")
mast.down <- mast[mast$avg_logFC < 0,]
write(rownames(mast.down),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.mast.tumor_vs_paratumor.down.txt")
# CD14 Mono.
CD14Mono <- mye.pdac.tissue.tumor_vs_paratumor.DEgenes[mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes=="CD14 Mono.",]
CD14Mono.up <- CD14Mono[CD14Mono$avg_logFC > 0,]
write(rownames(CD14Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.CD14Mono.tumor_vs_paratumor.up.txt")
CD14Mono.down <- CD14Mono[CD14Mono$avg_logFC < 0,]
write(rownames(CD14Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.CD14Mono.tumor_vs_paratumor.down.txt")
# CD16 Mono.
CD16Mono <- mye.pdac.tissue.tumor_vs_paratumor.DEgenes[mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes=="CD16 Mono.",]
CD16Mono.up <- CD16Mono[CD16Mono$avg_logFC > 0,]
write(rownames(CD16Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.CD16Mono.tumor_vs_paratumor.up.txt")
CD16Mono.down <- CD16Mono[CD16Mono$avg_logFC < 0,]
write(rownames(CD16Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.CD16Mono.tumor_vs_paratumor.down.txt")
# DCs
DC <- mye.pdac.tissue.tumor_vs_paratumor.DEgenes[mye.pdac.tissue.tumor_vs_paratumor.DEgenes$celltypes=="DCs",]
DC.up <- DC[DC$avg_logFC > 0,]
write(rownames(DC.up),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.DC.tumor_vs_paratumor.up.txt")
DC.down <- DC[DC$avg_logFC < 0,]
write(rownames(DC.down),file = "output/myeloid/DEgenes/txt/mye.pdac.tissue.DC.tumor_vs_paratumor.down.txt")

rm(list = ls());gc()

# 3.mye.pdac.tumor_vs_paratumor
load("output/myeloid/DEgenes/mye.pdac.tissue_vs_pbmc.DEgenes.rda")
# macro
macro <- mye.pdac.tissue_vs_pbmc.DEgenes[mye.pdac.tissue_vs_pbmc.DEgenes$celltypes=="Macrophage",]
macro.up <- macro[macro$avg_logFC > 0,]
write(rownames(macro.up),file = "output/myeloid/DEgenes/txt/mye.pdac.macro.tissue_vs_pbmc.up.txt")
macro.down <- macro[macro$avg_logFC < 0,]
write(rownames(macro.down),file = "output/myeloid/DEgenes/txt/mye.pdac.macro.tissue_vs_pbmc.down.txt")
# mast
mast <- mye.pdac.tissue_vs_pbmc.DEgenes[mye.pdac.tissue_vs_pbmc.DEgenes$celltypes=="Mast cells",]
mast.up <- mast[mast$avg_logFC > 0,]
write(rownames(mast.up),file = "output/myeloid/DEgenes/txt/mye.pdac.mast.tissue_vs_pbmc.up.txt")
mast.down <- mast[mast$avg_logFC < 0,]
write(rownames(mast.down),file = "output/myeloid/DEgenes/txt/mye.pdac.mast.tissue_vs_pbmc.down.txt")
# CD14 Mono.
CD14Mono <- mye.pdac.tissue_vs_pbmc.DEgenes[mye.pdac.tissue_vs_pbmc.DEgenes$celltypes=="CD14 Mono.",]
CD14Mono.up <- CD14Mono[CD14Mono$avg_logFC > 0,]
write(rownames(CD14Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.CD14Mono.tissue_vs_pbmc.up.txt")
CD14Mono.down <- CD14Mono[CD14Mono$avg_logFC < 0,]
write(rownames(CD14Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.CD14Mono.tissue_vs_pbmc.down.txt")
# CD16 Mono.
CD16Mono <- mye.pdac.tissue_vs_pbmc.DEgenes[mye.pdac.tissue_vs_pbmc.DEgenes$celltypes=="CD16 Mono.",]
CD16Mono.up <- CD16Mono[CD16Mono$avg_logFC > 0,]
write(rownames(CD16Mono.up),file = "output/myeloid/DEgenes/txt/mye.pdac.CD16Mono.tissue_vs_pbmc.up.txt")
CD16Mono.down <- CD16Mono[CD16Mono$avg_logFC < 0,]
write(rownames(CD16Mono.down),file = "output/myeloid/DEgenes/txt/mye.pdac.CD16Mono.tissue_vs_pbmc.down.txt")
# DCs
DC <- mye.pdac.tissue_vs_pbmc.DEgenes[mye.pdac.tissue_vs_pbmc.DEgenes$celltypes=="DCs",]
DC.up <- DC[DC$avg_logFC > 0,]
write(rownames(DC.up),file = "output/myeloid/DEgenes/txt/mye.pdac.DC.tissue_vs_pbmc.up.txt")
DC.down <- DC[DC$avg_logFC < 0,]
write(rownames(DC.down),file = "output/myeloid/DEgenes/txt/mye.pdac.DC.tissue_vs_pbmc.down.txt")

rm(list = ls());gc()


### ---------------->>>>>>>>>>>>> Joint Plot <<<<<<<<<<<<<<<----------------------
load("output/myeloid/DEgenes/mye.pbmc.macro.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pbmc.mast.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pbmc.CD14.mono.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pbmc.CD16.mono.DEgenes.rda")
load("output/myeloid/DEgenes/mye.pbmc.DC.DEgenes.rda")

# 1.pbmc.pdac_vs_normal
mye.pbmc.macro.pdac_vs_normal$celltypes <- "Macrophage"
mye.pbmc.macro.pdac_vs_normal$gene <- c(rownames(head(mye.pbmc.macro.pdac_vs_normal,10)),
                                        rep("",nrow(mye.pbmc.macro.pdac_vs_normal)-20),
                                        rownames(tail(mye.pbmc.macro.pdac_vs_normal,10)))
mye.pbmc.mast.pdac_vs_normal$celltypes <- "Mast cells"
mye.pbmc.mast.pdac_vs_normal$gene <- c(rownames(head(mye.pbmc.mast.pdac_vs_normal,10)),
                                        rep("",nrow(mye.pbmc.mast.pdac_vs_normal)-20),
                                        rownames(tail(mye.pbmc.mast.pdac_vs_normal,10)))
mye.pbmc.CD14.mono.pdac_vs_normal$celltypes <- "CD14 Mono."
mye.pbmc.CD14.mono.pdac_vs_normal$gene <- c(rownames(head(mye.pbmc.CD14.mono.pdac_vs_normal,10)),
                                            rep("",nrow(mye.pbmc.CD14.mono.pdac_vs_normal)-20),
                                            rownames(tail(mye.pbmc.CD14.mono.pdac_vs_normal,10)))
mye.pbmc.CD16.mono.pdac_vs_normal$celltypes <- "CD16 Mono."
mye.pbmc.CD16.mono.pdac_vs_normal$gene <- c(rownames(head(mye.pbmc.CD16.mono.pdac_vs_normal,10)),
                                            rep("",nrow(mye.pbmc.CD16.mono.pdac_vs_normal)-20),
                                            rownames(tail(mye.pbmc.CD16.mono.pdac_vs_normal,10)))
mye.pbmc.DC.pdac_vs_normal$celltypes <- "DCs"
mye.pbmc.DC.pdac_vs_normal$gene <- c(rownames(head(mye.pbmc.DC.pdac_vs_normal,10)),
                                     rep("",nrow(mye.pbmc.DC.pdac_vs_normal)-20),
                                     rownames(tail(mye.pbmc.DC.pdac_vs_normal,10)))
mye.pbmc.pdac_vs_normal.DEgenes <- rbind(mye.pbmc.macro.pdac_vs_normal,
                                         mye.pbmc.mast.pdac_vs_normal,
                                         mye.pbmc.CD14.mono.pdac_vs_normal,
                                         mye.pbmc.CD16.mono.pdac_vs_normal,
                                         mye.pbmc.DC.pdac_vs_normal)

save(mye.pbmc.pdac_vs_normal.DEgenes,file = "output/myeloid/DEgenes/mye.pbmc.pdac_vs_normal.DEgenes.rda")

mye.pbmc.pdac_vs_normal.DEgenes$bar.1 <- 0
mye.pbmc.pdac_vs_normal.DEgenes$bar.2 <- 0
for(i in 1:length(unique(mye.pbmc.pdac_vs_normal.DEgenes$celltypes))){
  tmp <- which(mye.pbmc.pdac_vs_normal.DEgenes$celltypes==unique(mye.pbmc.pdac_vs_normal.DEgenes$celltypes)[i])[1]
  mye.pbmc.pdac_vs_normal.DEgenes$bar.1[tmp] <- 0.2
  mye.pbmc.pdac_vs_normal.DEgenes$bar.2[tmp] <- -0.2
}
mye.pbmc.pdac_vs_normal.DEgenes$celltypes <- factor(mye.pbmc.pdac_vs_normal.DEgenes$celltypes,levels = c("Macrophage","Mast cells","CD14 Mono.","CD16 Mono.","DCs"))
p <- ggplot(mye.pbmc.pdac_vs_normal.DEgenes, aes(celltypes,avg_logFC, color=Significance)) +
  geom_point(aes(col = Significance)) +
  geom_jitter() + scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(aes(celltypes, avg_logFC, label = gene),color="black") +
  scale_y_continuous(limits = c(-1.5,2)) +
  theme_classic() +
  ylab("average log2FC") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size=5))) 

p2 <- p + geom_bar(aes(celltypes,bar.1,fill=celltypes),stat = 'identity',color=NA) +
  geom_bar(aes(celltypes,bar.2,fill=celltypes),stat = 'identity',color=NA) +
  scale_fill_manual(values = c("#DF78C1","#00B050","#DB4840","#FF7E07","#7030A0")) +
  theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.pdac_vs_normal.DEgenes.pdf"),p2,width = 18,height = 12)


#### ----->>> SELECT INTERESTED GENES FOR VISUALIZATION <<<-----
myeloid.pdac.tissue <- readRDS("input/mye.with.mast.pdac.tissue.sub.rds")
# macrophage
pdac.macro <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "Macrophage")))
pbmc.macro <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "Macrophage")))
pdac.tiss.macro <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "Macrophage")))

## -------->>> Up-regulated Genes <<<-----------

# GLUL
# 1) pdac.tumor_vs_paratumor
p.pdac.GLUL <- VlnPlot(pdac.macro,features = "GLUL",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.GLUL.pdf"),p.pdac.GLUL,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.GLUL <- VlnPlot(pdac.macro,features = "GLUL",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.tiss.GLUL.pdf"),p.pdac.tiss.GLUL,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.GLUL <- VlnPlot(pbmc.macro,features = "GLUL",pt.size = 0.0001,group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pbmc.GLUL.pdf"),p.pbmc.GLUL,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.GLUL.2 <- VlnPlot(pdac.tiss.macro,features = "GLUL",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.tiss.GLUL.2.pdf"),p.pdac.tiss.GLUL.2,width = 4,height = 3)

# SQSTM1
# 1) pdac.tumor_vs_paratumor
p.pdac.SQSTM1 <- VlnPlot(pdac.macro,features = "SQSTM1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.SQSTM1.pdf"),p.pdac.SQSTM1,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.SQSTM1 <- VlnPlot(pdac.macro,features = "SQSTM1",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.tiss.SQSTM1.pdf"),p.pdac.tiss.SQSTM1,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.SQSTM1 <- VlnPlot(pbmc.macro,features = "SQSTM1",pt.size = 0.0001,group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pbmc.SQSTM1.pdf"),p.pbmc.SQSTM1,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.SQSTM1.2 <- VlnPlot(pdac.tiss.macro,features = "SQSTM1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","macro.pdac.tiss.SQSTM1.2.pdf"),p.pdac.tiss.SQSTM1.2,width = 4,height = 3)
rm(pbmc.macro,pdac.macro,pdac.tiss.macro);gc()

# mast
pdac.mast <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "Mast cells")))
pbmc.mast <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "Mast cells")))
pdac.tiss.mast <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "Mast cells")))
# HSP90AB1
# 1) pdac.tumor_vs_paratumor
p.pdac.HSP90AB1 <- VlnPlot(pdac.mast,features = "HSP90AB1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.HSP90AB1.pdf"),p.pdac.HSP90AB1,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.HSP90AB1 <- VlnPlot(pdac.mast,features = "HSP90AB1",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.tiss.HSP90AB1.pdf"),p.pdac.tiss.HSP90AB1,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.HSP90AB1 <- VlnPlot(pbmc.mast,features = "HSP90AB1",group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pbmc.HSP90AB1.pdf"),p.pbmc.HSP90AB1,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.HSP90AB1.2 <- VlnPlot(pdac.tiss.mast,features = "HSP90AB1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.tiss.HSP90AB1.2.pdf"),p.pdac.tiss.HSP90AB1.2,width = 4,height = 3)

# HSP90AA1
# 1) pdac.tumor_vs_paratumor
p.pdac.HSP90AA1 <- VlnPlot(pdac.mast,features = "HSP90AA1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.HSP90AA1.pdf"),p.pdac.HSP90AA1,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.HSP90AA1 <- VlnPlot(pdac.mast,features = "HSP90AA1",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.tiss.HSP90AA1.pdf"),p.pdac.tiss.HSP90AA1,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.HSP90AA1 <- VlnPlot(pbmc.mast,features = "HSP90AA1",group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pbmc.HSP90AA1.pdf"),p.pbmc.HSP90AA1,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.HSP90AA1.2 <- VlnPlot(pdac.tiss.mast,features = "HSP90AA1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","mast.pdac.tiss.HSP90AA1.2.pdf"),p.pdac.tiss.HSP90AA1.2,width = 4,height = 3)

rm(pdac.mast,pbmc.mast,pdac.tiss.mast);gc()

# CD14 mono
pdac.CD14mono <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "CD14 Mono.")))
pbmc.CD14mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "CD14 Mono.")))
pdac.tiss.CD14mono <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "CD14 Mono.")))
# BIRC3
# 1) pdac.tumor_vs_paratumor
p.pdac.BIRC3 <- VlnPlot(pdac.CD14mono,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CD14mono.pdac.BIRC3.pdf"),p.pdac.BIRC3,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.BIRC3 <- VlnPlot(pdac.CD14mono,features = "BIRC3",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CD14mono.pdac.tiss.BIRC3.pdf"),p.pdac.tiss.BIRC3,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.BIRC3 <- VlnPlot(pbmc.CD14mono,features = "BIRC3",pt.size = 0.01,group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CD14mono.pbmc.BIRC3.pdf"),p.pbmc.BIRC3,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.BIRC3.2 <- VlnPlot(pdac.tiss.CD14mono,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CD14mono.pdac.tiss.BIRC3.2.pdf"),p.pdac.tiss.BIRC3.2,width = 4,height = 3)

rm(pdac.CD14mono,pbmc.CD14mono);gc()

# CD16 mono
pdac.CD16mono <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "CD16 Mono.")))
pbmc.CD16mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "CD16 Mono.")))
pdac.tiss.CD16mono <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "CD16 Mono.")))
# BIRC3
# 1) pdac.tumor_vs_paratumor
p.pdac.BIRC3 <- VlnPlot(pdac.CD16mono,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 5) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.BIRC3.pdf"),p.pdac.BIRC3,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.BIRC3 <- VlnPlot(pdac.CD16mono,features = "BIRC3",pt.size = 0.0001,group.by = "Tissue",y.max = 5) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.BIRC3.pdf"),p.pdac.tiss.BIRC3,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.BIRC3 <- VlnPlot(pbmc.CD16mono,features = "BIRC3",pt.size = 0.01,group.by = "Patient",y.max = 5) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pbmc.BIRC3.pdf"),p.pbmc.BIRC3,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.BIRC3.2 <- VlnPlot(pdac.tiss.CD16mono,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 5) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.BIRC3.2.pdf"),p.pdac.tiss.BIRC3.2,width = 4,height = 3)

# JAK3
# 1) pdac.tumor_vs_paratumor
p.pdac.JAK3 <- VlnPlot(pdac.CD16mono,features = "JAK3",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.JAK3.pdf"),p.pdac.JAK3,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.JAK3 <- VlnPlot(pdac.CD16mono,features = "JAK3",pt.size = 0.0001,group.by = "Tissue",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.JAK3.pdf"),p.pdac.tiss.JAK3,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.JAK3 <- VlnPlot(pbmc.CD16mono,features = "JAK3",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pbmc.JAK3.pdf"),p.pbmc.JAK3,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.JAK3.2 <- VlnPlot(pdac.tiss.CD16mono,features = "JAK3",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.JAK3.2.pdf"),p.pdac.tiss.JAK3.2,width = 4,height = 3)

# PPIA
# 1) pdac.tumor_vs_paratumor
p.pdac.PPIA <- VlnPlot(pdac.CD16mono,features = "PPIA",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.PPIA.pdf"),p.pdac.PPIA,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.PPIA <- VlnPlot(pdac.CD16mono,features = "PPIA",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.PPIA.pdf"),p.pdac.tiss.PPIA,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.PPIA <- VlnPlot(pbmc.CD16mono,features = "PPIA",pt.size = 0.01,group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pbmc.PPIA.pdf"),p.pbmc.PPIA,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.PPIA.2 <- VlnPlot(pdac.tiss.CD16mono,features = "PPIA",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.PPIA.2.pdf"),p.pdac.tiss.PPIA.2,width = 4,height = 3)

# TLR4
# 1) pdac.tumor_vs_paratumor
p.pdac.TLR4 <- VlnPlot(pdac.CD16mono,features = "TLR4",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.TLR4.pdf"),p.pdac.TLR4,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.TLR4 <- VlnPlot(pdac.CD16mono,features = "TLR4",pt.size = 0.0001,group.by = "Tissue",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.TLR4.pdf"),p.pdac.tiss.TLR4,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.TLR4 <- VlnPlot(pbmc.CD16mono,features = "TLR4",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pbmc.TLR4.pdf"),p.pbmc.TLR4,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.TLR4.2 <- VlnPlot(pdac.tiss.CD16mono,features = "TLR4",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","CD16mono.pdac.tiss.TLR4.2.pdf"),p.pdac.tiss.TLR4.2,width = 4,height = 3)
rm(pdac.CD16mono,pbmc.CD16mono,pdac.tiss.CD16mono);gc()

# DC
pdac.DC <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "DCs")))
pbmc.DC <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "DCs")))
pdac.tiss.DC <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "DCs")))

# CHMP1B
# 1) pdac.tumor_vs_paratumor
p.pdac.CHMP1B <- VlnPlot(pdac.DC,features = "CHMP1B",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","DC.pdac.CHMP1B.pdf"),p.pdac.CHMP1B,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.CHMP1B <- VlnPlot(pdac.DC,features = "CHMP1B",pt.size = 0.0001,group.by = "Tissue",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","DC.pdac.tiss.CHMP1B.pdf"),p.pdac.tiss.CHMP1B,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.CHMP1B <- VlnPlot(pbmc.DC,features = "CHMP1B",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","DC.pbmc.CHMP1B.pdf"),p.pbmc.CHMP1B,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.CHMP1B.2 <- VlnPlot(pdac.tiss.DC,features = "CHMP1B",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Upregulated/","DC.pdac.tiss.CHMP1B.2.pdf"),p.pdac.tiss.CHMP1B.2,width = 4,height = 3)


## -------->>> Down-regulated Genes <<<-----------
myeloid.pdac.tissue <- readRDS("input/mye.with.mast.pdac.tissue.sub.rds")
# macrophage
pdac.macro <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "Macrophage")))
pbmc.macro <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "Macrophage")))
pdac.tiss.macro <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "Macrophage")))
# SLC25A6
# 1) pdac.tumor_vs_paratumor
p.pdac.SLC25A6 <- VlnPlot(pdac.macro,features = "SLC25A6",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","macro.pdac.SLC25A6.pdf"),p.pdac.SLC25A6,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.SLC25A6 <- VlnPlot(pdac.macro,features = "SLC25A6",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","macro.pdac.tiss.SLC25A6.pdf"),p.pdac.tiss.SLC25A6,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.SLC25A6 <- VlnPlot(pbmc.macro,features = "SLC25A6",pt.size = 0.0001,group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","macro.pbmc.SLC25A6.pdf"),p.pbmc.SLC25A6,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.SLC25A6.2 <- VlnPlot(pdac.tiss.macro,features = "SLC25A6",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","macro.pdac.tiss.SLC25A6.2.pdf"),p.pdac.tiss.SLC25A6.2,width = 4,height = 3)

# mast
pdac.mast <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "Mast cells")))
pbmc.mast <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "Mast cells")))
pdac.tiss.mast <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "Mast cells")))
# BIRC3
# 1) pdac.tumor_vs_paratumor
p.pdac.BIRC3 <- VlnPlot(pdac.mast,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","mast.pdac.BIRC3.pdf"),p.pdac.BIRC3,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.BIRC3 <- VlnPlot(pdac.mast,features = "BIRC3",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","mast.pdac.tiss.BIRC3.pdf"),p.pdac.tiss.BIRC3,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.BIRC3 <- VlnPlot(pbmc.mast,features = "BIRC3",group.by = "Patient",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","mast.pbmc.BIRC3.pdf"),p.pbmc.BIRC3,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.BIRC3.2 <- VlnPlot(pdac.tiss.mast,features = "BIRC3",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","mast.pdac.tiss.BIRC3.2.pdf"),p.pdac.tiss.BIRC3.2,width = 4,height = 3)

# CD16 mono
pdac.CD16mono <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "CD16 Mono.")))
pbmc.CD16mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "CD16 Mono.")))
pdac.tiss.CD16mono <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "CD16 Mono.")))
# IFNGR1
# 1) pdac.tumor_vs_paratumor
p.pdac.IFNGR1 <- VlnPlot(pdac.CD16mono,features = "IFNGR1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","CD16mono.pdac.IFNGR1.pdf"),p.pdac.IFNGR1,width = 4,height = 3)
# 2) pdac.tissue_vs_pbmc
p.pdac.tiss.IFNGR1 <- VlnPlot(pdac.CD16mono,features = "IFNGR1",pt.size = 0.0001,group.by = "Tissue",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Tissue","PBMC")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","CD16mono.pdac.tiss.IFNGR1.pdf"),p.pdac.tiss.IFNGR1,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.IFNGR1.2 <- VlnPlot(pdac.CD16mono,features = "IFNGR1",pt.size = 0.0001,group.by = "State",y.max = 6) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","CD16mono.pdac.tiss.IFNGR1.2.pdf"),p.pdac.tiss.IFNGR1.2,width = 4,height = 3)

# DC
pdac.DC <- subset(myeloid.pdac,cells = rownames(subset(myeloid.pdac@meta.data, CellTypes == "DCs")))
pbmc.DC <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "DCs")))
pdac.tiss.DC <- subset(myeloid.pdac.tissue, cells = rownames(subset(myeloid.pdac.tissue@meta.data, CellTypes == "DCs")))
# PARP1
# 1) pdac.tumor_vs_paratumor
p.pdac.PARP1 <- VlnPlot(pdac.DC,features = "PARP1",group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","DC.pdac.PARP1.pdf"),p.pdac.PARP1,width = 4,height = 3)
# 3) pbmc.pdac_vs_normal
p.pbmc.PARP1 <- VlnPlot(pbmc.DC,features = "PARP1",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","DC.pbmc.PARP1.pdf"),p.pbmc.PARP1,width = 4,height = 3)
# 4) pdac.tissue.tumor_vs_paratumor
p.pdac.tiss.PARP1.2 <- VlnPlot(pdac.tiss.DC,features = "PARP1",pt.size = 0.0001,group.by = "State",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("Para-tumor","tumor")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","DC.pdac.tiss.PARP1.2.pdf"),p.pdac.tiss.PARP1.2,width = 4,height = 3)

# PBMC
pbmc.CD14mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes == "CD14 Mono.")))
# macro
p.pbmc.macro.EIF2AK2 <- VlnPlot(pbmc.macro,features = "EIF2AK2",group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","macro.pbmc.EIF2AK2.pdf"),p.pbmc.macro.EIF2AK2,width = 4,height = 3)
# mast
p.pbmc.mast.EIF2AK2 <- VlnPlot(pbmc.mast,features = "EIF2AK2",group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","mast.pbmc.EIF2AK2.pdf"),p.pbmc.mast.EIF2AK2,width = 4,height = 3)
# CD14 mono
p.pbmc.CD14mono.EIF2AK2 <- VlnPlot(pbmc.CD14mono,features = "EIF2AK2",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","CD14mono.pbmc.EIF2AK2.pdf"),p.pbmc.CD14mono.EIF2AK2,width = 4,height = 3)
# CD16 mono
p.pbmc.CD16mono.EIF2AK2 <- VlnPlot(pbmc.CD16mono,features = "EIF2AK2",pt.size = 0.01,group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","CD16mono.pbmc.EIF2AK2.pdf"),p.pbmc.CD16mono.EIF2AK2,width = 4,height = 3)
# DC
p.pbmc.DC.EIF2AK2 <- VlnPlot(pbmc.DC,features = "EIF2AK2",group.by = "Patient",y.max = 4) + NoLegend() +
  ggpubr::stat_compare_means(comparisons = list(c("PDAC","Normal")), label = "p.signif") +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 0,hjust = 0.4,vjust = 1))
ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/Downregulated/","DC.pbmc.EIF2AK2.pdf"),p.pbmc.DC.EIF2AK2,width = 4,height = 3)

rm(list = ls());gc()

# p.MMP9 <- VlnPlot(myeloid.pdac, features = "MMP9",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","MMP9.pdf"),p.MMP9,width = 4,height = 3)
# p.ACP5 <- VlnPlot(myeloid.pdac, features = "ACP5",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","ACP5.pdf"),p.ACP5,width = 4,height = 3)
# p.HSPB1 <- VlnPlot(myeloid.pdac, features = "HSPB1",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","HSPB1.pdf"),p.HSPB1,width = 4,height = 3)
# # mast cell
# p.CPA3 <- VlnPlot(myeloid.pdac, features = "CPA3",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CPA3.pdf"),p.CPA3,width = 4,height = 3)
# p.RUNX1 <- VlnPlot(myeloid.pdac, features = "RUNX1",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","RUNX1.pdf"),p.RUNX1,width = 4,height = 3)
# p.IGLL5 <- VlnPlot(myeloid.pdac, features = "IGLL5",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","IGLL5.pdf"),p.IGLL5,width = 4,height = 3)
# # CD14 Mono
# p.CRIP1 <- VlnPlot(myeloid.pdac, features = "CRIP1",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CRIP1.pdf"),p.CRIP1,width = 4,height = 3)
# p.TMEM176B <- VlnPlot(myeloid.pdac, features = "TMEM176B",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","TMEM176B.pdf"),p.TMEM176B,width = 4,height = 3)
# p.RAB11FIP1 <- VlnPlot(myeloid.pdac, features = "RAB11FIP1",group.by = "State",pt.size = 0) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","RAB11FIP1.pdf"),p.RAB11FIP1,width = 4,height = 3)
# # CD16 Mono
# p.CX3CR1 <- VlnPlot(myeloid.pdac, features = "CX3CR1",group.by = "State",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","CX3CR1.pdf"),p.CX3CR1,width = 4,height = 3)
# 
# p.VCAN <- VlnPlot(myeloid.pbmc, features = "VCAN",group.by = "Patient",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","VCAN.pdf"),p.VCAN,width = 4,height = 3)
# p.APOE <- VlnPlot(myeloid.pbmc, features = "APOE",group.by = "Patient",pt.size = 1e-5) + NoLegend() + 
#   theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0,hjust = 0.5))
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/genesOfInteret/","APOE.pdf"),p.APOE,width = 4,height = 3)
# rm(list = ls());gc()

## ----------- 6.Gene Enrichment Anlysis (GSEA/KEGG/GSVA) ---------------

if(FALSE){
  ## ---- GSEA ----
  ## PBMC (Tumor vs Normal): Macrophage
  myeloid.pbmc <- readRDS("input/mye.with.mast.pbmc.sub.rds")
  ## UP
  load("output/myeloid/DEgenes/mye.pbmc.macro.DEgenes.rda")
  up.genes <- rownames(mye.pbmc.macro.pdac_vs_normal[mye.pbmc.macro.pdac_vs_normal$avg_logFC > 0,])
  write(up.genes,"output/myeloid/DEgenes/mye.pbmc.marco.up.txt")
  macro.up.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[up.genes,]
  macro.up.expr <- data.frame(NAME = rownames(macro.up.expr), Description = rep("na",nrow(macro.up.expr)),macro.up.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.expr.gct",ncolumns = 1)
  write(c(nrow(macro.up.expr),(ncol(macro.up.expr)-2)), "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(macro.up.expr, "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(macro.up.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/macrophage/pdac_vs_normal/output/up/")
  
  # DOWN
  down.genes <- rownames(mye.pbmc.macro.pdac_vs_normal[mye.pbmc.macro.pdac_vs_normal$avg_logFC < 0,])
  write(down.genes,"output/myeloid/DEgenes/mye.pbmc.marco.down.txt")
  macro.down.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[down.genes,]
  macro.down.expr <- data.frame(NAME = rownames(macro.down.expr), Description = rep("na",nrow(macro.down.expr)),macro.down.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.expr.gct",ncolumns = 1)
  write(c(nrow(macro.down.expr),(ncol(macro.down.expr)-2)), "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(macro.down.expr, "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(macro.down.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/macrophage/pdac_vs_normal/input/down/macro.pdac_vs_normal.down.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/macrophage/pdac_vs_normal/input/up/macro.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/macrophage/pdac_vs_normal/output/up/")
  
  ## PBMC (Tumor vs Normal): Mast cells
  # UP
  load("output/myeloid/DEgenes/mye.pbmc.mast.DEgenes.rda")
  up.genes <- rownames(mye.pbmc.mast.pdac_vs_normal[mye.pbmc.mast.pdac_vs_normal$avg_logFC > 0,])
  write(up.genes,"output/myeloid/DEgenes/mye.pbmc.mast.up.txt")
  mast.up.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[up.genes,]
  mast.up.expr <- data.frame(NAME = rownames(mast.up.expr), Description = rep("na",nrow(mast.up.expr)),mast.up.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.expr.gct",ncolumns = 1)
  write(c(nrow(mast.up.expr),(ncol(mast.up.expr)-2)), "output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(mast.up.expr, "output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(mast.up.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/mast/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/macst/pdac_vs_normal/input/up/mast.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/mast/pdac_vs_normal/output/up/")
  
  # DOWN
  down.genes <- rownames(mye.pbmc.mast.pdac_vs_normal[mye.pbmc.mast.pdac_vs_normal$avg_logFC < 0,])
  write(down.genes,"output/myeloid/DEgenes/mye.pbmc.mast.down.txt")
  mast.down.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[down.genes,]
  mast.down.expr <- data.frame(NAME = rownames(mast.down.expr), Description = rep("na",nrow(mast.down.expr)),mast.down.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.expr.gct",ncolumns = 1)
  write(c(nrow(mast.down.expr),(ncol(mast.down.expr)-2)), "output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(mast.down.expr, "output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(mast.down.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/mast/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.expr.gct",
             input.cls = "output/myeloid/GSEA/macst/pdac_vs_normal/input/down/mast.pdac_vs_normal.down.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/mast/pdac_vs_normal/output/down/")
  
  ## PBMC (Tumor vs Normal): CD14 Mono.
  # UP
  load("output/myeloid/DEgenes/mye.pbmc.CD14.mono.DEgenes.rda")
  up.genes <- rownames(mye.pbmc.CD14.mono.pdac_vs_normal[mye.pbmc.CD14.mono.pdac_vs_normal$avg_logFC > 0,])
  write(up.genes,"output/myeloid/DEgenes/mye.pbmc.CD14.mono.up.txt")
  CD14.mono.up.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[up.genes,]
  CD14.mono.up.expr <- data.frame(NAME = rownames(CD14.mono.up.expr), Description = rep("na",nrow(CD14.mono.up.expr)),CD14.mono.up.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.expr.gct",ncolumns = 1)
  write(c(nrow(CD14.mono.up.expr),(ncol(CD14.mono.up.expr)-2)), "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(CD14.mono.up.expr, "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(CD14.mono.up.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/up/CD14.mono.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/output/up/")
  
  # DOWN
  down.genes <- rownames(mye.pbmc.CD14.mono.pdac_vs_normal[mye.pbmc.CD14.mono.pdac_vs_normal$avg_logFC < 0,])
  write(down.genes,"output/myeloid/DEgenes/mye.pbmc.CD14.mono.down.txt")
  CD14.mono.down.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[down.genes,]
  CD14.mono.down.expr <- data.frame(NAME = rownames(CD14.mono.down.expr), Description = rep("na",nrow(CD14.mono.down.expr)),CD14.mono.down.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.expr.gct",ncolumns = 1)
  write(c(nrow(CD14.mono.down.expr),(ncol(CD14.mono.down.expr)-2)), "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(CD14.mono.down.expr, "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(CD14.mono.down.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.expr.gct",
             input.cls = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/input/down/CD14.mono.pdac_vs_normal.down.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/CD14_Mono/pdac_vs_normal/output/down/")
  
  ## PBMC (Tumor vs Normal): CD16 Mono.
  # UP
  load("output/myeloid/DEgenes/mye.pbmc.CD16.mono.DEgenes.rda")
  up.genes <- rownames(mye.pbmc.CD16.mono.pdac_vs_normal[mye.pbmc.CD16.mono.pdac_vs_normal$avg_logFC > 0,])
  write(up.genes,"output/myeloid/DEgenes/mye.pbmc.CD16.mono.up.txt")
  CD16.mono.up.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[up.genes,]
  CD16.mono.up.expr <- data.frame(NAME = rownames(CD16.mono.up.expr), Description = rep("na",nrow(CD16.mono.up.expr)),CD16.mono.up.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.expr.gct",ncolumns = 1)
  write(c(nrow(CD16.mono.up.expr),(ncol(CD16.mono.up.expr)-2)), "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(CD16.mono.up.expr, "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(CD16.mono.up.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/up/CD16.mono.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/output/up/")
  
  # DOWN
  down.genes <- rownames(mye.pbmc.CD16.mono.pdac_vs_normal[mye.pbmc.CD16.mono.pdac_vs_normal$avg_logFC < 0,])
  write(down.genes,"output/myeloid/DEgenes/mye.pbmc.CD16.mono.down.txt")
  CD16.mono.down.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[down.genes,]
  CD16.mono.down.expr <- data.frame(NAME = rownames(CD16.mono.down.expr), Description = rep("na",nrow(CD16.mono.down.expr)),CD16.mono.down.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.expr.gct",ncolumns = 1)
  write(c(nrow(CD16.mono.down.expr),(ncol(CD16.mono.down.expr)-2)), "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(CD16.mono.down.expr, "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(CD16.mono.down.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.expr.gct",
             input.cls = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/input/down/CD16.mono.pdac_vs_normal.down.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/CD16_Mono/pdac_vs_normal/output/down/")
  
  ## PBMC (Tumor vs Normal): DCs
  # UP
  load("output/myeloid/DEgenes/mye.pbmc.DC.DEgenes.rda")
  up.genes <- rownames(mye.pbmc.DC.pdac_vs_normal[mye.pbmc.DC.pdac_vs_normal$avg_logFC > 0,])
  write(up.genes,"output/myeloid/DEgenes/mye.pbmc.DC.up.txt")
  DC.up.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[up.genes,]
  DC.up.expr <- data.frame(NAME = rownames(DC.up.expr), Description = rep("na",nrow(DC.up.expr)),DC.up.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.expr.gct",ncolumns = 1)
  write(c(nrow(DC.up.expr),(ncol(DC.up.expr)-2)), "output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(DC.up.expr, "output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(DC.up.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.expr.gct",
             input.cls = "output/myeloid/GSEA/DCs/pdac_vs_normal/input/up/DC.pdac_vs_normal.up.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/DCs/pdac_vs_normal/output/up/")
  
  # DOWN
  down.genes <- rownames(mye.pbmc.DC.pdac_vs_normal[mye.pbmc.DC.pdac_vs_normal$avg_logFC < 0,])
  write(down.genes,"output/myeloid/DEgenes/mye.pbmc.DC.down.txt")
  DC.down.expr <- GetAssayData(myeloid.pbmc, assay = "RNA", slot = "counts")[down.genes,]
  DC.down.expr <- data.frame(NAME = rownames(DC.down.expr), Description = rep("na",nrow(DC.down.expr)),DC.down.expr,stringsAsFactors = F)
  write("#1.2", "output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.expr.gct",ncolumns = 1)
  write(c(nrow(DC.down.expr),(ncol(DC.down.expr)-2)), "output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(DC.down.expr, "output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(DC.down.expr)-2),2,1)
  tmp <- table(myeloid.pbmc@meta.data$Patient)
  line.2 <- c("#",names(tmp))
  line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.group.cls",ncolumns = length(line.1),sep = "\t")
  write(line.2,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.group.cls",ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3,"output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.group.cls",ncolumns = length(line.3),append = T,sep = "\t")
  
  GSEA::GSEA(input.ds = "output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.expr.gct",
             input.cls = "output/myeloid/GSEA/DCs/pdac_vs_normal/input/down/DC.pdac_vs_normal.down.group.cls",
             gs.db = "input/h.all.v7.4.symbols.gmt",
             output.directory = "output/myeloid/GSEA/DCs/pdac_vs_normal/output/down/")
}

## ---- KEGG ----
pacman::p_load(clusterProfiler,org.Hs.eg.db,dplyr,patchwork,enrichplot,ggplot2)

EnrichmentAnalysis <- function(DEgenes, enrich.method = c("GO","KEGG")) {
  enrich.method <- match.arg(enrich.method)
  if(enrich.method == "GO") {
    ego_ALL <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL",
                        pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_CC <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_MF <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_BP <- enrichGO(gene = DEgenes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_CC@result$Description <- substring(ego_CC@result$Description, 1, 70)
    ego_BP@result$Description <- substring(ego_BP@result$Description, 1, 70)
    ego_MF@result$Description <- substring(ego_MF@result$Description, 1, 70)
    p_CC <- barplot(ego_BP, showCategory = 10, title = "Cellular Component")
    p_BP <- barplot(ego_CC, showCategory = 10, title = "Biological Process")
    p_MF <- barplot(ego_MF, showCategory = 10, title = "Molecular Function")
    pc <- p_BP | p_CC | p_MF
    ggsave("GO.pdf",plot = pc, device = "pdf",width = 12, height = 10)
    
    ego_results <- list(ego_ALL = ego_ALL, ego_BP = ego_BP, ego_CC = ego_CC, ego_MF = ego_MF)
    return(ego_results)
  } else if(enrich.method == "KEGG") {
    genelist <- bitr(DEgenes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    genelist <- pull(genelist, ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = "hsa")
    # p <- dotplot(ekegg, showCategory = 20)
    # ggsave("KEGG.pdf",plot = p, device = "pdf", width = 12, height = 10)
    return(ekegg)
  }
}

#### ---------->>>>>>>>>>>>> PBMC (PDAC vs Normal) <<<<<<<<<<<<<<<------------
# macrophage
load("output/myeloid/DEgenes/mye.pbmc.macro.DEgenes.rda")
up.genes <- mye.pbmc.macro.pdac_vs_normal[mye.pbmc.macro.pdac_vs_normal$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "macro.pdac.up",row.names = F)
down.genes <- mye.pbmc.macro.pdac_vs_normal[mye.pbmc.macro.pdac_vs_normal$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "macro.pdac.down",row.names = F, append = T)
# mast cell
load("output/myeloid/DEgenes/mye.pbmc.mast.DEgenes.rda")
up.genes <- mye.pbmc.mast.pdac_vs_normal[mye.pbmc.mast.pdac_vs_normal$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "mast.pdac.up",row.names = F, append = T)
down.genes <- mye.pbmc.mast.pdac_vs_normal[mye.pbmc.mast.pdac_vs_normal$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "mast.pdac.down",row.names = F, append = T)
# CD14 mono
load("output/myeloid/DEgenes/mye.pbmc.CD14.mono.DEgenes.rda")
up.genes <- mye.pbmc.CD14.mono.pdac_vs_normal[mye.pbmc.CD14.mono.pdac_vs_normal$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "CD14.mono.pdac.up",row.names = F, append = T)
down.genes <- mye.pbmc.CD14.mono.pdac_vs_normal[mye.pbmc.CD14.mono.pdac_vs_normal$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "CD14.mono.pdac.down",row.names = F, append = T)
# CD16 mono
load("output/myeloid/DEgenes/mye.pbmc.CD16.mono.DEgenes.rda")
up.genes <- mye.pbmc.CD16.mono.pdac_vs_normal[mye.pbmc.CD16.mono.pdac_vs_normal$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "CD16.mono.pdac.up",row.names = F, append = T)
down.genes <- mye.pbmc.CD16.mono.pdac_vs_normal[mye.pbmc.CD16.mono.pdac_vs_normal$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "CD16.mono.pdac.down",row.names = F, append = T)
# DC
load("output/myeloid/DEgenes/mye.pbmc.DC.DEgenes.rda")
up.genes <- mye.pbmc.DC.pdac_vs_normal[mye.pbmc.DC.pdac_vs_normal$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "DC.pdac.up",row.names = F, append = T)
down.genes <- mye.pbmc.DC.pdac_vs_normal[mye.pbmc.DC.pdac_vs_normal$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pbmc.tumor_vs_normal.diff.pathway.xlsx",sheetName = "DC.pdac.down",row.names = F, append = T)


#### ---->>> PDAC (tumor vs. para-tumor && tumor.tissue vs. paratumor.tissue && Tissue vs. PBMC) <<<----

# --->>> 1.Tumor vs. Para-Tumor <<<---
# macrophage
load("output/myeloid/DEgenes/mye.pdac.marco.DEgenes.rda")
up.genes <- mye.pdac.macro.tumor_vs_paratumor[mye.pdac.macro.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "macro.pdac.up",row.names = F)
down.genes <- mye.pdac.macro.tumor_vs_paratumor[mye.pdac.macro.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "macro.pdac.down",row.names = F, append = T)
# mast
load("output/myeloid/DEgenes/mye.pdac.mast.DEgenes.rda")
up.genes <- mye.pdac.mast.tumor_vs_paratumor[mye.pdac.mast.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "mast.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.mast.tumor_vs_paratumor[mye.pdac.mast.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "mast.pdac.down",row.names = F, append = T)
# CD14 mono
load("output/myeloid/DEgenes/mye.pdac.CD14.Mono.DEgenes.rda")
up.genes <- mye.pdac.CD14.mono.tumor_vs_paratumor[mye.pdac.CD14.mono.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD14mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.CD14.mono.tumor_vs_paratumor[mye.pdac.CD14.mono.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD14mono.pdac.down",row.names = F, append = T)
# CD16 mono
load("output/myeloid/DEgenes/mye.pdac.CD16.Mono.DEgenes.rda")
up.genes <- mye.pdac.CD16.mono.tumor_vs_paratumor[mye.pdac.CD16.mono.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD16mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.CD16.mono.tumor_vs_paratumor[mye.pdac.CD16.mono.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD16mono.pdac.down",row.names = F, append = T)
# DC
load("output/myeloid/DEgenes/mye.pdac.DC.DEgenes.rda")
up.genes <- mye.pdac.DC.tumor_vs_paratumor[mye.pdac.DC.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "DC.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.DC.tumor_vs_paratumor[mye.pdac.DC.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "DC.pdac.down",row.names = F, append = T)

# --->>> 2.tumor.tissue vs. paratumor.tissue <<<---
# macrophage
up.genes <- mye.pdac.tissue.macro.tumor_vs_paratumor[mye.pdac.tissue.macro.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "macro.pdac.up",row.names = F)
down.genes <- mye.pdac.tissue.macro.tumor_vs_paratumor[mye.pdac.tissue.macro.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "macro.pdac.down",row.names = F, append = T)
# mast cells
up.genes <- mye.pdac.tissue.mast.tumor_vs_paratumor[mye.pdac.tissue.mast.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "mast.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.tissue.mast.tumor_vs_paratumor[mye.pdac.tissue.mast.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "mast.pdac.down",row.names = F, append = T)
# CD14 mono
up.genes <- mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD14mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.tissue.CD14.mono.tumor_vs_paratumor[mye.pdac.tissue.CD14.mono.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD14mono.pdac.down",row.names = F, append = T)
# CD16 mono
up.genes <- mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD16mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.tissue.CD16.mono.tumor_vs_paratumor[mye.pdac.tissue.CD16.mono.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "CD16mono.pdac.down",row.names = F, append = T)
# DC
up.genes <- mye.pdac.tissue.DC.tumor_vs_paratumor[mye.pdac.tissue.DC.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "DC.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.tissue.DC.tumor_vs_paratumor[mye.pdac.tissue.DC.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "mye.pdac.tissue.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "DC.pdac.down",row.names = F, append = T)

# --->>> 3.Tissue vs. PBMC <<<---
# macrophage
up.genes <- mye.pdac.macro.tissue_vs_pbmc[mye.pdac.macro.tissue_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "macro.pdac.up",row.names = F)
down.genes <- mye.pdac.macro.tissue_vs_pbmc[mye.pdac.macro.tissue_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "macro.pdac.down",row.names = F, append = T)
# mast
up.genes <- mye.pdac.mast.tissue_vs_pbmc[mye.pdac.mast.tissue_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "mast.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.mast.tissue_vs_pbmc[mye.pdac.mast.tissue_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "mast.pdac.down",row.names = F, append = T)
# CD14 mono
up.genes <- mye.pdac.CD14.mono.tissue_vs_pbmc[mye.pdac.CD14.mono.tissue_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "CD14mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.CD14.mono.tissue_vs_pbmc[mye.pdac.CD14.mono.tissue_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "CD14mono.pdac.down",row.names = F, append = T)
# CD16 mono
up.genes <- mye.pdac.CD16.mono.tissue_vs_pbmc[mye.pdac.CD16.mono.tissue_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "CD16mono.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.CD16.mono.tissue_vs_pbmc[mye.pdac.CD16.mono.tissue_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "CD16mono.pdac.down",row.names = F, append = T)
# DC
up.genes <- mye.pdac.DC.tissue_vs_pbmc[mye.pdac.DC.tissue_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "DC.pdac.up",row.names = F,append = T)
down.genes <- mye.pdac.DC.tissue_vs_pbmc[mye.pdac.DC.tissue_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/myeloid/pathway/mye.pdac.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "DC.pdac.down",row.names = F, append = T)

rm(list = ls());gc()


# --->>> visualize selected pathway <<<---

# --->> 1.myeloid.pbmc <<---
## macro
sel.macro <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pbmc.xlsx",sheetIndex = 1)
sel.macro$value <- -log10(sel.macro$pvalue)
sel.macro$value[sel.macro$Group=="macro_normal"] <- -sel.macro$value[sel.macro$Group=="macro_normal"]
sel.macro$hjust <- ifelse(sel.macro$Group=="macro_pdac",1,0)
p.macro <- ggplot(data = sel.macro, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("macro_pdac"="#CA73AA","macro_normal"="#F7CBBF")) +
  scale_fill_manual(values = c("macro_pdac"="#CA73AA","macro_normal"="#F7CBBF")) +
  geom_text(data = sel.macro, mapping = aes(label=Description, y = 0), hjust = sel.macro$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.macro$value))),ceiling(max(abs(sel.macro$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pbmc.macro.pdf"),p.macro,width = 6,height = 6)

## mast
sel.mast <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pbmc.xlsx",sheetIndex = 2)
sel.mast$value <- -log10(sel.mast$pvalue)
sel.mast$value[sel.mast$Group=="mast_normal"] <- -sel.mast$value[sel.mast$Group=="mast_normal"]
sel.mast$hjust <- ifelse(sel.mast$Group=="mast_pdac",1,0)
p.mast <- ggplot(data = sel.mast, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("mast_pdac"="#CA73AA","mast_normal"="#F7CBBF")) +
  scale_fill_manual(values = c("mast_pdac"="#CA73AA","mast_normal"="#F7CBBF")) +
  geom_text(data = sel.mast, mapping = aes(label=Description, y = 0), hjust = sel.mast$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.mast$value))),ceiling(max(abs(sel.mast$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pbmc.mast.pdf"),p.mast,width = 6,height = 7)

## CD14mono
sel.CD14mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pbmc.xlsx",sheetIndex = 3)
sel.CD14mono$value <- -log10(sel.CD14mono$pvalue)
sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_normal"] <- -sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_normal"]
sel.CD14mono$hjust <- ifelse(sel.CD14mono$Group=="CD14mono_pdac",1,0)
p.CD14mono <- ggplot(data = sel.CD14mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD14mono_pdac"="#CA73AA","CD14mono_normal"="#F7CBBF")) +
  scale_fill_manual(values = c("CD14mono_pdac"="#CA73AA","CD14mono_normal"="#F7CBBF")) +
  geom_text(data = sel.CD14mono, mapping = aes(label=Description, y = 0), hjust = sel.CD14mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD14mono$value))),ceiling(max(abs(sel.CD14mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pbmc.CD14mono.pdf"),p.CD14mono,width = 6,height = 6)

## CD16mono
sel.CD16mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pbmc.xlsx",sheetIndex = 4)
sel.CD16mono$value <- -log10(sel.CD16mono$pvalue)
sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_normal"] <- -sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_normal"]
sel.CD16mono$hjust <- ifelse(sel.CD16mono$Group=="CD16mono_pdac",1,0)
p.CD16mono <- ggplot(data = sel.CD16mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD16mono_pdac"="#CA73AA","CD16mono_normal"="#F7CBBF")) +
  scale_fill_manual(values = c("CD16mono_pdac"="#CA73AA","CD16mono_normal"="#F7CBBF")) +
  geom_text(data = sel.CD16mono, mapping = aes(label=Description, y = 0), hjust = sel.CD16mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD16mono$value))),ceiling(max(abs(sel.CD16mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pbmc.CD16mono.pdf"),p.CD16mono,width = 6,height = 6)

## DC
sel.DC <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pbmc.xlsx",sheetIndex = 5)
sel.DC$value <- -log10(sel.DC$pvalue)
sel.DC$value[sel.DC$Group=="DC_normal"] <- -sel.DC$value[sel.DC$Group=="DC_normal"]
sel.DC$hjust <- ifelse(sel.DC$Group=="DC_pdac",1,0)
p.DC <- ggplot(data = sel.DC, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("DC_pdac"="#CA73AA","DC_normal"="#F7CBBF")) +
  scale_fill_manual(values = c("DC_pdac"="#CA73AA","DC_normal"="#F7CBBF")) +
  geom_text(data = sel.DC, mapping = aes(label=Description, y = 0), hjust = sel.DC$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.DC$value))),ceiling(max(abs(sel.DC$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pbmc.DC.pdf"),p.DC,width = 6,height = 6)

# --->> 2.myeloid.pdac.tumor_vs_paratumor <<---
## macro
sel.macro <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tumor_vs_paratumor.xlsx",sheetIndex = 1)
sel.macro$value <- -log10(sel.macro$pvalue)
sel.macro$value[sel.macro$Group=="macro_paratumor"] <- -sel.macro$value[sel.macro$Group=="macro_paratumor"]
sel.macro$hjust <- ifelse(sel.macro$Group=="macro_tumor",1,0)
p.macro <- ggplot(data = sel.macro, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("macro_tumor"="#CA73AA","macro_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("macro_tumor"="#CA73AA","macro_paratumor"="#F7CBBF")) +
  geom_text(data = sel.macro, mapping = aes(label=Description, y = 0), hjust = sel.macro$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.macro$value))),ceiling(max(abs(sel.macro$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.macro.tumor_vs_paratumor.pdf"),p.macro,width = 6,height = 6)

## mast
sel.mast <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tumor_vs_paratumor.xlsx",sheetIndex = 2)
sel.mast$value <- -log10(sel.mast$pvalue)
sel.mast$value[sel.mast$Group=="mast_paratumor"] <- -sel.mast$value[sel.mast$Group=="mast_paratumor"]
sel.mast$hjust <- ifelse(sel.mast$Group=="mast_tumor",1,0)
p.mast <- ggplot(data = sel.mast, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("mast_tumor"="#CA73AA","mast_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("mast_tumor"="#CA73AA","mast_paratumor"="#F7CBBF")) +
  geom_text(data = sel.mast, mapping = aes(label=Description, y = 0), hjust = sel.mast$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.mast$value))),ceiling(max(abs(sel.mast$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.mast.tumor_vs_paratumor.pdf"),p.mast,width = 6,height = 7)

## CD14mono
sel.CD14mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tumor_vs_paratumor.xlsx",sheetIndex = 3)
sel.CD14mono$value <- -log10(sel.CD14mono$pvalue)
sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_paratumor"] <- -sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_paratumor"]
sel.CD14mono$hjust <- ifelse(sel.CD14mono$Group=="CD14mono_tumor",1,0)
p.CD14mono <- ggplot(data = sel.CD14mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD14mono_tumor"="#CA73AA","CD14mono_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("CD14mono_tumor"="#CA73AA","CD14mono_paratumor"="#F7CBBF")) +
  geom_text(data = sel.CD14mono, mapping = aes(label=Description, y = 0), hjust = sel.CD14mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD14mono$value))),ceiling(max(abs(sel.CD14mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.CD14mono.tumor_vs_paratumor.pdf"),p.CD14mono,width = 6,height = 6)

## CD16mono
sel.CD16mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tumor_vs_paratumor.xlsx",sheetIndex = 4)
sel.CD16mono$value <- -log10(sel.CD16mono$pvalue)
sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_paratumor"] <- -sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_paratumor"]
sel.CD16mono$hjust <- ifelse(sel.CD16mono$Group=="CD16mono_tumor",1,0)
p.CD16mono <- ggplot(data = sel.CD16mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD16mono_tumor"="#CA73AA","CD16mono_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("CD16mono_tumor"="#CA73AA","CD16mono_paratumor"="#F7CBBF")) +
  geom_text(data = sel.CD16mono, mapping = aes(label=Description, y = 0), hjust = sel.CD16mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD16mono$value))),ceiling(max(abs(sel.CD16mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.CD16mono.tumor_vs_paratumor.pdf"),p.CD16mono,width = 6,height = 6)

## DC
sel.DC <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tumor_vs_paratumor.xlsx",sheetIndex = 5)
sel.DC$value <- -log10(sel.DC$pvalue)
sel.DC$value[sel.DC$Group=="DC_paratumor"] <- -sel.DC$value[sel.DC$Group=="DC_paratumor"]
sel.DC$hjust <- ifelse(sel.DC$Group=="DC_tumor",1,0)
p.DC <- ggplot(data = sel.DC, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("DC_tumor"="#CA73AA","DC_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("DC_tumor"="#CA73AA","DC_paratumor"="#F7CBBF")) +
  geom_text(data = sel.DC, mapping = aes(label=Description, y = 0), hjust = sel.DC$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.DC$value))),ceiling(max(abs(sel.DC$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.DC.tumor_vs_paratumor.pdf"),p.DC,width = 6,height = 6)


# --->> 3.myeloid.pdac.tissue.tumor_vs_paratumor <<---
## macro
sel.macro <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue.tumor_vs_paratumor.xlsx",sheetIndex = 1)
sel.macro$value <- -log10(sel.macro$pvalue)
sel.macro$value[sel.macro$Group=="macro_paratumor"] <- -sel.macro$value[sel.macro$Group=="macro_paratumor"]
sel.macro$hjust <- ifelse(sel.macro$Group=="macro_tumor",1,0)
p.macro <- ggplot(data = sel.macro, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("macro_tumor"="#CA73AA","macro_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("macro_tumor"="#CA73AA","macro_paratumor"="#F7CBBF")) +
  geom_text(data = sel.macro, mapping = aes(label=Description, y = 0), hjust = sel.macro$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.macro$value))),ceiling(max(abs(sel.macro$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.tiss.macro.tumor_vs_paratumor.pdf"),p.macro,width = 6,height = 6)

## mast
sel.mast <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue.tumor_vs_paratumor.xlsx",sheetIndex = 2)
sel.mast$value <- -log10(sel.mast$pvalue)
sel.mast$value[sel.mast$Group=="mast_paratumor"] <- -sel.mast$value[sel.mast$Group=="mast_paratumor"]
sel.mast$hjust <- ifelse(sel.mast$Group=="mast_tumor",1,0)
p.mast <- ggplot(data = sel.mast, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("mast_tumor"="#CA73AA","mast_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("mast_tumor"="#CA73AA","mast_paratumor"="#F7CBBF")) +
  geom_text(data = sel.mast, mapping = aes(label=Description, y = 0), hjust = sel.mast$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.mast$value))),ceiling(max(abs(sel.mast$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.tiss.mast.tumor_vs_paratumor.pdf"),p.mast,width = 6,height = 7)

## CD14mono
sel.CD14mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue.tumor_vs_paratumor.xlsx",sheetIndex = 3)
sel.CD14mono$value <- -log10(sel.CD14mono$pvalue)
sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_paratumor"] <- -sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_paratumor"]
sel.CD14mono$hjust <- ifelse(sel.CD14mono$Group=="CD14mono_tumor",1,0)
p.CD14mono <- ggplot(data = sel.CD14mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD14mono_tumor"="#CA73AA","CD14mono_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("CD14mono_tumor"="#CA73AA","CD14mono_paratumor"="#F7CBBF")) +
  geom_text(data = sel.CD14mono, mapping = aes(label=Description, y = 0), hjust = sel.CD14mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD14mono$value))),ceiling(max(abs(sel.CD14mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.tiss.CD14mono.tumor_vs_paratumor.pdf"),p.CD14mono,width = 6,height = 6)

## CD16mono
sel.CD16mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue.tumor_vs_paratumor.xlsx",sheetIndex = 4)
sel.CD16mono$value <- -log10(sel.CD16mono$pvalue)
sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_paratumor"] <- -sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_paratumor"]
sel.CD16mono$hjust <- ifelse(sel.CD16mono$Group=="CD16mono_tumor",1,0)
p.CD16mono <- ggplot(data = sel.CD16mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD16mono_tumor"="#CA73AA","CD16mono_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("CD16mono_tumor"="#CA73AA","CD16mono_paratumor"="#F7CBBF")) +
  geom_text(data = sel.CD16mono, mapping = aes(label=Description, y = 0), hjust = sel.CD16mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD16mono$value))),ceiling(max(abs(sel.CD16mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.tiss.CD16mono.tumor_vs_paratumor.pdf"),p.CD16mono,width = 6,height = 6)

## DC
sel.DC <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue.tumor_vs_paratumor.xlsx",sheetIndex = 5)
sel.DC$value <- -log10(sel.DC$pvalue)
sel.DC$value[sel.DC$Group=="DC_paratumor"] <- -sel.DC$value[sel.DC$Group=="DC_paratumor"]
sel.DC$hjust <- ifelse(sel.DC$Group=="DC_tumor",1,0)
p.DC <- ggplot(data = sel.DC, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("DC_tumor"="#CA73AA","DC_paratumor"="#F7CBBF")) +
  scale_fill_manual(values = c("DC_tumor"="#CA73AA","DC_paratumor"="#F7CBBF")) +
  geom_text(data = sel.DC, mapping = aes(label=Description, y = 0), hjust = sel.DC$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.DC$value))),ceiling(max(abs(sel.DC$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.tiss.DC.tumor_vs_paratumor.pdf"),p.DC,width = 6,height = 6)


# --->> 4.myeloid.pdac.tissue_vs_pbmc <<---
## macro
sel.macro <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue_vs_pbmc.xlsx",sheetIndex = 1)
sel.macro$value <- -log10(sel.macro$pvalue)
sel.macro$value[sel.macro$Group=="macro_pbmc"] <- -sel.macro$value[sel.macro$Group=="macro_pbmc"]
sel.macro$hjust <- ifelse(sel.macro$Group=="macro_tissue",1,0)
p.macro <- ggplot(data = sel.macro, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("macro_tissue"="#CA73AA","macro_pbmc"="#F7CBBF")) +
  scale_fill_manual(values = c("macro_tissue"="#CA73AA","macro_pbmc"="#F7CBBF")) +
  geom_text(data = sel.macro, mapping = aes(label=Description, y = 0), hjust = sel.macro$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.macro$value))),ceiling(max(abs(sel.macro$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.macro.tissue_vs_pbmc.pdf"),p.macro,width = 6,height = 6)

## mast
sel.mast <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue_vs_pbmc.xlsx",sheetIndex = 2)
sel.mast$value <- -log10(sel.mast$pvalue)
sel.mast$value[sel.mast$Group=="mast_pbmc"] <- -sel.mast$value[sel.mast$Group=="mast_pbmc"]
sel.mast$hjust <- ifelse(sel.mast$Group=="mast_tissue",1,0)
p.mast <- ggplot(data = sel.mast, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("mast_tissue"="#CA73AA","mast_pbmc"="#F7CBBF")) +
  scale_fill_manual(values = c("mast_tissue"="#CA73AA","mast_pbmc"="#F7CBBF")) +
  geom_text(data = sel.mast, mapping = aes(label=Description, y = 0), hjust = sel.mast$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.mast$value))),ceiling(max(abs(sel.mast$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.mast.tissue_vs_pbmc.pdf"),p.mast,width = 6,height = 7)

## CD14mono
sel.CD14mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue_vs_pbmc.xlsx",sheetIndex = 3)
sel.CD14mono$value <- -log10(sel.CD14mono$pvalue)
sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_pbmc"] <- -sel.CD14mono$value[sel.CD14mono$Group=="CD14mono_pbmc"]
sel.CD14mono$hjust <- ifelse(sel.CD14mono$Group=="CD14mono_tissue",1,0)
p.CD14mono <- ggplot(data = sel.CD14mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD14mono_tissue"="#CA73AA","CD14mono_pbmc"="#F7CBBF")) +
  scale_fill_manual(values = c("CD14mono_tissue"="#CA73AA","CD14mono_pbmc"="#F7CBBF")) +
  geom_text(data = sel.CD14mono, mapping = aes(label=Description, y = 0), hjust = sel.CD14mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD14mono$value))),ceiling(max(abs(sel.CD14mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.CD14mono.tissue_vs_pbmc.pdf"),p.CD14mono,width = 6,height = 6)

## CD16mono
sel.CD16mono <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue_vs_pbmc.xlsx",sheetIndex = 4)
sel.CD16mono$value <- -log10(sel.CD16mono$pvalue)
sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_pbmc"] <- -sel.CD16mono$value[sel.CD16mono$Group=="CD16mono_pbmc"]
sel.CD16mono$hjust <- ifelse(sel.CD16mono$Group=="CD16mono_tissue",1,0)
p.CD16mono <- ggplot(data = sel.CD16mono, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("CD16mono_tissue"="#CA73AA","CD16mono_pbmc"="#F7CBBF")) +
  scale_fill_manual(values = c("CD16mono_tissue"="#CA73AA","CD16mono_pbmc"="#F7CBBF")) +
  geom_text(data = sel.CD16mono, mapping = aes(label=Description, y = 0), hjust = sel.CD16mono$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.CD16mono$value))),ceiling(max(abs(sel.CD16mono$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.CD16mono.tissue_vs_pbmc.pdf"),p.CD16mono,width = 6,height = 6)

## DC
sel.DC <- xlsx::read.xlsx("output/myeloid/pathway/forVis.mye.pdac.tissue_vs_pbmc.xlsx",sheetIndex = 5)
sel.DC$value <- -log10(sel.DC$pvalue)
sel.DC$value[sel.DC$Group=="DC_pbmc"] <- -sel.DC$value[sel.DC$Group=="DC_pbmc"]
sel.DC$hjust <- ifelse(sel.DC$Group=="DC_tissue",1,0)
p.DC <- ggplot(data = sel.DC, aes(factor(Description,levels = Description), value, color = Group, fill = Group)) +
  geom_bar(stat = "identity",width = 0.5) + theme_bw() + coord_flip() +
  # geom_hline(yintercept = 0, color = 'grey50', linetype = 'solid') +
  scale_color_manual(values = c("DC_tissue"="#CA73AA","DC_pbmc"="#F7CBBF")) +
  scale_fill_manual(values = c("DC_tissue"="#CA73AA","DC_pbmc"="#F7CBBF")) +
  geom_text(data = sel.DC, mapping = aes(label=Description, y = 0), hjust = sel.DC$hjust, size = 5,color="black") +
  labs(y = "-Log10(pvalue)", color = NULL, fill = NULL) +
  scale_y_continuous(limits = c(-ceiling(max(abs(sel.DC$value))),ceiling(max(abs(sel.DC$value))))) +
  theme(
    axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.text.x = element_text(color = "black", size = 18),
    # legend.title = element_text(color = "black",size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.position = "top",
    panel.grid = element_blank(),panel.border = element_blank(),
    axis.line.x = element_line(size = 1.5, colour = "black"),
    axis.ticks.x = element_line(size = 1.5, colour = "black")
  )
ggsave(paste0("figures/celltypes/myeloid/pathway/","mye.pdac.DC.tissue_vs_pbmc.pdf"),p.DC,width = 6,height = 6)


if(F){
  ## ---- GSVA ----
  library(Seurat)
  library(GSVA)
  library(tidyverse)
  
  myeloid.pbmc <- readRDS("input/mye.with.mast.pbmc.sub.rds")
  
  gmt2list <- function(gmt.file){
    sets <- as.list(read_lines(gmt.file))
    for(i in 1:length(sets)){
      tmp <- str_split(sets[[i]],'\t')
      n <- length(tmp[[1]])
      names(sets)[i] <- tmp[[1]][1]
      sets[[i]] <- tmp[[1]][3:n]
      rm(tmp,n)
    }
    return(sets)
  }
  
  s.sets <- gmt2list(gmt.file = "input/msigdb.v7.4.symbols.gmt")
  
  # macro
  macro <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes=="Macrophage")))
  macro.expr <- as.matrix(macro@assays$RNA@counts)
  macro.expr <- macro.expr[rowSums(macro.expr) >= 1,]
  macro.meta <- macro@meta.data[,c("Patient")];names(macro.meta) <- colnames(macro)
  es.mat <- GSVA::gsva(macro.expr,s.sets,kcdf="Poisson")
  write.table(es.mat, file = paste0("output/macro.gsva.txt",row.names = T, col.names = NA, sep = "\t"))
  pheatmap::pheatmap(es.mat, show_rownames = 1, show_colnames = 0, annotation_col = macro.meta,
                     fontsize_row = 5, filename = "output/macro.gsva.pdf",width=15, height=12)
  # mast
  mast <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes=="Mast cells")))
  mast.expr <- as.matrix(mast@assays$RNA@counts)
  mast.expr <- mast.expr[rowSums(mast.expr) >= 1,]
  mast.meta <- mast@meta.data[,c("Patient")];names(mast.meta) <- colnames(mast)
  es.mat <- GSVA::gsva(mast.expr,s.sets,kcdf="Poisson")
  write.table(es.mat, file = paste0("output/mast.gsva.txt",row.names = T, col.names = NA, sep = "\t"))
  pheatmap::pheatmap(es.mat, show_rownames = 1, show_colnames = 0, annotation_col = mast.meta,
                     fontsize_row = 5, filename = "output/mast.gsva.pdf",width=15, height=12)
  # CD14 mono
  CD14mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes=="CD14 Mono.")))
  CD14mono.expr <- as.matrix(CD14mono@assays$RNA@counts)
  CD14mono.expr <- CD14mono.expr[rowSums(CD14mono.expr) >= 1,]
  CD14mono.meta <- CD14mono@meta.data[,c("Patient")];names(CD14mono.meta) <- colnames(CD14mono)
  es.mat <- GSVA::gsva(CD14mono.expr,s.sets,kcdf="Poisson")
  write.table(es.mat, file = paste0("output/CD14mono.gsva.txt",row.names = T, col.names = NA, sep = "\t"))
  pheatmap::pheatmap(es.mat, show_rownames = 1, show_colnames = 0, annotation_col = CD14mono.meta,
                     fontsize_row = 5, filename = "output/CD14mono.gsva.pdf",width=15, height=12)
  # CD16 mono
  CD16mono <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes=="CD16 Mono.")))
  CD16mono.expr <- as.matrix(CD16mono@assays$RNA@counts)
  CD16mono.expr <- CD16mono.expr[rowSums(CD16mono.expr) >= 1,]
  CD16mono.meta <- CD16mono@meta.data[,c("Patient")];names(CD16mono.meta) <- colnames(CD16mono)
  es.mat <- GSVA::gsva(CD16mono.expr,s.sets,kcdf="Poisson")
  write.table(es.mat, file = paste0("output/CD16mono.gsva.txt",row.names = T, col.names = NA, sep = "\t"))
  pheatmap::pheatmap(es.mat, show_rownames = 1, show_colnames = 0, annotation_col = CD16mono.meta,
                     fontsize_row = 5, filename = "output/CD16mono.gsva.pdf",width=15, height=12)
  # DC
  DC <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data, CellTypes=="DCs")))
  DC.expr <- as.matrix(DC@assays$RNA@counts)
  DC.expr <- DC.expr[rowSums(DC.expr) >= 1,]
  DC.meta <- DC@meta.data[,c("Patient")];names(DC.meta) <- colnames(DC)
  es.mat <- GSVA::gsva(DC.expr,s.sets,kcdf="Poisson")
  write.table(es.mat, file = paste0("output/DC.gsva.txt",row.names = T, col.names = NA, sep = "\t"))
  pheatmap::pheatmap(es.mat, show_rownames = 1, show_colnames = 0, annotation_col = DC.meta,
                     fontsize_row = 5, filename = "output/DC.gsva.pdf",width=15, height=12)
}


