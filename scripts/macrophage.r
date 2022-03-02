## ------------------ Mcrophage subset ---------------------------

## ------------------ 1. Cell Annotation ---------------------------

### --------- PDAC subset -----------------
myeloid.pdac <- readRDS("input/mye.with.mast.pdac.sub.rds")
mye.pdac.macro <- subset(myeloid.pdac, cells = rownames(subset(myeloid.pdac@meta.data,CellTypes == "Macrophage")))
mye.pdac.macro$CellTypes <- as.character(mye.pdac.macro$CellTypes)
mye.pdac.macro <- NormalizeData(object = mye.pdac.macro, normalization.method = "LogNormalize", scale.factor = 10000)
mye.pdac.macro <- FindVariableFeatures(mye.pdac.macro,selection.method="vst",nfeatures=2000)
DefaultAssay(mye.pdac.macro) <- "integrated"
mye.pdac.macro <- ScaleData(mye.pdac.macro, features = rownames(mye.pdac.macro))
mye.pdac.macro <- RunPCA(object = mye.pdac.macro, features = VariableFeatures(mye.pdac.macro))
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
PC.num <- PCDeterminators(mye.pdac.macro)
# Find neighbors and clusters with harmony batch correction
mye.pdac.macro <- FindNeighbors(mye.pdac.macro, dims = 1:PC.num, reduction = "pca")
mye.pdac.macro <- FindClusters(mye.pdac.macro, resolution = 0.1, reduction = "pca")
# myeloid <- myeloid %>% RunHarmony("SampleName",plot_convergence = F)
mye.pdac.macro <- RunUMAP(mye.pdac.macro, dims = 1:PC.num, reduction = "pca")
mye.pdac.macro <- RunTSNE(mye.pdac.macro,dims = 1:PC.num)

saveRDS(mye.pdac.macro, file = "input/mye.pdac.macro.sub.dr.rds")

mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "right",
                 plot.title = element_blank(),
                 legend.text = element_text(size = 14))
mycols <- c('#1f77b4', '#39d486', '#43315c', '#00c8ff', '#de4343', '#ff7f0e', '#e377c2')
p.mye.pdac.macro <- DimPlot(mye.pdac.macro,reduction = "tsne",label = T,cols = mycols) + mytheme
ggsave(paste0("figures/celltypes/macrophage/","macrophage.pdf"),p.mye.pdac.macro,width = 3,height = 2.6)
p.mye.pdac.macro.splitBy.state <- DimPlot(mye.pdac.macro,reduction = "tsne",split.by = "State",cols = mycols)
ggsave(paste0("figures/celltypes/macrophage/","macro.splitBy.State.pdf"),p.mye.pdac.macro.splitBy.state,width = 5,height = 2.8)

## ---->> M1 M2 markers
# M1:
# iNOS,IL12,CD64(FcγR1A),CD64(FcγR1B),CD64(FcγR1C),CD80(B7-1),CXCL10,IL23,CXCL9,CXCL11
# CD86(B7-2),IL1A,IL1B,IL6,TNFa,MHCII,CCL5,IRF5,IRF1,CD40,IDO1,KYNU,CCR7
# 
# m1.fea <- c("NOS2","IL12A","FCGR1A","FCGR1B","CD80","CXCL10","IL23A","CXCL9","CXCL11",
#             "CD86","IL1A","IL1B","IL6","TNF","HLA-DPB1","HLA-DMA","HLA-DMB","HLA-DOA",
#             "HLA-DOB","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DQB1-AS1","HLA-DRA","HLA-DRB1","HLA-DRB5",
#             "CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7")
# M2:
# ARG1/2(Arginase),IL10,CD32,CD163,CD23(FCER2),CD200R1,PD-L2(PDCD1LG2),PD-L1(CD274),MARCO,CSF1R
# CD206(MRC1),Il1RA(IL1RN),Il1R2,IL4R,CCL4,CCL13,CCL20,CCL17,CCL18,CCL22,CCL24,LYVE1,VEGFA,VEGFB
# VEGFC,VEGFD,EGF,Cathepsin A(CTSA),CTSB,CSTC,CTSD,TGFB1,TGFB2,TGFB3,MMP14,MMP19,MMP9,CLEC7A
# WNT7b,TNFSF12,TNFSF8,CD276(B7-H3),VTCN1(BH-H4),MSR1(CD204),FN1,IRF4
# 
# m2.fea <- c("ARG1","ARG2","IL10","FCGR2A","CD163","FCER2","CD200R1","PDCD1LG2","CD274","MARCO",
#             "MRC1","IL1RN","IL1R2","IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22",
#             "CCL24","LYVE1","VEGFA","VEGFB","VEGFC","FIGF","FASLG","EGF","CTSA","CTSB","CST3","CTSD",
#             "TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","TNFSF12","TNFSF8","CD276",
#             "VTCN1","MSR1","FN1","IRF4")
# RTM (Tissue resident macrophage):
# RGS1, HLA-DPA1, CD74, C1QC, FCGBP, SEPP1, CSF1R, TREM2, S100A4
# rtm.fea <- c("RGS1","HLA-DPA1","CD74","C1QC","FCGBP","SEPP1","CSF1R","TREM2","S100A4")

# features <- c(m1.fea, m2.fea, rtm.fea)
# features <- unique(features)

DefaultAssay(mye.pdac.macro) <- "RNA"

### Option I:
# M1: HLA-DR, CD11C, CD86, iNOS, pSTAT1
# M2: CD163, CD204, CD206, VEGFA/B, cMAF
# 
p.m1 <- VlnPlot(mye.pdac.macro,features = c("HLA-DRA","HLA-DRB5","HLA-DRB1","ITGAX","CD86","NOS2","STAT1"),ncol = 1,pt.size = 0)
ggsave("p.m1.pdf",p.m1,width = 6,height = 12)
p.m2 <- VlnPlot(mye.pdac.macro,features = c("CD163","MSR1","MRC1","VEGFA","VEGFB","MAF"),ncol = 1,pt.size = 0)
ggsave("p.m2.pdf",p.m2,width = 6,height = 12)
p.rtm <- VlnPlot(mye.pdac.macro,features = c("RGS1","HLA-DPA1","CD74","C1QC","FCGBP","SEPP1","CSF1R","TREM2","S100A4"),ncol = 1,pt.size = 0)
ggsave("p.rmt.pdf",p.rtm,width = 6,height = 16)

### Option II:
m1.fea <- c("HLA-DRA","HLA-DRB5","HLA-DRB1","ITGAX","CD86","NOS2","STAT1")
rtm.fea <- c("RGS1","HLA-DPA1","CD74","C1QC","FCGBP","SEPP1","CSF1R","TREM2","S100A4")
m2.fea <- c("CD40","SPP1","MARCO","APOE","CHIT1","FABP5","CCL18","HLA-DQA1","COX17","LY6E",
            "LAMP1","HAVCR2","SERPINB6","IL18","CCL2","ATF5","CXCL3","VEGFB","SLC2A1","CD163","MSR1","VEGFA","MAF")
features <- c(m1.fea,m2.fea,rtm.fea)

p.macro.dot <- DotPlot(mye.pdac.macro,features = features) + NoLegend() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1,colour = "black")) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlBu")))
ggsave(paste0("figures/celltypes/macrophage/","fea.dot.pdf"),p.macro.dot,width = 20,height = 3)

curr.ids <- c(0:6)
new.ids <- c("M2","RTM","M2","M1","M2","RTM","M2")
mye.pdac.macro$CellTypes <- plyr::mapvalues(mye.pdac.macro$seurat_clusters, from = curr.ids, to = new.ids)
saveRDS(mye.pdac.macro,file = "input/cellSub/mye.pdac.macro.sub.cell.annot.rds")
mycols <- c('#39d486', '#00c8ff', '#de4343')
p.ct <- DimPlot(mye.pdac.macro, reduction = "tsne",cols = mycols,group.by = "CellTypes") + mytheme
ggsave(paste0("figures/celltypes/macrophage/","celltypes.pdf"),p.ct,width = 5,height = 4.4)

p.mye.pdac.macro.celltype.splitBy.state <- DimPlot(mye.pdac.macro,reduction = "tsne",split.by = "State",cols = mycols,group.by = "CellTypes")
ggsave(paste0("figures/celltypes/macrophage/","macro.celltype.splitBy.State.pdf"),p.mye.pdac.macro.celltype.splitBy.state,width = 5,height = 2.8)

## ---->>> cell proportion calculation

mytheme.1 <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                   panel.grid = element_blank(),
                   axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
                   axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
                   legend.position = "top", legend.direction = "horizontal",
                   legend.text = element_text(size = 14),legend.title = element_text(size = 15))

p.pdac.macro.tumor_vs_paratumor <- ggplot(mye.pdac.macro@meta.data, aes(x = State, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mycols) +
  mytheme.1 +
  guides(fill = guide_legend(nrow = 1)) +
  labs(y = "Cell Fraction")
ggsave(paste0("figures/celltypes/macrophage/","mye.pdac.macro.cell.prop_vs_State.pdf"),p.pdac.macro.tumor_vs_paratumor, width = 5,height = 2)


### --------- PBMC subset -----------------
# myeloid.pbmc <- readRDS("input/mye.with.mast.pbmc.sub.rds")
# mye.pbmc.macro <- subset(myeloid.pbmc, cells = rownames(subset(myeloid.pbmc@meta.data,CellTypes == "Macrophage")))
# mye.pbmc.macro$CellTypes <- as.character(mye.pbmc.macro$CellTypes)
# mye.pbmc.macro <- NormalizeData(object = mye.pbmc.macro, normalization.method = "LogNormalize", scale.factor = 10000)
# mye.pbmc.macro <- FindVariableFeatures(mye.pbmc.macro,selection.method="vst",nfeatures=2000)
# DefaultAssay(mye.pbmc.macro) <- "integrated"
# mye.pbmc.macro <- ScaleData(mye.pbmc.macro, features = rownames(mye.pbmc.macro))
# mye.pbmc.macro <- RunPCA(object = mye.pbmc.macro, features = VariableFeatures(mye.pbmc.macro))


## ------------------ 2. Survival Analysis ---------------------------

dataProcess <- function(TCGA.exp.mtx, clinical.data, signatures){
  expr <- read.table(TCGA.exp.mtx, sep = "\t", header = T)
  expr <- expr[!duplicated(expr$Hugo_Symbol),]
  rownames(expr) <- expr$Hugo_Symbol
  expr <- expr[, -c(1,2)]
  colnames(expr) <- gsub("\\.", "-", colnames(expr))
  expr <- data.frame(t(expr))
  colnames(expr) <- gsub("\\.", "-", colnames(expr))
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


# M1
m1.fea <- c("ITGAX","CD86","NOS2","STAT1")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = m1.fea)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "M1 signatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","M1.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","M1.surv.2.pdf"),p.surv.2,width = 4,height = 4)
# RTM
rtm.fea <- c("RGS1","CD74","C1QC","FCGBP","SEPP1","CSF1R","TREM2","S100A4")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = rtm.fea)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "RTM signatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","RTM.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","RTM.surv.2.pdf"),p.surv.2,width = 4,height = 4)

# M2
m2.fea <- c("CD40","SPP1","MARCO","APOE","CHIT1","FABP5","CCL18","COX17","LY6E",
            "LAMP1","HAVCR2","SERPINB6","IL18","CCL2","ATF5","CXCL3","VEGFB","SLC2A1","CD163","MSR1","VEGFA","MAF")
merged.data <- dataProcess(TCGA.exp.mtx = "input/tcga/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
                           clinical.data = "input/tcga/data_clinical_patient.txt",
                           signatures = m2.fea)
mysurv <- survival::survfit(survival::Surv(clin.OS_MONTHS,Dead) ~ group, data = merged.data)

p.Surv <- survminer::ggsurvplot(mysurv)
p.surv <- p.Surv$plot + 
  annotate("text", label = paste0("p = ", signif(as.numeric(survminer::surv_pvalue(mysurv)[2]),digits = 3)),
           x = 9, y = 0.15, size = 5) +
  annotate("text", label = paste0("n = 178"), x = 4, y = 0.22, size = 5) +
  guides(color = guide_legend(title = "M2 signatures")) +
  mytheme.3 +
  scale_color_manual(values = rev(ggsci::pal_aaas(palette = c("default"), alpha = 1)(2)))
p.surv.2 <- p.surv + theme(legend.position = "none")

ggsave(paste0("figures/celltypes/myeloid/survival/","M2.surv.pdf"),p.surv,width = 5.5,height = 4)
ggsave(paste0("figures/celltypes/myeloid/survival/","M2.surv.2.pdf"),p.surv.2,width = 4,height = 4)

rm(list = ls());gc()

## ----------- Visualize Interesting Genes
## for macrophage
mye.pdac.macro <- readRDS("input/cellSub/mye.pdac.macro.sub.cell.annot.rds")

mypal <- c("grey88","DarkCyan")
mytheme <- theme(panel.border = element_rect(size = 1.5,colour = "black"),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   legend.position = "bottom")

f.GLUL <- FeaturePlot(mye.pdac.macro,reduction = "tsne",features = "GLUL",cols = mypal,pt.size = 0.05) + mytheme
ggsave(paste0("figures/celltypes/features/","GLUL.pdf"),f.GLUL,width = 4,height = 4.4)

f.SQSTM1 <- FeaturePlot(mye.pdac.macro,reduction = "tsne",features = "SQSTM1",cols = mycols) + mytheme
ggsave(paste0("figures/celltypes/features/","SQSTM1.pdf"),f.SQSTM1,width = 4,height = 4.4)

library(MySeuratWrappers)
mycols <- c('#1f77b4', '#39d486', '#43315c', '#00c8ff', '#de4343', '#ff7f0e', '#e377c2')
p.vln <- VlnPlot(mye.pdac.macro,features = c("GLUL","SQSTM1"),stacked = T,cols = mycols,
                 pt.size = 0,combine = T,direction = "horizontal",x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/features/","vlnplot.pdf"),p.vln,width = 2,height = 5)

Idents(mye.pdac.macro) <- mye.pdac.macro$seurat_clusters
all.markers <- FindAllMarkers(mye.pdac.macro,assay = "RNA")

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# PRSS1, CTRB1, CLPS, PLA2G1B, PNLIP, CPA1
p.vln <- VlnPlot(mye.pdac.macro,features = c("PRSS1","CTRB1","CLPS","PLA2G1B","PNLIP","CPA1"),
                 stacked = T,cols = mycols,pt.size = 0,combine = T,direction = "horizontal",
                 x.lab = '',y.lab = '') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste0("figures/celltypes/features/","c5.specific.pdf"),p.vln,width = 4,height = 4)

##------ validation in this corhort

## correlation with CD8+T cells
Tsub <- readRDS("input/cellSub/Tsub.cell.annot.rds")
gene.1 <- c("CD8A","CD8B")
CD8T <- subset(Tsub,cells = rownames(subset(Tsub@meta.data,CellTypes=="CD8T")))
Idents(CD8T) <- CD8T$SampleName
saveRDS(CD8T,file = "input/cellSub/CD8T.rds")
avgExp.T <- AverageExpression(CD8T,assays = "RNA",slot = "data",features = gene.1)
avgExp.T <- data.frame(t(avgExp.T$RNA),stringsAsFactors = F)
avgExp.T$Sample <- rownames(avgExp.T)
gene.2 <- c("PRSS1","CTRB1","CLPS","PLA2G1B","PNLIP","CPA1")
Idents(mye.pdac.macro) <- mye.pdac.macro$SampleName
avgExp.T2 <- AverageExpression(mye.pdac.macro,assays = "RNA",slot = "data",features = gene.2)
avgExp.T2 <- data.frame(t(avgExp.T2$RNA),stringsAsFactors = F)
avgExp.T2$Sample <- rownames(avgExp.T2)

dplyr::left_join(avgExp.T2,avgExp.T,by = "Sample") -> avgExp
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL

for(i in 1:length(gene.2)){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[7])) +
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
  ggsave(paste0("figures/celltypes/macrophage/corr.c5/thisCorhort/",colnames(avgExp)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

## correlation with epithelial cells
Epigene <- c("EPCAM","ACTA2","KRT7","KRT8","KRT18","KRT19","CDH1","PRSS1","CTRB2",
                     "REG1A","CLU","MKI67","SPINK1","TFF1","MUC1")
Epi <- readRDS("input/cellSub/Epithelial.rds")
Idents(Epi) <- Epi$SampleName
avgExp.E <- AverageExpression(Epi,assays = "RNA",slot = "data",features = Epigene)
avgExp.E <- data.frame(t(avgExp.E$RNA),stringsAsFactors = F)
avgExp.E$Sample <- rownames(avgExp.E)

avgExp <- dplyr::left_join(avgExp.T2,avgExp.E,by = "Sample")
rownames(avgExp) <- avgExp$Sample
avgExp$Sample <- NULL
colnames(avgExp)[1] <- "PRSS1"

for(i in 1:6){
  p <- ggplot(avgExp,aes_string(x = colnames(avgExp)[i], y = colnames(avgExp)[7])) +
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
  ggsave(paste0("figures/celltypes/macrophage/corr.c5/thisCorhort/Epi/",colnames(avgExp)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}

##------ validation in TCGA
# correlation with CD8T
load("input/tcga/TCGA_36cancers.rnaseq.RData")
PAAD.rnaseq <- TCGA_36cancers.rnaseq$PAAD.rnaseq
in_use <- colnames(PAAD.rnaseq[which(colnames(PAAD.rnaseq) %in% gene.2)])
dt <- log2(PAAD.rnaseq[,c("CD8A",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
dt2 <- dt[-grep("11A",rownames(dt)),]
for(i in 7:2){
  p <- ggplot(dt2,aes_string(x = colnames(dt2)[i], y = colnames(dt2)[1])) +
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
  ggsave(paste0("figures/celltypes/macrophage/corr.c5/TCGA/CD8T/",colnames(dt)[i],"_vs_","CD8A.pdf"),p,width = 5,height = 5)
}

# correlation with epithelial cell
dt <- log2(PAAD.rnaseq[,c("EPCAM",in_use)]+1)
rownames(dt) <- PAAD.rnaseq$bcr_patient_barcode
dt2 <- dt[-grep("11A",rownames(dt)),]
for(i in 7:2){
  p <- ggplot(dt2,aes_string(x = colnames(dt2)[i], y = colnames(dt2)[1])) +
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
  ggsave(paste0("figures/celltypes/macrophage/corr.c5/TCGA/Epi/",colnames(dt)[i],"_vs_","EPCAM.pdf"),p,width = 5,height = 5)
}

## visualize macrophage subset splitBy samples
macro <- readRDS("input/cellSub/mye.pdac.macro.sub.cell.annot.rds")
curr.ids <- c(0:6)
new.ids <- c("M2","RTM","M2","M1","M2","GLUL-SQSTM1-RTM","M2")
macro$CellTypes.2 <- plyr::mapvalues(macro$seurat_clusters,from = curr.ids,to = new.ids)
p.macro <- DimPlot(macro, reduction = "tsne", label = T, pt.size = 0.2, split.by = "SampleName",
                   group.by = "CellTypes.2",ncol = 5)
ggsave("figures/celltypes/macrophage/splitBySampleName.pdf",p.macro,width = 15,height = 15)

# unused
# avgExprPerSample <- function(object,gene){
#   expr <- GetAssayData(object, assay = "RNA", splot = "data")[gene,]
#   expr <- as.data.frame(expr, stringAsFactors = F)
#   expr.T <- data.frame(t(expr))
#   expr.T$Sample <- plyr::mapvalues(rownames(expr.T),from = colnames(object), to = object$SampleName)
#   colnames(expr.T) <- c(gene,"Sample")
#   
#   patient <- NULL
#   for(j in 1:(ncol(expr.T)-1)){
#     patient[j] <- data.frame(tapply(expr.T[,j],expr.T$Sample,sum))
#   }
#   names(patient) <- colnames(expr.T)[-ncol(expr.T)]
#   patient <- data.frame(patient)
#   rownames(patient) <- names(tapply(expr.T[,1],expr.T$Sample,sum))
#   Patient_sample <- data.frame(table(expr.T$Sample),stringsAsFactors = F)
#   patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1, to = Patient_sample$Freq)
#   patient <- data.frame(patient, stringsAsFactors = F)
#   patient$Sample <- as.character(patient$Sample)
#   for(k in 1:(ncol(patient)-1)){
#     for(l in 1:nrow(patient)){
#       patient[l,k] <- patient[l,k]/patient$Sample[l]
#     }
#   }
#   patient$Sample <- NULL
#   colnames(patient) <- gene
#   return(patient)
# }


## ----------- 3.Differential Gene Analysis ---------------
mye.pdac.macro <- readRDS("input/mye.pdac.macro.sub.cell.annot.rds")
library(ggplot2)
library(ggrepel)
mytheme<- theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 16),
                axis.text = element_text(size = 12),
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 12))

# PDAC (tumor vs. para-tumor && tumor.tissue vs. paratumor.tissue && Tissue vs. PBMC)
mye.pdac.macro.tiss <- subset(mye.pdac.macro, cells = rownames(subset(mye.pdac.macro@meta.data, Tissue == "Tissue")))
Idents(mye.pdac.macro.tiss) <- mye.pdac.macro.tiss$CellTypes
saveRDS(mye.pdac.macro.tiss, file = "input/mye.with.mast.pdac.macro.tissue.sub.rds")
mye.pdac.macro.tumor <- subset(mye.pdac.macro, cells = rownames(subset(mye.pdac.macro@meta.data, State == "tumor")))

## --- M1 ---
## 1.tumor vs para-tumor
Idents(mye.pdac.macro) <- mye.pdac.macro$CellTypes
mye.pdac.m1.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro, ident.1 = "tumor",group.by = "State",subset.ident = "M1",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.m1.tumor_vs_paratumor <- mye.pdac.m1.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.m1.tumor_vs_paratumor)),]
mye.pdac.m1.tumor_vs_paratumor <- mye.pdac.m1.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.m1.tumor_vs_paratumor)),]
mye.pdac.m1.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.m1.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.m1.tumor_vs_paratumor <- mye.pdac.m1.tumor_vs_paratumor[order(mye.pdac.m1.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.m1.tumor_vs_paratumor,10),tail(mye.pdac.m1.tumor_vs_paratumor,10))

p.mye.pdac.m1.tumor_vs_paratumor <- ggplot(mye.pdac.m1.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.m1.tumor_vs_paratumor.pdf"),p.mye.pdac.m1.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tiss.m1.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro.tiss,ident.1 = "tumor",group.by = "State",subset.ident = "M1")
# remove ribosomal and mitochondrial genes
mye.pdac.tiss.m1.tumor_vs_paratumor <- mye.pdac.tiss.m1.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tiss.m1.tumor_vs_paratumor)),]
mye.pdac.tiss.m1.tumor_vs_paratumor <- mye.pdac.tiss.m1.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tiss.m1.tumor_vs_paratumor)),]
mye.pdac.tiss.m1.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tiss.m1.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tiss.m1.tumor_vs_paratumor <- mye.pdac.tiss.m1.tumor_vs_paratumor[order(mye.pdac.tiss.m1.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tiss.m1.tumor_vs_paratumor,10), tail(mye.pdac.tiss.m1.tumor_vs_paratumor,10))
p.mye.pdac.tiss.m1.tumor_vs_paratumor <- ggplot(mye.pdac.tiss.m1.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.tissue.m1.tumor_vs_paratumor.pdf"),p.mye.pdac.tiss.m1.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.m1.tiss_vs_pbmc <- FindMarkers(mye.pdac.macro.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "M1",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.m1.tiss_vs_pbmc <- mye.pdac.m1.tiss_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.m1.tiss_vs_pbmc)),]
mye.pdac.m1.tiss_vs_pbmc <- mye.pdac.m1.tiss_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.m1.tiss_vs_pbmc)),]
mye.pdac.m1.tiss_vs_pbmc$Significance <- ifelse(mye.pdac.m1.tiss_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.m1.tiss_vs_pbmc <- mye.pdac.m1.tiss_vs_pbmc[order(mye.pdac.m1.tiss_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.m1.tiss_vs_pbmc,10), tail(mye.pdac.m1.tiss_vs_pbmc,10))
p.mye.pdac.m1.tiss_vs_pbmc <- ggplot(mye.pdac.m1.tiss_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.m1.tissue_vs_pbmc.pdf"),p.mye.pdac.m1.tiss_vs_pbmc,width = 7,height = 5)

save(mye.pdac.m1.tumor_vs_paratumor, 
     mye.pdac.tiss.m1.tumor_vs_paratumor,
     mye.pdac.m1.tiss_vs_pbmc, 
     file = "output/macrophage/DEgenes/mye.pdac.m1.DEgenes.rda")

## --- M2 ---
## 1.tumor vs para-tumor
mye.pdac.m2.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro, ident.1 = "tumor",group.by = "State",subset.ident = "M2",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.m2.tumor_vs_paratumor <- mye.pdac.m2.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.m2.tumor_vs_paratumor)),]
mye.pdac.m2.tumor_vs_paratumor <- mye.pdac.m2.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.m2.tumor_vs_paratumor)),]
mye.pdac.m2.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.m2.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.m2.tumor_vs_paratumor <- mye.pdac.m2.tumor_vs_paratumor[order(mye.pdac.m2.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.m2.tumor_vs_paratumor,10),tail(mye.pdac.m2.tumor_vs_paratumor,10))

p.mye.pdac.m2.tumor_vs_paratumor <- ggplot(mye.pdac.m2.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.m2.tumor_vs_paratumor.pdf"),p.mye.pdac.m2.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tiss.m2.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro.tiss,ident.1 = "tumor",group.by = "State",subset.ident = "M2")
# remove ribosomal and mitochondrial genes
mye.pdac.tiss.m2.tumor_vs_paratumor <- mye.pdac.tiss.m2.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tiss.m2.tumor_vs_paratumor)),]
mye.pdac.tiss.m2.tumor_vs_paratumor <- mye.pdac.tiss.m2.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tiss.m2.tumor_vs_paratumor)),]
mye.pdac.tiss.m2.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tiss.m2.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tiss.m2.tumor_vs_paratumor <- mye.pdac.tiss.m2.tumor_vs_paratumor[order(mye.pdac.tiss.m2.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tiss.m2.tumor_vs_paratumor,10), tail(mye.pdac.tiss.m2.tumor_vs_paratumor,10))
p.mye.pdac.tiss.m2.tumor_vs_paratumor <- ggplot(mye.pdac.tiss.m2.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.tissue.m2.tumor_vs_paratumor.pdf"),p.mye.pdac.tiss.m2.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.m2.tiss_vs_pbmc <- FindMarkers(mye.pdac.macro.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "M2",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.m2.tiss_vs_pbmc <- mye.pdac.m2.tiss_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.m2.tiss_vs_pbmc)),]
mye.pdac.m2.tiss_vs_pbmc <- mye.pdac.m2.tiss_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.m2.tiss_vs_pbmc)),]
mye.pdac.m2.tiss_vs_pbmc$Significance <- ifelse(mye.pdac.m2.tiss_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.m2.tiss_vs_pbmc <- mye.pdac.m2.tiss_vs_pbmc[order(mye.pdac.m2.tiss_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.m2.tiss_vs_pbmc,10), tail(mye.pdac.m2.tiss_vs_pbmc,10))
p.mye.pdac.m2.tiss_vs_pbmc <- ggplot(mye.pdac.m2.tiss_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.m2.tissue_vs_pbmc.pdf"),p.mye.pdac.m2.tiss_vs_pbmc,width = 7,height = 5)

save(mye.pdac.m2.tumor_vs_paratumor, 
     mye.pdac.tiss.m2.tumor_vs_paratumor,
     mye.pdac.m2.tiss_vs_pbmc, 
     file = "output/macrophage/DEgenes/mye.pdac.m2.DEgenes.rda")


## --- RTM ---
## 1.tumor vs para-tumor
mye.pdac.rtm.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro, ident.1 = "tumor",group.by = "State",subset.ident = "RTM",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.rtm.tumor_vs_paratumor <- mye.pdac.rtm.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.rtm.tumor_vs_paratumor)),]
mye.pdac.rtm.tumor_vs_paratumor <- mye.pdac.rtm.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.rtm.tumor_vs_paratumor)),]
mye.pdac.rtm.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.rtm.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.rtm.tumor_vs_paratumor <- mye.pdac.rtm.tumor_vs_paratumor[order(mye.pdac.rtm.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.rtm.tumor_vs_paratumor,10),tail(mye.pdac.rtm.tumor_vs_paratumor,10))

p.mye.pdac.rtm.tumor_vs_paratumor <- ggplot(mye.pdac.rtm.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.rtm.tumor_vs_paratumor.pdf"),p.mye.pdac.rtm.tumor_vs_paratumor,width = 7,height = 5)

## 2.tumor.tissue vs. paratumor.tissue
mye.pdac.tiss.rtm.tumor_vs_paratumor <- FindMarkers(mye.pdac.macro.tiss,ident.1 = "tumor",group.by = "State",subset.ident = "RTM")
# remove ribosomal and mitochondrial genes
mye.pdac.tiss.rtm.tumor_vs_paratumor <- mye.pdac.tiss.rtm.tumor_vs_paratumor[!grepl("^RP[SL]",rownames(mye.pdac.tiss.rtm.tumor_vs_paratumor)),]
mye.pdac.tiss.rtm.tumor_vs_paratumor <- mye.pdac.tiss.rtm.tumor_vs_paratumor[!grepl("^MT-",rownames(mye.pdac.tiss.rtm.tumor_vs_paratumor)),]
mye.pdac.tiss.rtm.tumor_vs_paratumor$Significance <- ifelse(mye.pdac.tiss.rtm.tumor_vs_paratumor$p_val < 0.01, TRUE, FALSE)
mye.pdac.tiss.rtm.tumor_vs_paratumor <- mye.pdac.tiss.rtm.tumor_vs_paratumor[order(mye.pdac.tiss.rtm.tumor_vs_paratumor[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.tiss.rtm.tumor_vs_paratumor,10), tail(mye.pdac.tiss.rtm.tumor_vs_paratumor,10))
p.mye.pdac.tiss.rtm.tumor_vs_paratumor <- ggplot(mye.pdac.tiss.rtm.tumor_vs_paratumor, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.tissue.rtm.tumor_vs_paratumor.pdf"),p.mye.pdac.tiss.rtm.tumor_vs_paratumor,width = 7,height = 5)

## 3.Tissue vs PBMC
mye.pdac.rtm.tiss_vs_pbmc <- FindMarkers(mye.pdac.macro.tumor, ident.1 = "Tissue",group.by = "Tissue",subset.ident = "RTM",logfc.threshold = 0.1)
# remove ribosomal and mitochondrial genes
mye.pdac.rtm.tiss_vs_pbmc <- mye.pdac.rtm.tiss_vs_pbmc[!grepl("^RP[SL]",rownames(mye.pdac.rtm.tiss_vs_pbmc)),]
mye.pdac.rtm.tiss_vs_pbmc <- mye.pdac.rtm.tiss_vs_pbmc[!grepl("^MT-",rownames(mye.pdac.rtm.tiss_vs_pbmc)),]
mye.pdac.rtm.tiss_vs_pbmc$Significance <- ifelse(mye.pdac.rtm.tiss_vs_pbmc$p_val < 0.01, TRUE, FALSE)
mye.pdac.rtm.tiss_vs_pbmc <- mye.pdac.rtm.tiss_vs_pbmc[order(mye.pdac.rtm.tiss_vs_pbmc[,2],decreasing = T),]
top10 <- rbind(head(mye.pdac.rtm.tiss_vs_pbmc,10), tail(mye.pdac.rtm.tiss_vs_pbmc,10))
p.mye.pdac.rtm.tiss_vs_pbmc <- ggplot(mye.pdac.rtm.tiss_vs_pbmc, aes(avg_logFC,-log10(p_val))) +
  geom_point(aes(col = Significance)) +
  scale_color_manual(values = c("#BEBBBD","#C72719")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
                  box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
                  segment.color = "black", show.legend = F) +
  theme_classic() + ylab("-log10(P value)") + mytheme
ggsave(paste0("figures/celltypes/macrophage/DEgenes/","mye.pdac.rtm.tissue_vs_pbmc.pdf"),p.mye.pdac.rtm.tiss_vs_pbmc,width = 7,height = 5)

save(mye.pdac.rtm.tumor_vs_paratumor, 
     mye.pdac.tiss.rtm.tumor_vs_paratumor,
     mye.pdac.rtm.tiss_vs_pbmc, 
     file = "output/macrophage/DEgenes/mye.pdac.rtm.DEgenes.rda")


#### >>>>>>>>>>>>> PBMC (PDAC vs Normal)
# mye.pbmc.macro <- subset(mye.pdac.macro, cells = rownames(subset(mye.pdac.macro@meta.data, Tissue == "PBMC")))
# saveRDS(mye.pbmc.macro,file = "input/mye.with.mast.macro.pbmc.sub.rds")
# 
# Idents(mye.pbmc.macro) <- mye.pbmc.macro$CellTypes
# # M1
# mye.pbmc.m1.pdac_vs_normal <- FindMarkers(mye.pbmc.macro, ident.1 = "PDAC", group.by = "Patient", subset.ident = "M1",logfc.threshold = 0.1)
# # remove ribosomal and mitochondrial genes
# mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.macro.pdac_vs_normal)),]
# mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.macro.pdac_vs_normal)),]
# mye.pbmc.macro.pdac_vs_normal$Significance <- ifelse(mye.pbmc.macro.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
# mye.pbmc.macro.pdac_vs_normal <- mye.pbmc.macro.pdac_vs_normal[order(mye.pbmc.macro.pdac_vs_normal[,2],decreasing = T),]
# top10 <- rbind(head(mye.pbmc.macro.pdac_vs_normal,10), tail(mye.pbmc.macro.pdac_vs_normal,10))
# p.mye.pbmc.macro.pdac_vs_normal <- ggplot(mye.pbmc.macro.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
#   geom_point(aes(col = Significance)) +
#   scale_color_manual(values = c("#BEBBBD","#C72719")) +
#   geom_text_repel(data = top10, aes(label = rownames(top10)),
#                   box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
#                   segment.color = "black", show.legend = F) +
#   theme_classic() + ylab("-log10(P value)") + mytheme.4
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.macro.pdac_vs_normal.pdf"),p.mye.pbmc.macro.pdac_vs_normal,width = 7,height = 5)
# save(mye.pbmc.macro.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.macro.DEgenes.rda")
# 
# # M2
# mye.pbmc.mast.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "Mast cells", logfc.threshold = 0.1)
# # remove ribosomal and mitochondrial genes
# mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.mast.pdac_vs_normal)),]
# mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.mast.pdac_vs_normal)),]
# mye.pbmc.mast.pdac_vs_normal$Significance <- ifelse(mye.pbmc.mast.pdac_vs_normal$p_val < 0.05, TRUE, FALSE)
# mye.pbmc.mast.pdac_vs_normal <- mye.pbmc.mast.pdac_vs_normal[order(mye.pbmc.mast.pdac_vs_normal[,2],decreasing = T),]
# top10 <- rbind(head(mye.pbmc.mast.pdac_vs_normal,10), tail(mye.pbmc.mast.pdac_vs_normal,10))
# p.mye.pbmc.mast.pdac_vs_normal <- ggplot(mye.pbmc.mast.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
#   geom_point(aes(col = Significance)) +
#   scale_color_manual(values = c("#BEBBBD","#C72719")) +
#   geom_text_repel(data = top10, aes(label = rownames(top10)),
#                   box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
#                   segment.color = "black", show.legend = F) +
#   theme_classic() + ylab("-log10(P value)") + mytheme.4
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.mast.pdac_vs_normal.pdf"),p.mye.pbmc.mast.pdac_vs_normal,width = 7,height = 5)
# save(mye.pbmc.mast.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.mast.DEgenes.rda")
# 
# # RTM
# mye.pbmc.CD14.mono.pdac_vs_normal <- FindMarkers(myeloid.pbmc, ident.1 = "PDAC", group.by = "Patient", subset.ident = "CD14 Mono.", logfc.threshold = 0.1)
# # remove ribosomal and mitochondrial genes
# mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[!grepl("^RP[SL]",rownames(mye.pbmc.CD14.mono.pdac_vs_normal)),]
# mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[!grepl("^MT-",rownames(mye.pbmc.CD14.mono.pdac_vs_normal)),]
# mye.pbmc.CD14.mono.pdac_vs_normal$Significance <- ifelse(mye.pbmc.CD14.mono.pdac_vs_normal$p_val < 0.01, TRUE, FALSE)
# mye.pbmc.CD14.mono.pdac_vs_normal <- mye.pbmc.CD14.mono.pdac_vs_normal[order(mye.pbmc.CD14.mono.pdac_vs_normal[,2],decreasing = T),]
# top10 <- rbind(head(mye.pbmc.CD14.mono.pdac_vs_normal,10), tail(mye.pbmc.CD14.mono.pdac_vs_normal,10))
# p.mye.pbmc.CD14.mono.pdac_vs_normal <- ggplot(mye.pbmc.CD14.mono.pdac_vs_normal, aes(avg_logFC,-log10(p_val))) +
#   geom_point(aes(col = Significance)) +
#   scale_color_manual(values = c("#BEBBBD","#C72719")) +
#   scale_y_continuous(limits = c(0,50)) +
#   geom_text_repel(data = top10, aes(label = rownames(top10)),
#                   box.padding = unit(0.5,'lines'), point.padding = unit(0.8,'lines'),
#                   segment.color = "black", show.legend = F) +
#   theme_classic() + ylab("-log10(P value)") + mytheme.4
# ggsave(paste0("figures/celltypes/myeloid/DEgenes/","mye.pbmc.CD14.mono.pdac_vs_normal.pdf"),p.mye.pbmc.CD14.mono.pdac_vs_normal,width = 7,height = 5)
# save(mye.pbmc.CD14.mono.pdac_vs_normal, file = "output/myeloid/DEgenes/mye.pbmc.CD14.mono.DEgenes.rda")

## ----------- 4.Gene Enrichment Anlysis (GSEA/KEGG/GSVA)
rm(list = ls());gc()
## ---- KEGG
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

#### ---->>> PDAC (tumor vs. para-tumor && tumor.tissue vs. paratumor.tissue && Tissue vs. PBMC) 
# --->>> 1.Tumor vs. Para-Tumor <<<---
# M1
load("output/macrophage/DEgenes/mye.pdac.m1.DEgenes.rda")
up.genes <- mye.pdac.m1.tumor_vs_paratumor[mye.pdac.m1.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M1.pdac.up",row.names = F)
down.genes <- mye.pdac.m1.tumor_vs_paratumor[mye.pdac.m1.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M1.pdac.down",row.names = F, append = T)

# M2
load("output/macrophage/DEgenes/mye.pdac.m2.DEgenes.rda")
up.genes <- mye.pdac.m2.tumor_vs_paratumor[mye.pdac.m2.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M2.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.m2.tumor_vs_paratumor[mye.pdac.m2.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M2.pdac.down",row.names = F, append = T)

# RTM
load("output/macrophage/DEgenes/mye.pdac.rtm.DEgenes.rda")
up.genes <- mye.pdac.rtm.tumor_vs_paratumor[mye.pdac.rtm.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "RTM.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.rtm.tumor_vs_paratumor[mye.pdac.rtm.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "RTM.pdac.down",row.names = F, append = T)

# --->>> 2.tumor.tissue vs. paratumor.tissue <<<---
# M1
up.genes <- mye.pdac.tiss.m1.tumor_vs_paratumor[mye.pdac.tiss.m1.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M1.pdac.up",row.names = F)
down.genes <- mye.pdac.tiss.m1.tumor_vs_paratumor[mye.pdac.tiss.m1.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M1.pdac.down",row.names = F, append = T)
# M2
up.genes <- mye.pdac.tiss.m2.tumor_vs_paratumor[mye.pdac.tiss.m2.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M2.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.tiss.m2.tumor_vs_paratumor[mye.pdac.tiss.m2.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "M2.pdac.down",row.names = F, append = T)
# RTM
up.genes <- mye.pdac.tiss.rtm.tumor_vs_paratumor[mye.pdac.tiss.rtm.tumor_vs_paratumor$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "RTM.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.tiss.rtm.tumor_vs_paratumor[mye.pdac.tiss.rtm.tumor_vs_paratumor$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.tissue.macro.tumor_vs_paratumor.diff.pathway.xlsx",sheetName = "RTM.pdac.down",row.names = F, append = T)

# --->>> 3.Tissue vs. PBMC <<<---
# M1
up.genes <- mye.pdac.m1.tiss_vs_pbmc[mye.pdac.m1.tiss_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
openxlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "M1.pdac.up",row.names = F)
down.genes <- mye.pdac.m1.tiss_vs_pbmc[mye.pdac.m1.tiss_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "M1.pdac.down",row.names = F, append = T)
# M2
up.genes <- mye.pdac.m2.tiss_vs_pbmc[mye.pdac.m2.tiss_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "M2.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.m2.tiss_vs_pbmc[mye.pdac.m2.tiss_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "M2.pdac.down",row.names = F, append = T)
# RTM
up.genes <- mye.pdac.rtm.tiss_vs_pbmc[mye.pdac.rtm.tiss_vs_pbmc$avg_logFC > 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(up.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "RTM.pdac.up",row.names = F, append = T)
down.genes <- mye.pdac.rtm.tiss_vs_pbmc[mye.pdac.rtm.tiss_vs_pbmc$avg_logFC < 0,]
kegg.result <- EnrichmentAnalysis(DEgenes = rownames(down.genes),enrich.method = "KEGG")
xlsx::write.xlsx(kegg.result@result, file = "output/macrophage/pathway/mye.pdac.macro.tissue_vs_pbmc.diff.pathway.xlsx",sheetName = "RTM.pdac.down",row.names = F, append = T)


# --->>> visualize selected pathway <<<---

