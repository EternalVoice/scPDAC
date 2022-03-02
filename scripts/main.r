# GSE155698 -- Main script
library(Seurat)
library(dplyr)

AdjNorm_TISSUE_1 <- Read10X("data/AdjNorm_TISSUE_1/filtered_feature_bc_matrix/")
AdjNorm_TISSUE_1 <- CreateSeuratObject(AdjNorm_TISSUE_1, min.cells = 3)
AdjNorm_TISSUE_1$SampleName <- "AdjNorm_TISSUE_1"
AdjNorm_TISSUE_1$Group <- "PDACII-AdjNorm"
AdjNorm_TISSUE_1$ClinicalStage <- "II"
AdjNorm_TISSUE_1$Class <- "Paracancerous"
AdjNorm_TISSUE_1$Tissue <- "Tissue"
AdjNorm_TISSUE_1$Patient <- "PDAC"
AdjNorm_TISSUE_1$Age <- 63
AdjNorm_TISSUE_1$Gender <- "Female"
AdjNorm_TISSUE_1$GEO <- "GSE155698"
saveRDS(AdjNorm_TISSUE_1, file="AdjNorm_TISSUE_1.rds")

AdjNorm_TISSUE_2 <- Read10X("data/AdjNorm_TISSUE_2/filtered_feature_bc_matrix/")
AdjNorm_TISSUE_2 <- CreateSeuratObject(AdjNorm_TISSUE_2, min.cells = 3)
AdjNorm_TISSUE_2$SampleName <- "AdjNorm_TISSUE_2"
AdjNorm_TISSUE_2$Group <- "PDAC-AdjNorm"
AdjNorm_TISSUE_2$ClinicalStage <- "II"
AdjNorm_TISSUE_2$Class <- "Paracancerous"
AdjNorm_TISSUE_2$Tissue <- "Tissue"
AdjNorm_TISSUE_2$Patient <- "PDAC"
AdjNorm_TISSUE_2$Age <- 46
AdjNorm_TISSUE_2$Gender <- "Female"
AdjNorm_TISSUE_2$GEO <- "GSE155698"
saveRDS(AdjNorm_TISSUE_2, file="AdjNorm_TISSUE_2.rds")

AdjNorm_TISSUE_3 <- Read10X("data/AdjNorm_TISSUE_3/filtered_feature_bc_matrix/")
AdjNorm_TISSUE_3 <- CreateSeuratObject(AdjNorm_TISSUE_3, min.cells = 3)
AdjNorm_TISSUE_3$SampleName <- "AdjNorm_TISSUE_3"
AdjNorm_TISSUE_3$Group <- "PDACII-AdjNorm"
AdjNorm_TISSUE_3$ClinicalStage <- "II"
AdjNorm_TISSUE_3$Class <- "Paracancerous"
AdjNorm_TISSUE_3$Tissue <- "Tissue"
AdjNorm_TISSUE_3$Patient <- "PDAC"
AdjNorm_TISSUE_3$Age <- 68
AdjNorm_TISSUE_3$Gender <- "Male"
AdjNorm_TISSUE_3$GEO <- "GSE155698"
saveRDS(AdjNorm_TISSUE_3, file="AdjNorm_TISSUE_3.rds")

Healthy_PBMC_1 <- Read10X("data/Healthy_PBMC_1/filtered_feature_bc_matrix/")
Healthy_PBMC_1 <- CreateSeuratObject(Healthy_PBMC_1, min.cells = 3)
Healthy_PBMC_1$SampleName <- "Healthy_PBMC_1"
Healthy_PBMC_1$Group <- "Healthy"
Healthy_PBMC_1$Patient <- "Normal"
Healthy_PBMC_1$Tissue <- "PBMC"
Healthy_PBMC_1$Age <- 63
Healthy_PBMC_1$Gender <- "Female"
Healthy_PBMC_1$GEO <- "GSE155698"
saveRDS(Healthy_PBMC_1, file="Healthy_PBMC_1.rds")

Healthy_PBMC_2 <- Read10X("data/Healthy_PBMC_2/filtered_feature_bc_matrix/")
Healthy_PBMC_2 <- CreateSeuratObject(Healthy_PBMC_2, min.cells = 3)
Healthy_PBMC_2$SampleName <- "Healthy_PBMC_2"
Healthy_PBMC_2$Group <- "Healthy"
Healthy_PBMC_2$Patient <- "Normal"
Healthy_PBMC_2$Tissue <- "PBMC"
Healthy_PBMC_2$Age <- 64
Healthy_PBMC_2$Gender <- "Female"
Healthy_PBMC_2$GEO <- "GSE155698"
saveRDS(Healthy_PBMC_2, file="Healthy_PBMC_2.rds")

Healthy_PBMC_3 <- Read10X("data/Healthy_PBMC_3/filtered_feature_bc_matrix/")
Healthy_PBMC_3 <- CreateSeuratObject(Healthy_PBMC_3, min.cells = 3)
Healthy_PBMC_3$SampleName <- "Healthy_PBMC_3"
Healthy_PBMC_3$Group <- "Healthy"
Healthy_PBMC_3$Patient <- "Normal"
Healthy_PBMC_3$Tissue <- "PBMC"
Healthy_PBMC_3$Age <- 71
Healthy_PBMC_3$Gender <- "Male"
Healthy_PBMC_3$GEO <- "GSE155698"
saveRDS(Healthy_PBMC_3, file="Healthy_PBMC_3.rds")

Healthy_PBMC_4 <- Read10X("data/Healthy_PBMC_4/filtered_feature_bc_matrix/")
Healthy_PBMC_4 <- CreateSeuratObject(Healthy_PBMC_4, min.cells = 3)
Healthy_PBMC_4$SampleName <- "Healthy_PBMC_4"
Healthy_PBMC_4$Group <- "Healthy"
Healthy_PBMC_4$Patient <- "Normal"
Healthy_PBMC_4$Tissue <- "PBMC"
Healthy_PBMC_4$Age <- 63
Healthy_PBMC_4$Gender <- "Male"
Healthy_PBMC_4$GEO <- "GSE155698"
saveRDS(Healthy_PBMC_4, file="Healthy_PBMC_4.rds")

PDAC_PBMC_1 <- Read10X("data/PDAC_PBMC_1/filtered_feature_bc_matrix/")
PDAC_PBMC_1 <- CreateSeuratObject(PDAC_PBMC_1, min.cells = 3)
PDAC_PBMC_1$SampleName <- "PDAC_PBMC_1"
PDAC_PBMC_1$Group <- "PDACIV-Met-PBMC"
PDAC_PBMC_1$ClinicalStage <- "IV"
PDAC_PBMC_1$Class <- "Metastatic"
PDAC_PBMC_1$Tissue <- "PBMC"
PDAC_PBMC_1$Patient <- "PDAC"
PDAC_PBMC_1$Age <- 71
PDAC_PBMC_1$Gender <- "Female"
PDAC_PBMC_1$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_1, file="PDAC_PBMC_1.rds")

PDAC_PBMC_2 <- Read10X("data/PDAC_PBMC_2/filtered_feature_bc_matrix/")
PDAC_PBMC_2 <- CreateSeuratObject(PDAC_PBMC_2, min.cells = 3)
PDAC_PBMC_2$SampleName <- "PDAC_PBMC_2"
PDAC_PBMC_2$Group <- "PDACIV-Met-PBMC"
PDAC_PBMC_2$ClinicalStage <- "IV"
PDAC_PBMC_2$Class <- "Metastatic"
PDAC_PBMC_2$Tissue <- "PBMC"
PDAC_PBMC_2$Patient <- "PDAC"
PDAC_PBMC_2$Age <- 67
PDAC_PBMC_2$Gender <- "Male"
PDAC_PBMC_2$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_2, file="PDAC_PBMC_2.rds")

PDAC_PBMC_3 <- Read10X("data/PDAC_PBMC_3/filtered_feature_bc_matrix/")
PDAC_PBMC_3 <- CreateSeuratObject(PDAC_PBMC_3, min.cells = 3)
PDAC_PBMC_3$SampleName <- "PDAC_PBMC_3"
PDAC_PBMC_3$Group <- "PDACIV-Met-PBMC"
PDAC_PBMC_3$ClinicalStage <- "IV"
PDAC_PBMC_3$Class <- "Metastatic"
PDAC_PBMC_3$Tissue <- "PBMC"
PDAC_PBMC_3$Patient <- "PDAC"
PDAC_PBMC_3$Age <- 53
PDAC_PBMC_3$Gender <- "Female"
PDAC_PBMC_3$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_3, file="PDAC_PBMC_3.rds")

PDAC_PBMC_4 <- Read10X("data/PDAC_PBMC_4/filtered_feature_bc_matrix/")
PDAC_PBMC_4 <- CreateSeuratObject(PDAC_PBMC_4, min.cells = 3)
PDAC_PBMC_4$SampleName <- "PDAC_PBMC_4"
PDAC_PBMC_4$Group <- "PDACIII-Pri-PBMC"
PDAC_PBMC_4$ClinicalStage <- "III"
PDAC_PBMC_4$Class <- "Primary"
PDAC_PBMC_4$Tissue <- "PBMC"
PDAC_PBMC_4$Patient <- "PDAC"
PDAC_PBMC_4$Age <- 62
PDAC_PBMC_4$Gender <- "Female"
PDAC_PBMC_4$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_4, file="PDAC_PBMC_4.rds")

PDAC_PBMC_5 <- Read10X("data/PDAC_PBMC_5/filtered_feature_bc_matrix/")
PDAC_PBMC_5 <- CreateSeuratObject(PDAC_PBMC_5, min.cells = 3)
PDAC_PBMC_5$SampleName <- "PDAC_PBMC_5"
PDAC_PBMC_5$Group <- "PDACIII-Pri-PBMC"
PDAC_PBMC_5$ClinicalStage <- "III"
PDAC_PBMC_5$Class <- "Primary"
PDAC_PBMC_5$Tissue <- "PBMC"
PDAC_PBMC_5$Patient <- "PDAC"
PDAC_PBMC_5$Age <- 66
PDAC_PBMC_5$Gender <- "Male"
PDAC_PBMC_5$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_5, file="PDAC_PBMC_5.rds")

PDAC_PBMC_6 <- Read10X("data/PDAC_PBMC_6/filtered_feature_bc_matrix/")
PDAC_PBMC_6 <- CreateSeuratObject(PDAC_PBMC_6, min.cells = 3)
PDAC_PBMC_6$SampleName <- "PDAC_PBMC_6"
PDAC_PBMC_6$Group <- "PDACIII-Pri-PBMC"
PDAC_PBMC_6$ClinicalStage <- "III"
PDAC_PBMC_6$Class <- "Primary"
PDAC_PBMC_6$Tissue <- "PBMC"
PDAC_PBMC_6$Patient <- "PDAC"
PDAC_PBMC_6$Age <- 63
PDAC_PBMC_6$Gender <- "Male"
PDAC_PBMC_6$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_6, file="PDAC_PBMC_6.rds")

PDAC_PBMC_7 <- Read10X("data/PDAC_PBMC_7/filtered_feature_bc_matrix/")
PDAC_PBMC_7 <- CreateSeuratObject(PDAC_PBMC_7, min.cells = 3)
PDAC_PBMC_7$SampleName <- "PDAC_PBMC_7"
PDAC_PBMC_7$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_7$ClinicalStage <- "II"
PDAC_PBMC_7$Class <- "Primary"
PDAC_PBMC_7$Tissue <- "PBMC"
PDAC_PBMC_7$Patient <- "PDAC"
PDAC_PBMC_7$Age <- 49
PDAC_PBMC_7$Gender <- "Female"
PDAC_PBMC_7$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_7, file="PDAC_PBMC_7.rds")

PDAC_PBMC_8 <- Read10X("data/PDAC_PBMC_8/filtered_feature_bc_matrix/")
PDAC_PBMC_8 <- CreateSeuratObject(PDAC_PBMC_8, min.cells = 3)
PDAC_PBMC_8$SampleName <- "PDAC_PBMC_8"
PDAC_PBMC_8$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_8$ClinicalStage <- "II"
PDAC_PBMC_8$Class <- "Primary"
PDAC_PBMC_8$Tissue <- "PBMC"
PDAC_PBMC_8$Patient <- "PDAC"
PDAC_PBMC_8$Age <- 69
PDAC_PBMC_8$Gender <- "Female"
PDAC_PBMC_8$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_8, file="PDAC_PBMC_8.rds")

PDAC_PBMC_9 <- Read10X("data/PDAC_PBMC_9/filtered_feature_bc_matrix/")
PDAC_PBMC_9 <- CreateSeuratObject(PDAC_PBMC_9, min.cells = 3)
PDAC_PBMC_9$SampleName <- "PDAC_PBMC_9"
PDAC_PBMC_9$Group <- "PDACI-Pri-PBMC"
PDAC_PBMC_9$ClinicalStage <- "I"
PDAC_PBMC_9$Class <- "Primary"
PDAC_PBMC_9$Tissue <- "PBMC"
PDAC_PBMC_9$Patient <- "PDAC"
PDAC_PBMC_9$Age <- 62
PDAC_PBMC_9$Gender <- "Male"
PDAC_PBMC_9$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_9, file="PDAC_PBMC_9.rds")

PDAC_PBMC_10A <- Read10X("data/PDAC_PBMC_10A/filtered_feature_bc_matrix/")
PDAC_PBMC_10A <- CreateSeuratObject(PDAC_PBMC_10A, min.cells = 3)
PDAC_PBMC_10A$SampleName <- "PDAC_PBMC_10"
PDAC_PBMC_10A$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_10A$Patient <- "PDAC"
PDAC_PBMC_10A$Age <- 84
PDAC_PBMC_10A$Gender <- "Female"
saveRDS(PDAC_PBMC_10A, file="PDAC_PBMC_10A.rds")

PDAC_PBMC_10B <- Read10X("data/PDAC_PBMC_10B/filtered_feature_bc_matrix/")
PDAC_PBMC_10B <- CreateSeuratObject(PDAC_PBMC_10B, min.cells = 3)
PDAC_PBMC_10B$SampleName <- "PDAC_PBMC_10"
PDAC_PBMC_10B$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_10B$Patient <- "PDAC"
PDAC_PBMC_10B$Age <- 84
PDAC_PBMC_10B$Gender <- "Female"
saveRDS(PDAC_PBMC_10B, file="PDAC_PBMC_10B.rds")

PDAC_PBMC_10 <- merge(PDAC_PBMC_10A, PDAC_PBMC_10B)
PDAC_PBMC_10$GEO <- "GSE155698"
PDAC_PBMC_10$ClinicalStage <- "II"
PDAC_PBMC_10$Class <- "Primary"
PDAC_PBMC_10$Tissue <- "PBMC"
saveRDS(PDAC_PBMC_10, file="PDAC_PBMC_10.rds")

PDAC_PBMC_11 <- Read10X("data/PDAC_PBMC_11/filtered_feature_bc_matrix/")
PDAC_PBMC_11 <- CreateSeuratObject(PDAC_PBMC_11, min.cells = 3)
PDAC_PBMC_11$SampleName <- "PDAC_PBMC_11"
PDAC_PBMC_11$Group <- "PDACIII-Pri-PBMC"
PDAC_PBMC_11$ClinicalStage <- "III"
PDAC_PBMC_11$Class <- "Primary"
PDAC_PBMC_11$Tissue <- "PBMC"
PDAC_PBMC_11$Patient <- "PDAC"
PDAC_PBMC_11$Age <- 60
PDAC_PBMC_11$Gender <- "Female"
PDAC_PBMC_11$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_11, file="PDAC_PBMC_11.rds")

PDAC_PBMC_12 <- Read10X("data/PDAC_PBMC_12/filtered_feature_bc_matrix/")
PDAC_PBMC_12 <- CreateSeuratObject(PDAC_PBMC_12, min.cells = 3)
PDAC_PBMC_12$SampleName <- "PDAC_PBMC_12"
PDAC_PBMC_12$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_12$ClinicalStage <- "II"
PDAC_PBMC_12$Class <- "Primary"
PDAC_PBMC_12$Tissue <- "PBMC"
PDAC_PBMC_12$Patient <- "PDAC"
PDAC_PBMC_12$Age <- 47
PDAC_PBMC_12$Gender <- "Female"
PDAC_PBMC_12$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_12, file="PDAC_PBMC_12.rds")

PDAC_PBMC_13 <- Read10X("data/PDAC_PBMC_13/filtered_feature_bc_matrix/")
PDAC_PBMC_13 <- CreateSeuratObject(PDAC_PBMC_13, min.cells = 3)
PDAC_PBMC_13$SampleName <- "PDAC_PBMC_13"
PDAC_PBMC_13$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_13$ClinicalStage <- "II"
PDAC_PBMC_13$Class <- "Primary"
PDAC_PBMC_13$Tissue <- "PBMC"
PDAC_PBMC_13$Patient <- "PDAC"
PDAC_PBMC_13$Age <- 70
PDAC_PBMC_13$Gender <- "Female"
PDAC_PBMC_13$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_13, file="PDAC_PBMC_13.rds")

PDAC_PBMC_14 <- Read10X("data/PDAC_PBMC_14/filtered_feature_bc_matrix/")
PDAC_PBMC_14 <- CreateSeuratObject(PDAC_PBMC_14, min.cells = 3)
PDAC_PBMC_14$SampleName <- "PDAC_PBMC_14"
PDAC_PBMC_14$Group <- "PDACIV-Met-PBMC"
PDAC_PBMC_14$ClinicalStage <- "IV"
PDAC_PBMC_14$Class <- "Metastatic"
PDAC_PBMC_14$Tissue <- "PBMC"
PDAC_PBMC_14$Patient <- "PDAC"
PDAC_PBMC_14$Age <- 74
PDAC_PBMC_14$Gender <- "Female"
PDAC_PBMC_14$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_14, file="PDAC_PBMC_14.rds")

PDAC_PBMC_15 <- Read10X("data/PDAC_PBMC_15/filtered_feature_bc_matrix/")
PDAC_PBMC_15 <- CreateSeuratObject(PDAC_PBMC_15, min.cells = 3)
PDAC_PBMC_15$SampleName <- "PDAC_PBMC_15"
PDAC_PBMC_15$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_15$ClinicalStage <- "II"
PDAC_PBMC_15$Class <- "Primary"
PDAC_PBMC_15$Tissue <- "PBMC"
PDAC_PBMC_15$Patient <- "PDAC"
PDAC_PBMC_15$Age <- 60
PDAC_PBMC_15$Gender <- "Female"
PDAC_PBMC_15$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_15, file="PDAC_PBMC_15.rds")

PDAC_PBMC_16 <- Read10X("data/PDAC_PBMC_16/filtered_feature_bc_matrix/")
PDAC_PBMC_16 <- CreateSeuratObject(PDAC_PBMC_16, min.cells = 3)
PDAC_PBMC_16$SampleName <- "PDAC_PBMC_16"
PDAC_PBMC_16$Group <- "PDACII-Pri-PBMC"
PDAC_PBMC_16$ClinicalStage <- "II"
PDAC_PBMC_16$Class <- "Primary"
PDAC_PBMC_16$Tissue <- "PBMC"
PDAC_PBMC_16$Patient <- "PDAC"
PDAC_PBMC_16$Age <- 74
PDAC_PBMC_16$Gender <- "Male"
PDAC_PBMC_16$GEO <- "GSE155698"
saveRDS(PDAC_PBMC_16, file="PDAC_PBMC_16.rds")

PDAC_TISSUE_1 <- Read10X("data/PDAC_TISSUE_1/filtered_feature_bc_matrix/")
PDAC_TISSUE_1 <- CreateSeuratObject(PDAC_TISSUE_1, min.cells = 3)
PDAC_TISSUE_1$SampleName <- "PDAC_TISSUE_1"
PDAC_TISSUE_1$Group <- "PDACIV-Metastasis"
PDAC_TISSUE_1$ClinicalStage <- "IV"
PDAC_TISSUE_1$Class <- "Metastatic"
PDAC_TISSUE_1$Tissue <- "Tissue"
PDAC_TISSUE_1$Patient <- "PDAC"
PDAC_TISSUE_1$Age <- 71
PDAC_TISSUE_1$Gender <- "Female"
PDAC_TISSUE_1$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_1, file="PDAC_TISSUE_1.rds")

PDAC_TISSUE_2 <- Read10X("data/PDAC_TISSUE_2/filtered_feature_bc_matrix/")
PDAC_TISSUE_2 <- CreateSeuratObject(PDAC_TISSUE_2, min.cells = 3)
PDAC_TISSUE_2$SampleName <- "PDAC_TISSUE_2"
PDAC_TISSUE_2$Group <- "PDACIV-Metastasis"
PDAC_TISSUE_2$ClinicalStage <- "IV"
PDAC_TISSUE_2$Class <- "Metastatic"
PDAC_TISSUE_2$Tissue <- "Tissue"
PDAC_TISSUE_2$Patient <- "PDAC"
PDAC_TISSUE_2$Age <- 67
PDAC_TISSUE_2$Gender <- "Male"
PDAC_TISSUE_2$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_2, file="PDAC_TISSUE_2.rds")

PDAC_TISSUE_3 <- Read10X("data/PDAC_TISSUE_3/filtered_feature_bc_matrix/")
PDAC_TISSUE_3 <- CreateSeuratObject(PDAC_TISSUE_3, min.cells = 3)
PDAC_TISSUE_3$SampleName <- "PDAC_TISSUE_3"
PDAC_TISSUE_3$Group <- "PDACIV-Metastasis"
PDAC_TISSUE_3$ClinicalStage <- "IV"
PDAC_TISSUE_3$Class <- "Metastatic"
PDAC_TISSUE_3$Tissue <- "Tissue"
PDAC_TISSUE_3$Patient <- "PDAC"
PDAC_TISSUE_3$Age <- 53
PDAC_TISSUE_3$Gender <- "Male"
PDAC_TISSUE_3$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_3, file="PDAC_TISSUE_3.rds")

PDAC_TISSUE_4 <- Read10X("data/PDAC_TISSUE_4/filtered_feature_bc_matrix/")
PDAC_TISSUE_4 <- CreateSeuratObject(PDAC_TISSUE_4, min.cells = 3)
PDAC_TISSUE_4$SampleName <- "PDAC_TISSUE_4"
PDAC_TISSUE_4$Group <- "PDACIII-Primary"
PDAC_TISSUE_4$ClinicalStage <- "III"
PDAC_TISSUE_4$Class <- "Primary"
PDAC_TISSUE_4$Tissue <- "Tissue"
PDAC_TISSUE_4$Patient <- "PDAC"
PDAC_TISSUE_4$Age <- 62
PDAC_TISSUE_4$Gender <- "Female"
PDAC_TISSUE_4$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_4, file="PDAC_TISSUE_4.rds")

PDAC_TISSUE_5 <- Read10X("data/PDAC_TISSUE_5/filtered_feature_bc_matrix/")
PDAC_TISSUE_5 <- CreateSeuratObject(PDAC_TISSUE_5, min.cells = 3)
PDAC_TISSUE_5$SampleName <- "PDAC_TISSUE_5"
PDAC_TISSUE_5$Group <- "PDACIII-Primary"
PDAC_TISSUE_5$ClinicalStage <- "III"
PDAC_TISSUE_5$Class <- "Primary"
PDAC_TISSUE_5$Tissue <- "Tissue"
PDAC_TISSUE_5$Patient <- "PDAC"
PDAC_TISSUE_5$Age <- 66
PDAC_TISSUE_5$Gender <- "Male"
PDAC_TISSUE_5$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_5, file="PDAC_TISSUE_5.rds")

PDAC_TISSUE_6 <- Read10X("data/PDAC_TISSUE_6/filtered_feature_bc_matrix/")
PDAC_TISSUE_6 <- CreateSeuratObject(PDAC_TISSUE_6, min.cells = 3)
PDAC_TISSUE_6$SampleName <- "PDAC_TISSUE_6"
PDAC_TISSUE_6$Group <- "PDACIII-Primary"
PDAC_TISSUE_6$ClinicalStage <- "III"
PDAC_TISSUE_6$Class <- "Primary"
PDAC_TISSUE_6$Tissue <- "Tissue"
PDAC_TISSUE_6$Patient <- "PDAC"
PDAC_TISSUE_6$Age <- 63
PDAC_TISSUE_6$Gender <- "Male"
PDAC_TISSUE_6$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_6, file="PDAC_TISSUE_6.rds")

PDAC_TISSUE_7 <- Read10X("data/PDAC_TISSUE_7/filtered_feature_bc_matrix/")
PDAC_TISSUE_7 <- CreateSeuratObject(PDAC_TISSUE_7, min.cells = 3)
PDAC_TISSUE_7$SampleName <- "PDAC_TISSUE_7"
PDAC_TISSUE_7$Group <- "PDACII-Primary"
PDAC_TISSUE_7$ClinicalStage <- "II"
PDAC_TISSUE_7$Class <- "Primary"
PDAC_TISSUE_7$Tissue <- "Tissue"
PDAC_TISSUE_7$Patient <- "PDAC"
PDAC_TISSUE_7$Age <- 49
PDAC_TISSUE_7$Gender <- "Female"
PDAC_TISSUE_7$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_7, file="PDAC_TISSUE_7.rds")

PDAC_TISSUE_8 <- Read10X("data/PDAC_TISSUE_8/filtered_feature_bc_matrix/")
PDAC_TISSUE_8 <- CreateSeuratObject(PDAC_TISSUE_8, min.cells = 3)
PDAC_TISSUE_8$SampleName <- "PDAC_TISSUE_8"
PDAC_TISSUE_8$Group <- "PDACII-Primary"
PDAC_TISSUE_8$ClinicalStage <- "II"
PDAC_TISSUE_8$Class <- "Primary"
PDAC_TISSUE_8$Tissue <- "Tissue"
PDAC_TISSUE_8$Patient <- "PDAC"
PDAC_TISSUE_8$Age <- 69
PDAC_TISSUE_8$Gender <- "Female"
PDAC_TISSUE_8$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_8, file="PDAC_TISSUE_8.rds")

PDAC_TISSUE_9 <- Read10X("data/PDAC_TISSUE_9/filtered_feature_bc_matrix/")
PDAC_TISSUE_9 <- CreateSeuratObject(PDAC_TISSUE_9, min.cells = 3)
PDAC_TISSUE_9$SampleName <- "PDAC_TISSUE_9"
PDAC_TISSUE_9$Group <- "PDACI-Primary"
PDAC_TISSUE_9$ClinicalStage <- "I"
PDAC_TISSUE_9$Class <- "Primary"
PDAC_TISSUE_9$Tissue <- "Tissue"
PDAC_TISSUE_9$Patient <- "PDAC"
PDAC_TISSUE_9$Age <- 69
PDAC_TISSUE_9$Gender <- "Male"
PDAC_TISSUE_9$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_9, file="PDAC_TISSUE_9.rds")

PDAC_TISSUE_10 <- Read10X("data/PDAC_TISSUE_10/filtered_feature_bc_matrix/")
PDAC_TISSUE_10 <- CreateSeuratObject(PDAC_TISSUE_10, min.cells = 3)
PDAC_TISSUE_10$SampleName <- "PDAC_TISSUE_10"
PDAC_TISSUE_10$Group <- "PDACII-Primary"
PDAC_TISSUE_10$ClinicalStage <- "II"
PDAC_TISSUE_10$Class <- "Primary"
PDAC_TISSUE_10$Tissue <- "Tissue"
PDAC_TISSUE_10$Patient <- "PDAC"
PDAC_TISSUE_10$Age <- 84
PDAC_TISSUE_10$Gender <- "Female"
PDAC_TISSUE_10$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_10, file="PDAC_TISSUE_10.rds")

PDAC_TISSUE_11A <- Read10X("data/PDAC_TISSUE_11A/filtered_feature_bc_matrix/")
PDAC_TISSUE_11A <- CreateSeuratObject(PDAC_TISSUE_11A, min.cells = 3)
PDAC_TISSUE_11A$SampleName <- "PDAC_TISSUE_11"
PDAC_TISSUE_11A$Group <- "PDACIII-Primary"
PDAC_TISSUE_11A$Patient <- "PDAC"
PDAC_TISSUE_11A$Age <- 60
PDAC_TISSUE_11A$Gender <- "Female"
saveRDS(PDAC_TISSUE_11A, file="PDAC_TISSUE_11A.rds")

PDAC_TISSUE_11B <- Read10X("data/PDAC_TISSUE_11B/filtered_feature_bc_matrix/")
PDAC_TISSUE_11B <- CreateSeuratObject(PDAC_TISSUE_11B, min.cells = 3)
PDAC_TISSUE_11B$SampleName <- "PDAC_TISSUE_11"
PDAC_TISSUE_11B$Group <- "PDACIII-Primary"
PDAC_TISSUE_11B$Patient <- "PDAC"
PDAC_TISSUE_11B$Age <- 60
PDAC_TISSUE_11B$Gender <- "Female"
saveRDS(PDAC_TISSUE_11B, file="PDAC_TISSUE_11B.rds")

PDAC_TISSUE_11 <- merge(PDAC_TISSUE_11A,PDAC_TISSUE_11B)
PDAC_TISSUE_11$GEO <- "GSE155698"
PDAC_TISSUE_11$ClinicalStage <- "III"
PDAC_TISSUE_11$Class <- "Primary"
PDAC_TISSUE_11$Tissue <- "Tissue"
saveRDS(PDAC_TISSUE_11, file="PDAC_TISSUE_11.rds")

PDAC_TISSUE_12 <- Read10X("data/PDAC_TISSUE_12/filtered_feature_bc_matrix/")
PDAC_TISSUE_12 <- CreateSeuratObject(PDAC_TISSUE_12, min.cells = 3)
PDAC_TISSUE_12$SampleName <- "PDAC_TISSUE_12"
PDAC_TISSUE_12$Group <- "PDACII-Primary"
PDAC_TISSUE_12$ClinicalStage <- "II"
PDAC_TISSUE_12$Class <- "Primary"
PDAC_TISSUE_12$Tissue <- "Tissue"
PDAC_TISSUE_12$Patient <- "PDAC"
PDAC_TISSUE_12$Age <- 47
PDAC_TISSUE_12$Gender <- "Female"
PDAC_TISSUE_12$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_12, file="PDAC_TISSUE_12.rds")

PDAC_TISSUE_13 <- Read10X("data/PDAC_TISSUE_13/filtered_feature_bc_matrix/")
PDAC_TISSUE_13 <- CreateSeuratObject(PDAC_TISSUE_13, min.cells = 3)
PDAC_TISSUE_13$SampleName <- "PDAC_TISSUE_13"
PDAC_TISSUE_13$Group <- "PDACII-Primary"
PDAC_TISSUE_13$ClinicalStage <- "II"
PDAC_TISSUE_13$Class <- "Primary"
PDAC_TISSUE_13$Tissue <- "Tissue"
PDAC_TISSUE_13$Patient <- "PDAC"
PDAC_TISSUE_13$Age <- 70
PDAC_TISSUE_13$Gender <- "Female"
PDAC_TISSUE_13$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_13, file="PDAC_TISSUE_13.rds")

PDAC_TISSUE_14 <- Read10X_h5("data/PDAC_TISSUE_14/filtered_feature_bc_matrix.h5")
PDAC_TISSUE_14 <- CreateSeuratObject(PDAC_TISSUE_14, min.cells = 3)
PDAC_TISSUE_14$SampleName <- "PDAC_TISSUE_14"
PDAC_TISSUE_14$Group <- "PDACIV-Metastasis"
PDAC_TISSUE_14$ClinicalStage <- "IV"
PDAC_TISSUE_14$Class <- "Metastatic"
PDAC_TISSUE_14$Tissue <- "Tissue"
PDAC_TISSUE_14$Patient <- "PDAC"
PDAC_TISSUE_14$Age <- 74
PDAC_TISSUE_14$Gender <- "Female"
PDAC_TISSUE_14$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_14, file="PDAC_TISSUE_14.rds")

PDAC_TISSUE_15 <- Read10X("data/PDAC_TISSUE_15/filtered_feature_bc_matrix/")
PDAC_TISSUE_15 <- CreateSeuratObject(PDAC_TISSUE_15, min.cells = 3)
PDAC_TISSUE_15$SampleName <- "PDAC_TISSUE_15"
PDAC_TISSUE_15$Group <- "PDACII-Primary"
PDAC_TISSUE_15$ClinicalStage <- "II"
PDAC_TISSUE_15$Class <- "Primary"
PDAC_TISSUE_15$Tissue <- "Tissue"
PDAC_TISSUE_15$Patient <- "PDAC"
PDAC_TISSUE_15$Age <- 60
PDAC_TISSUE_15$Gender <- "Male"
PDAC_TISSUE_15$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_15, file="PDAC_TISSUE_15.rds")

PDAC_TISSUE_16 <- Read10X("data/PDAC_TISSUE_16/filtered_feature_bc_matrix/")
PDAC_TISSUE_16 <- CreateSeuratObject(PDAC_TISSUE_16, min.cells = 3)
PDAC_TISSUE_16$SampleName <- "PDAC_TISSUE_16"
PDAC_TISSUE_16$Group <- "PDACII-Primary"
PDAC_TISSUE_16$ClinicalStage <- "II"
PDAC_TISSUE_16$Class <- "Primary"
PDAC_TISSUE_16$Tissue <- "Tissue"
PDAC_TISSUE_16$Patient <- "PDAC"
PDAC_TISSUE_16$Age <- 74
PDAC_TISSUE_16$Gender <- "Male"
PDAC_TISSUE_16$GEO <- "GSE155698"
saveRDS(PDAC_TISSUE_16, file="PDAC_TISSUE_16.rds")

library(future)
library(future.apply)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 200000 * 1024^2)

seurat.list <- list(
  AdjNorm_TISSUE_1 = AdjNorm_TISSUE_1,
  AdjNorm_TISSUE_2 = AdjNorm_TISSUE_2,
  AdjNorm_TISSUE_3 = AdjNorm_TISSUE_3,
  Healthy_PBMC_1 = Healthy_PBMC_1,
  Healthy_PBMC_2 = Healthy_PBMC_2,
  Healthy_PBMC_3 = Healthy_PBMC_3,
  Healthy_PBMC_4 = Healthy_PBMC_4,
  PDAC_PBMC_1 = PDAC_PBMC_1,
  PDAC_PBMC_2 = PDAC_PBMC_2,
  PDAC_PBMC_3 = PDAC_PBMC_3,
  PDAC_PBMC_4 = PDAC_PBMC_4,
  PDAC_PBMC_5 = PDAC_PBMC_5,
  PDAC_PBMC_6 = PDAC_PBMC_6,
  PDAC_PBMC_7 = PDAC_PBMC_7,
  PDAC_PBMC_8 = PDAC_PBMC_8,
  PDAC_PBMC_9 = PDAC_PBMC_9,
  PDAC_PBMC_10 = PDAC_PBMC_10,
  PDAC_PBMC_11 = PDAC_PBMC_11,
  PDAC_PBMC_12 = PDAC_PBMC_12,
  PDAC_PBMC_13 = PDAC_PBMC_13,
  PDAC_PBMC_14 = PDAC_PBMC_14,
  PDAC_PBMC_15 = PDAC_PBMC_15,
  PDAC_PBMC_16 = PDAC_PBMC_16,
  PDAC_TISSUE_1 = PDAC_TISSUE_1,
  PDAC_TISSUE_2 = PDAC_TISSUE_2,
  PDAC_TISSUE_3 = PDAC_TISSUE_3,
  PDAC_TISSUE_4 = PDAC_TISSUE_4,
  PDAC_TISSUE_5 = PDAC_TISSUE_5,
  PDAC_TISSUE_6 = PDAC_TISSUE_6,
  PDAC_TISSUE_7 = PDAC_TISSUE_7,
  PDAC_TISSUE_8 = PDAC_TISSUE_8,
  PDAC_TISSUE_9 = PDAC_TISSUE_9,
  PDAC_TISSUE_10 = PDAC_TISSUE_10,
  PDAC_TISSUE_11 = PDAC_TISSUE_11,
  PDAC_TISSUE_12 = PDAC_TISSUE_12,
  PDAC_TISSUE_13 = PDAC_TISSUE_13,
  PDAC_TISSUE_14 = PDAC_TISSUE_14,
  PDAC_TISSUE_15 = PDAC_TISSUE_15,
  PDAC_TISSUE_16 = PDAC_TISSUE_16
)

seurat.list <- lapply(seurat.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method="vst",nfeatures=2000)
})
saveRDS(seurat.list, file = "seurat.list.rds")

total.tissue.list <- list(
  AdjNorm_TISSUE_1 = seurat.list$AdjNorm_TISSUE_1,
  AdjNorm_TISSUE_2 = seurat.list$AdjNorm_TISSUE_2,
  AdjNorm_TISSUE_3 = seurat.list$AdjNorm_TISSUE_3,
  PDAC_TISSUE_1 = seurat.list$PDAC_TISSUE_1,
  PDAC_TISSUE_2 = seurat.list$PDAC_TISSUE_2,
  PDAC_TISSUE_3 = seurat.list$PDAC_TISSUE_3,
  PDAC_TISSUE_4 = seurat.list$PDAC_TISSUE_4,
  PDAC_TISSUE_5 = seurat.list$PDAC_TISSUE_5,
  PDAC_TISSUE_6 = seurat.list$PDAC_TISSUE_6,
  PDAC_TISSUE_7 = seurat.list$PDAC_TISSUE_7,
  PDAC_TISSUE_8 = seurat.list$PDAC_TISSUE_8,
  PDAC_TISSUE_9 = seurat.list$PDAC_TISSUE_9,
  PDAC_TISSUE_10 = seurat.list$PDAC_TISSUE_10,
  PDAC_TISSUE_11 = seurat.list$PDAC_TISSUE_11,
  PDAC_TISSUE_12 = seurat.list$PDAC_TISSUE_12,
  PDAC_TISSUE_13 = seurat.list$PDAC_TISSUE_13,
  PDAC_TISSUE_14 = seurat.list$PDAC_TISSUE_14,
  PDAC_TISSUE_15 = seurat.list$PDAC_TISSUE_15,
  PDAC_TISSUE_16 = seurat.list$PDAC_TISSUE_16
)

total.pbmc.list <- list(
  Healthy_PBMC_1 = seurat.list$Healthy_PBMC_1,
  Healthy_PBMC_2 = seurat.list$Healthy_PBMC_2,
  Healthy_PBMC_3 = seurat.list$Healthy_PBMC_3,
  Healthy_PBMC_4 = seurat.list$Healthy_PBMC_4,
  PDAC_PBMC_1 = seurat.list$PDAC_PBMC_1,
  PDAC_PBMC_2 = seurat.list$PDAC_PBMC_2,
  PDAC_PBMC_3 = seurat.list$PDAC_PBMC_3,
  PDAC_PBMC_4 = seurat.list$PDAC_PBMC_4,
  PDAC_PBMC_5 = seurat.list$PDAC_PBMC_5,
  PDAC_PBMC_6 = seurat.list$PDAC_PBMC_6,
  PDAC_PBMC_7 = seurat.list$PDAC_PBMC_7,
  PDAC_PBMC_8 = seurat.list$PDAC_PBMC_8,
  PDAC_PBMC_9 = seurat.list$PDAC_PBMC_9,
  PDAC_PBMC_10 = seurat.list$PDAC_PBMC_10,
  PDAC_PBMC_11 = seurat.list$PDAC_PBMC_11,
  PDAC_PBMC_12 = seurat.list$PDAC_PBMC_12,
  PDAC_PBMC_13 = seurat.list$PDAC_PBMC_13,
  PDAC_PBMC_14 = seurat.list$PDAC_PBMC_14,
  PDAC_PBMC_15 = seurat.list$PDAC_PBMC_15,
  PDAC_PBMC_16 = seurat.list$PDAC_PBMC_16
)

saveRDS(total.tissue.list, file = "total.tissue.list.rds")
saveRDS(total.pbmc.list,file = "total.pbmc.list.rds")
rm(seurat.list);gc()

total.tissue.merged <- merge(total.tissue.list[[1]], total.tissue.list[2:length(total.tissue.list)])
total.pbmc.merged <- merge(total.pbmc.list[[1]], total.pbmc.list[2:length(total.pbmc.list)])
saveRDS(total.tissue.merged,file = "total.tissue.merged.rds")
saveRDS(total.pbmc.merged, file = "total.pbmc.merged.rds")

# total.merged <- merge(total.tissue.merged, total.pbmc.merged)
# saveRDS(total.merged, file = "total.merged.rds")

seurat.list <- list(total.tissue.merged = total.tissue.merged, total.pbmc.merged = total.pbmc.merged)

scRNA.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
scRNA.anchors <- FindIntegrationAnchors(object.list = seurat.list,anchor.features = scRNA.features, verbose = F)
scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, verbose = F)
saveRDS(scRNA.anchors, "anchors.rds")
saveRDS(scRNA.integrated, "total.integrated.rds")

# --------------------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

if(!file.exists("figures")) dir.create("figures")

# preprocessing
total.integrated <- readRDS("total.integrated.rds")

DefaultAssay(total.integrated) <- "RNA"
total.integrated <- NormalizeData(total.integrated)

#Percent Mitochondrial Genes
#QC Metric used to remove cells with overabundant Mitochondrial genes, typically associated with nuclear wash out during sequencing
total.integrated[["percent.mt"]] <- PercentageFeatureSet(total.integrated, pattern = "^MT-")
p1 <- VlnPlot(total.integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 0)
ggsave(paste0("figures/","metrics.beforeQC.pdf"),p1,width = 8,height = 6)

## Quality Control
# QC criteria: 
#     i). remove cells had either lower than 200 or higher than 5000 expressed genes
#    ii). discarded cells with more than 30,000 UMIs
#   iii). mitochondria content higher than 30%
total.integrated <- subset(total.integrated, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 30000 & percent.mt < 30)
p2 <- VlnPlot(total.integrated, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3, pt.size = 0)
ggsave(paste0("figures/","metrics.afterQC.pdf"),p1,width = 8,height = 6)

# Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
total.integrated <- CellCycleScoring(total.integrated,s.features = s.genes, g2m.features = g2m.features, set.ident = T)
total.integrated$CC.Difference <- total.integrated$S.score - total.integrated$G2M.Score

# Switch assay to integrated
DefaultAssay(total.integrated) <- "integrated"

# scale data
total.integrated <- ScaleData(object = total.integrated, vars.to.regress = "CC.Difference", features = rownames(total.integrated))

## Dimensional reduction
# Run PCA and Determine Dimensions for 90% Variance
total.integrated <- RunPCA(object = total.integrated, features = VariableFeatures(total.integrated))

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

PC.num <- PCDeterminators(total.integrated)

# Find Neighbors + Find Clusters
total.integrated <- FindNeighbors(total.integrated, dims = 1:PC.num)
total.integrated <- FindClusters(total.integrated, resolution = 1.2)
# Run UMAP and get unlabelled cluster UMAP & tSNE and violin plot (without harmony batch correction)
total.integrated <- RunUMAP(total.integrated, dims = 1:PC.num)
total.integrated <- RunTSNE(total.integrated,dims = 1:PC.num)

mypal <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
           "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
           "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
           "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
           "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
           "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
           "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
           "#e598d7", "#343b20", "#bbb5d9", "#975223")
mypal2 <- c("#c14089","#6f5553","#E5C494","#738f4c","#bb6240","#66C2A5","#2dfd29","#0c0fdc")

# UMAP
p1 <- DimPlot(total.integrated, reduction = "umap", label = T, pt.size = 0.2)
ggsave(paste0("figures/","umap.pdf"),p1,width = 7,height = 6)
# split by sampleName
p1.1.1 <- DimPlot(total.integrated, reduction = "umap", label = T, pt.size = 0.7, split.by = "SampleName", ncol = 5)
p1.1.2 <- DimPlot(total.integrated,reduction = "umap",label = F,pt.size = 0.2,group.by = "SampleName",cols = mypal[1:39])
ggsave(paste0("figures/","umap.splitBy.SampleName.pdf"),p1.1.1,width = 21,height = 30)
ggsave(paste0("figures/","umap.groupBy.SampleName.pdf"),p1.1.2,width = 9,height = 5)
# split/group by clinical stage
p1.2.1 <- DimPlot(total.integrated,reduction = "umap", label = T, pt.size = 0.2, split.by = "ClinicalStage", ncol = 3)
p1.2.2 <- DimPlot(total.integrated,reduction = "umap", label = T, pt.size = 0.2, group.by = "ClinicalStage", cols = mypal2[1:5])
ggsave(paste0("figures/","umap.splitBy.ClinicalStage.pdf"),p1.2.1, width = 12,height = 8)
ggsave(paste0("figures/","umap.groupBy.ClinicalStage.pdf"),p1.2.2, width = 7, height = 6)
# split by class (primary, metastatic, paracancerous)
p1.3.1 <- DimPlot(total.integrated,reduction = "umap",label = T,pt.size = 0.2, split.by = "Class", ncol = 2)
p1.3.2 <- DimPlot(total.integrated,reduction = "umap",label = T,pt.size = 0.2, group.by = "Class", mypal2[1:4])
ggsave(paste0("figures/","umap.splitBy.Class.pdf"),p1.3.1, width = 9,height = 8)
ggsave(paste0("figures/","umap.groupBy.Class.pdf"),p1.3.2, width = 7,height = 6)
# split by tissue (tissue vs. pbmc)
p1.4.1 <- DimPlot(total.integrated,reduction = "umap",label = T,pt.size = 0.2,split.by = "Tissue",ncol = 2)
p1.4.2 <- DimPlot(total.integrated,reduction = "umap",label = T,pt.size = 0.2, group.by = "Tissue", mypal2[1:4])
ggsave(paste0("figures/","umap.splitBy.Tissue.pdf"),p1.4.1,width = 9,height = 4)
ggsave(paste0("figures/","umap.groupBy.Tissue.pdf"),p1.4.2, width = 7, height = 6)
# split by patient (PDAC vs Normal)
p1.5.1 <- DimPlot(total.integrated, reduction = "umap", label = T, pt.size = 0.2, split.by = "Patient", ncol = 2)
p1.5.2 <- DimPlot(total.integrated, reduction = "umap", label = T, pt.size = 0.2, group.by = "Patient", mypal2[1:4])
ggsave(paste0("figures/","umap.splitBy.Patient.pdf"),p1.5.1,width = 9,height = 4)
ggsave(paste0("figures/","umap.groupBy.Patient.pdf"),p1.5.2,width = 7,height = 6)

# tSNE
p2 <- DimPlot(total.integrated, reduction = "tsne", label = T, pt.size = 0.2)
ggsave(paste0("figures/","tsne.pdf"),p2,width = 7,height = 6)
# split by sampleName
p2.1.1 <- DimPlot(total.integrated, reduction = "tsne", label = T, pt.size = 0.7, split.by = "SampleName", ncol = 5)
p2.1.2 <- DimPlot(total.integrated,reduction = "tsne",label = F,pt.size = 0.2,group.by = "SampleName",cols = mypal[1:39])
ggsave(paste0("figures/","tsne.splitBy.SampleName.pdf"),p2.1.1,width = 21,height = 30)
ggsave(paste0("figures/","tsne.groupBy.SampleName.pdf"),p2.1.2,width = 9,height = 5)
# split/group by clinical stage
p2.2.1 <- DimPlot(total.integrated,reduction = "tsne", label = T, pt.size = 0.2, split.by = "ClinicalStage", ncol = 3)
p2.2.2 <- DimPlot(total.integrated,reduction = "tsne", label = T, pt.size = 0.2, group.by = "ClinicalStage", cols = mypal2[1:5])
ggsave(paste0("figures/","tsne.splitBy.ClinicalStage.pdf"),p2.2.1, width = 12,height = 8)
ggsave(paste0("figures/","tsne.groupBy.ClinicalStage.pdf"),p2.2.2, width = 7, height = 6)
# split by class (primary, metastatic, paracancerous)
p2.3.1 <- DimPlot(total.integrated,reduction = "tsne",label = T,pt.size = 0.2, split.by = "Class", ncol = 2)
p2.3.2 <- DimPlot(total.integrated,reduction = "tsne",label = T,pt.size = 0.2, group.by = "Class", cols = mypal2[1:4])
ggsave(paste0("figures/","tsne.splitBy.Class.pdf"),p2.3.1, width = 9,height = 8)
ggsave(paste0("figures/","tsne.groupBy.Class.pdf"),p2.3.2, width = 7,height = 6)
# split by tissue (tissue vs. pbmc)
p2.4.1 <- DimPlot(total.integrated,reduction = "tsne",label = T,pt.size = 0.2,split.by = "Tissue",ncol = 2)
p2.4.2 <- DimPlot(total.integrated,reduction = "tsne",label = T,pt.size = 0.2, group.by = "Tissue", cols = mypal2[1:2])
ggsave(paste0("figures/","tsne.splitBy.Tissue.pdf"),p2.4.1,width = 9,height = 4)
ggsave(paste0("figures/","tsne.groupBy.Tissue.pdf"),p2.4.2, width = 7, height = 6)
# split by patient (PDAC vs Normal)
p2.5.1 <- DimPlot(total.integrated, reduction = "tsne", label = T, pt.size = 0.2, split.by = "Patient", ncol = 2)
p2.5.2 <- DimPlot(total.integrated, reduction = "tsne", label = T, pt.size = 0.2, group.by = "Patient", cols = mypal2[1:2])
ggsave(paste0("figures/","tsne.splitBy.Patient.pdf"),p2.5.1,width = 9,height = 4)
ggsave(paste0("figures/","tsne.groupBy.Patient.pdf"),p2.5.2,width = 7,height = 6)

saveRDS(total.integrated, file = "total.integrated.rds")

# Find neighbors and clusters WITH harmony batch correction
Idents(total.integrated) <- total.integrated$SampleName
total.integrated <- RunHarmony(total.integrated, group.by.vars = "SampleName")
total.integrated.harmony <- FindNeighbors(total.integrated,dims = 1:PC.num, reduction = "harmony") %>%
  FindClusters(resolution = 1.2, reduction = "harmony")

saveRDS(total.integrated.harmony, file = "input/total.integrated.harmony.rds")

# raw vs. harmony batch corrected
# samplename
p.raw.sn <- DimPlot(total.integrated, reduction = "tsne",label = F, pt.size = 0.1, group.by = "SampleName", cols = mypal[1:39])
p.harmony.sn <- DimPlot(total.integrated.harmony, reduction = "tsne",label = F, pt.size = 0.2, group.by = "SampleName", cols = mypal[1:39])
ggsave(paste0("figures/","tsne.raw.groupBy.SampleName"),p.raw.sn, width = 7,height = 5)
ggsave(paste0("figures/","tsne.harmony.corrected.groupBy.SampleName"),p.harmony.sn,width = 7,height = 5)
# seurat_clusters
p.raw.cls <- DimPlot(total.integrated,reduction = "tsne",label = T,pt.size = 0.1,group.by = "seurat_clusters",cols = mypal)
p.harmony.cls <- DimPlot(total.integrated.harmony,reduction = "tsne",label = T,pt.size = 0.1,group.by = "seurat_clusters",cols = mypal[1:40])
ggsave(paste0("figures/","tsne.raw.groupBy.clusters"),p.raw.cls, width = 7,height = 6)
ggsave(paste0("figures/","tsne.harmony.corrected.groupBy.clusters"),p.harmony.cls,width = 7,height = 6)

### Stacked Barplot (cell proportion)
# SampleName vs. seurat_clusters
p.sn <- ggplot(total.integrated.harmony@meta.data, aes(x = SampleName, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal[1:40]) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.direction = "vertical") +
  labs(x = "Samples", y = "Cluster Fraction") +
  guides(fill = guide_legend(ncol = 1))
ggsave(paste0("figures/","cluster_prop.pdf"),p.sn,width = 8,height = 12)

# Patient vs. seurat_clusters
p.pt <- ggplot(total.integrated.harmony@meta.data, aes(x = Patient, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal[1:40]) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Patient", y = "Cluster Fraction")
ggsave(paste0("figures/","patient_prop.pdf"),p.pt, width = 14,height = 3)

# ClinicalStage vs. seurat_clusters
p.cs <- ggplot(total.integrated.harmony@meta.data, aes(x = ClinicalStage, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal[1:40]) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Patient", y = "Cluster Fraction")
ggsave(paste0("figures/","ClinicalStage_prop.pdf"),p.cs, width = 14,height = 3)

# Tissue vs. seurat_clusters
p.ts <- ggplot(total.integrated.harmony@meta.data, aes(x = Tissue, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal[1:40]) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Tissue", y = "Cluster Fraction")
ggsave(paste0("figures/","Tissue_prop.pdf"),p.ts, width = 14,height = 3)

# Tissue vs. seurat_clusters
p.cls <- ggplot(total.integrated.harmony@meta.data, aes(x = Class, fill = seurat_clusters)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal[1:40]) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Class", y = "Cluster Fraction")
ggsave(paste0("figures/","Class_prop.pdf"),p.cls, width = 14,height = 3)

# ---- FIND CLUSTER MARKER ----

# This funtion automates the FindMarkers function and uses the list of markers to broadly identify
# cell types on a preselected list of markers. Markers chosen for human samples.
# The output is a dataframe containing the FindMarkers output

all.markers <- FindAllMarkers(object = total.integrated.harmony, min.pct = 0.25, only.pos = TRUE)
save(all.markers,file = "output/markers/all.markers.rda")

## heatmap of top markers
library(Seurat)
library(dplyr)
library(ggplot2)
mypal <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
           "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
           "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
           "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
           "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
           "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
           "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
           "#e598d7", "#343b20", "#bbb5d9", "#975223")
load("output/markers/all.markers.rda")
total.integrated.harmony <- readRDS("input/total.integrated.harmony.rds")
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
library(ComplexHeatmap)
mat <- GetAssayData(total.integrated.harmony, assay = "RNA", slot = "data")
mat[mat > 5] <- 5
cluster.info <- sort(total.integrated.harmony$seurat_clusters)
mat <- as.matrix(mat[top10$gene, names(cluster.info)])
gene <- c("IL7R","CCR7","IL8","CXCR2","LTB","CRIP1","KRT19","KRT18","KRT8","KRT7","FCGR3B","S100A8","S100A9",
          "GZMK","CCL5","GNLY","GZMB","NKG7","APOE","SPP1","C1QA","C1QB","LYZ","VCAN","CCL5","MMP9","SPARC",
          "KRT18","EPCAM","MS4A1","CD79A","CD79B","PRSS1","CTRB1","TPSAB1","GATA2","FCGR3A","IGLL5","IGJ","MZB1",
          "STMN1","TUBA1B","IFIT3","ACTA2","RGS5","HLA-DRA","HLA-DPA1","HLA-DPB1","CLU","KRT18","KRT8","SPARCL1",
          "PPBP","IL32","IFITM1","PLD4","TCF4","ISG15","MX1","CPA2","GP2","SEC11C","KLRB1","CYP1B1","MMP9","MME",
          "APOC1","IGLL5","MMP7")
gene <- unique(gene)
gene.pos <- which(rownames(mat) %in% gene)
row.anno <- ComplexHeatmap::rowAnnotation(gene = ComplexHeatmap::anno_mark(at = gene.pos, labels = rownames(mat)[which(rownames(mat) %in% gene)]))
col <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
         "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
         "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
         "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
         "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
         "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
         "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
         "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b")
names(col) <- levels(cluster.info)
top.anno <- ComplexHeatmap::HeatmapAnnotation(
  cluster = ComplexHeatmap::anno_block(gp = gpar(fill=col),labels = levels(cluster.info),labels_gp = gpar(cex=0.5,col='white'))
)
col.fun <- circlize::colorRamp2(c(0,1,5),c('#377EB8','white','#E41A1C'))
pdf("diff.cluster.hetamap.pdf",width = 16,height = 24)
ComplexHeatmap::Heatmap(matrix = mat, cluster_rows = F, cluster_columns = F, show_column_names = F,
                        show_row_names = F, column_split = cluster.info, right_annotation = row.anno,
                        column_title = NULL, top_annotation = top.anno, 
                        heatmap_legend_param = list(title='Expression',title_position='leftcenter-rot'),
                        col = col.fun)
dev.off()
# pdf("diff.cluster.hetamap.pdf",width = 8,height = 24)
# DoHeatmap(total.integrated.harmony,features = top5$gene,group.by = "seurat_clusters", size = 3,group.colors = mypal) + 
#   NoLegend() + scale_fill_gradientn(colours = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100))
# dev.off()

# clusterMarkerTBL <- function(clust.num, seurat.obj, num.markers, pct){
#   clusterMarkers <- FindMarkers(object = seurat.obj, ident.1 = clust.num[1], min.pct = pct)
#   clusterTable <- as.data.frame(head(clusterMarkers, num.markers))
#   markerLableRow <- length(clusterTable$p_val) + 1
#   markerLabel <- c("Cluster", clust.num[1], "Markers", "Are", "Above")
#   clusterTable[markerLableRow,] <- markerLabel
#   
#   if(length(clust.num)>2){
#     for(i in clust.num[2:length(clust.num)]){
#       clusterMarkers <- FindMarkers(object = seurat.obj, ident.1 = i, min.pct = pct)
#       clusterTable <- rbind(clusterTable, head(clusterMarkers,num.markers))
#       markerLableRow <- length(clusterTable$p_val)+1
#       markerLabel <- c("Cluster",i,"Markers","Are","Above")
#       clusterTable[markerLableRow,] <- markerLabel
#     }
#   }
#   return(clusterTable)
# }
# # change cluster marker range to match desired cluster numbers from seurat_clusters
# total.integrated.harmony_cluster <- clusterMarkerTBL(clust.num = c(0:(length(levels(total.integrated.harmony$seurat_clusters))-1)),
#                                                      seurat.obj = total.integrated.harmony,
#                                                      num.markers = 100,
#                                                      pct = 0.25)
# write.csv(total.integrated.harmony_cluster,file=paste0("total.integrated.harmony.cluster0-",length(levels(total.integrated.harmony$seurat_clusters)),".csv"))

## ------------------------------------------ CellType -------------------------
rm(list = ls());gc()
total.integrated.harmony <- readRDS("rds/total.integrated.harmony.rds")

rownames(total.integrated.harmony)[grep("CD38",rownames(total.integrated.harmony))]

# markers
iCAF <- c("IL6","PDGFRA","CFD","PLA2G2A","HAS1","CXCL2","CCL2","CLU","EMP1","LMNA")
myCAF <- c("TAGLN","ACTA2","MMP11","PDGFRB","HOPX","POSTN")
ductal_cell <- c("AMBP","CFTR","MMP7","KRT19","KRT7","TSPAN8","SLPI")
Acinar_cell <- c("PRSS1","CTRB1","CTRB2","REG1B","SPINK1","AMY2A")
Endocrine_cell <- c("CHGB","CHGA","INS","IAPP")
Stellate_cell <- c("RGS5","ACTA2","PDGFRB","ADIRF")
Fibroblast_cell <- c("LUM","DCN","COL1A1","ACTA2","SPARC","CDH11","PDGFRA","PDGFRB","COL3A1",
                     "RGS5","IGFBP7","PDPN","MCAM","IL6","APOE","GLI1","GLI2","GLI3","PDGFA")
Endothelial_cell <- c("CDH5","PLVAP","VWF","CLDN5","KDR","PECAM1")
Mast_cell <- c("TPSAB1","CPA3")
Macrophage_cell <- c("AIF1","FCGR1A","CD14","CD68","CSF3")
T_cells <- c("CD3D","CD3E","CD4","CD8A","CD3G","IL7R","LEF1")
B_cell <- c("MS4A1","CD79A","CD79B","CD52","CD19","SDC1","IGJ","IGLL5","CXCR4","KIT","CD27","HLA-DRA")
Epithelial_cell <- c("EPCAM","ACTA2","KRT7","KRT8","KRT18","KRT19","CDH1","PRSS1","CTRB2",
                     "REG1A","CLU","MKI67","SPINK1","TFF1","MUC1")
Lymphocytes <- c("CD3D","IL7R","CD3G")
DCs <- c("FCER1A","CD1A","LYZ","CLEC9A","BATF3","IRF8","IDO1","CD207","CD1C","HLA-DRA","CCL22","LAMP3","IL22RA2","CD101")
EMT_like_cell <- c("CDH2","SNAI2","ZEB1")
Myeloid_cell <- c("CD14","ITGAM","MNDA","MPEG1","ITGAX","FCGR3A","FCGR3B","APOE","C1QA","MARCO","LYZ","HLA-DRA")
NK <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38")

# custom theme
mycols <- c("grey88","DarkCyan")
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank())

## ----------------- FeaturePlot ----------------------------------

DefaultAssay(total.integrated.harmony) <- "RNA"

# iCAF
p.feature.iCAF <- list()
iCAF.match <- match(iCAF,rownames(total.integrated.harmony@assays$RNA))
iCAF <- rownames(total.integrated.harmony@assays$RNA)[iCAF.match]
iCAF <- iCAF[!is.na(iCAF)]
for(i in seq_along(iCAF)){
  p.feature.iCAF[[i]] <- FeaturePlot(total.integrated.harmony,features = iCAF[i], cols = mycols, reduction = "tsne",pt.size = 0.05) + 
    NoLegend() + mytheme
}
p.feature.comb.iCAF <- cowplot::plot_grid(plotlist = p.feature.iCAF, ncol = 5, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.iCAF.pdf"),p.feature.comb.iCAF,width = 14,height = 6)
# myCAF
p.feature.myCAF <- list()
myCAF.match <- match(myCAF,rownames(total.integrated.harmony@assays$RNA))
myCAF <- rownames(total.integrated.harmony@assays$RNA)[myCAF.match]
myCAF <- myCAF[!is.na(myCAF)]
for(i in seq_along(myCAF)){
  p.feature.myCAF[[i]] <- FeaturePlot(total.integrated.harmony, features = myCAF[i], cols = mycols,reduction = "tsne", pt.size = .1) + 
    NoLegend() + mytheme
}
p.feature.comb.myCAF <- cowplot::plot_grid(plotlist = p.feature.myCAF, ncol = 3, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.myCAF.pdf"),p.feature.comb.myCAF,width = 8,height = 6)
# ductal_cell
p.feature.ductal_cell <- list()
ductal.match <- match(ductal_cell,rownames(total.integrated.harmony@assays$RNA))
ductal_cell <- rownames(total.integrated.harmony@assays$RNA)[ductal.match]
ductal_cell <- ductal_cell[!is.na(ductal_cell)]
for(i in seq_along(ductal_cell)){
  p.feature.ductal_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = ductal_cell[i],cols = mycols,reduction = "tsne",pt.size = .2) +
    NoLegend() + mytheme
}
p.feature.comb.ductal <- cowplot::plot_grid(plotlist = p.feature.ductal_cell, ncol = 4, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.ductal.pdf"),p.feature.comb.ductal,width = 11,height = 6)
# Acinar_cell
p.feature.Acinar_cell <- list()
Acinar.match <- match(Acinar_cell,rownames(total.integrated.harmony@assays$RNA))
Acinar_cell <- rownames(total.integrated.harmony@assays$RNA)[Acinar.match]
Acinar_cell <- Acinar_cell[!is.na(Acinar_cell)]
for(i in seq_along(Acinar_cell)){
  p.feature.Acinar_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Acinar_cell[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.acinar <- cowplot::plot_grid(plotlist = p.feature.Acinar_cell, ncol = 3, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.acinar.pdf"),p.feature.comb.acinar,width = 8,height = 6)
# Endocrine_cell
p.feature.Endocrine_cell <- list()
Endocrine.match <- match(Endocrine_cell,rownames(total.integrated.harmony@assays$RNA))
Endocrine_cell <- rownames(total.integrated.harmony@assays$RNA)[Endocrine.match]
Endocrine_cell <- Endocrine_cell[!is.na(Endocrine_cell)]
for(i in seq_along(Endocrine_cell)){
  p.feature.Endocrine_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Endocrine_cell[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.endocrine <- cowplot::plot_grid(plotlist = p.feature.Endocrine_cell, ncol = 2, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.endocrine.pdf"),p.feature.comb.endocrine,width = 6,height = 6)
# Stellate_cell
p.feature.Stellate_cell <- list()
Stellate.match <- match(Stellate_cell,rownames(total.integrated.harmony@assays$RNA))
Stellate_cell <- rownames(total.integrated.harmony@assays$RNA)[Stellate.match]
Stellate_cell <- Stellate_cell[!is.na(Stellate_cell)]
for(i in seq_along(Stellate_cell)){
  p.feature.Stellate_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Stellate_cell[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.stellate <- cowplot::plot_grid(plotlist = p.feature.Stellate_cell, ncol = 2, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.stellate.pdf"),p.feature.comb.stellate,width = 6,height = 6)
# Fibroblast_cell
p.feature.Fibroblast_cell <- list()
Fibroblast.match <- match(Fibroblast_cell,rownames(total.integrated.harmony@assays$RNA))
Fibroblast_cell <- rownames(total.integrated.harmony@assays$RNA)[Fibroblast.match]
Fibroblast_cell <- Fibroblast_cell[!is.na(Fibroblast_cell)]
for(i in seq_along(Fibroblast_cell)){
  p.feature.Fibroblast_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Fibroblast_cell[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.fibroblast <- cowplot::plot_grid(plotlist = p.feature.Fibroblast_cell, ncol = 5, nrow = 4, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.fibroblast.pdf"),p.feature.comb.fibroblast,width = 14,height = 12)
# Endothelial_cell
p.feature.Endothelial_cell <- list()
Endothelial.match <- match(Endothelial_cell,rownames(total.integrated.harmony@assays$RNA))
Endothelial_cell <- rownames(total.integrated.harmony@assays$RNA)[Endothelial.match]
Endothelial_cell <- Endothelial_cell[!is.na(Endothelial_cell)]
for(i in seq_along(Endothelial_cell)){
  p.feature.Endothelial_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Endothelial_cell[i],cols = mycols,
                                                reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.endothelial <- cowplot::plot_grid(plotlist = p.feature.Endothelial_cell, ncol = 3, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.endothelial.pdf"),p.feature.comb.endothelial,width = 8,height = 6)
# Mast_cell
p.feature.Mast_cell <- list()
Mast.match <- match(Mast_cell,rownames(total.integrated.harmony@assays$RNA))
Mast_cell <- rownames(total.integrated.harmony@assays$RNA)[Mast.match]
Mast_cell <- Mast_cell[!is.na(Mast_cell)]
for(i in seq_along(Mast_cell)){
  p.feature.Mast_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Mast_cell[i],cols = mycols,
                                          reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.mast <- cowplot::plot_grid(plotlist = p.feature.Mast_cell, ncol = 2, nrow = 1, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.mast.pdf"),p.feature.comb.mast,width = 6,height = 3)
# Macrophage_cell
p.feature.Macrophage_cell <- list()
Macrophage.match <- match(Macrophage_cell,rownames(total.integrated.harmony@assays$RNA))
Macrophage_cell <- rownames(total.integrated.harmony@assays$RNA)[Macrophage.match]
Macrophage_cell <- Macrophage_cell[!is.na(Macrophage_cell)]
for(i in seq_along(Macrophage_cell)){
  p.feature.Macrophage_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Macrophage_cell[i],
                                                cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.macrophage <- cowplot::plot_grid(plotlist = p.feature.Macrophage_cell, ncol = 3, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.macrophage.pdf"),p.feature.comb.macrophage,width = 8,height = 6)
# T_cells
p.feature.T_cells <- list()
T.match <- match(T_cells,rownames(total.integrated.harmony@assays$RNA))
T_cells <- rownames(total.integrated.harmony@assays$RNA)[T.match]
T_cells <- T_cells[!is.na(T_cells)]
for(i in seq_along(T_cells)){
  p.feature.T_cells[[i]] <- FeaturePlot(total.integrated.harmony,features = T_cells[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.T <- cowplot::plot_grid(plotlist = p.feature.T_cells, ncol = 4, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.T.pdf"),p.feature.comb.T,width = 11,height = 6)
# B_cell
p.feature.B_cell <- list()
B.match <- match(B_cell,rownames(total.integrated.harmony@assays$RNA))
B_cell <- rownames(total.integrated.harmony@assays$RNA)[B.match]
B_cell <- B_cell[!is.na(B_cell)]
for(i in seq_along(B_cell)){
  p.feature.B_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = B_cell[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.B <- cowplot::plot_grid(plotlist = p.feature.B_cell, ncol = 4, nrow = 3, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.B.pdf"),p.feature.comb.B,width = 11,height = 9)
# Epithelial_cell
p.feature.Epithelial_cell <- list()
Epithelial.match <- match(Epithelial_cell,rownames(total.integrated.harmony@assays$RNA))
Epithelial_cell <- rownames(total.integrated.harmony@assays$RNA)[Epithelial.match]
Epithelial_cell <- Epithelial_cell[!is.na(Epithelial_cell)]
for(i in seq_along(Epithelial_cell)){
  p.feature.Epithelial_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Epithelial_cell[i],
                                                cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.epithelial <- cowplot::plot_grid(plotlist = p.feature.Epithelial_cell, ncol = 5, nrow = 3, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.epithelial.pdf"),p.feature.comb.epithelial,width = 14,height = 9)
# Lymphocytes
p.feature.Lymphocytes <- list()
Lymphocytes.match <- match(Lymphocytes,rownames(total.integrated.harmony@assays$RNA))
Lymphocytes <- rownames(total.integrated.harmony@assays$RNA)[Lymphocytes.match]
Lymphocytes <- Lymphocytes[!is.na(Lymphocytes)]
for(i in seq_along(Lymphocytes)){
  p.feature.Lymphocytes[[i]] <- FeaturePlot(total.integrated.harmony,features = Lymphocytes[i],
                                            cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.lymphocytes <- cowplot::plot_grid(plotlist = p.feature.Lymphocytes, ncol = 3, nrow = 1, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.lymphocytes.pdf"),p.feature.comb.lymphocytes,width = 8,height = 3)
# DCs
p.feature.DCs <- list()
DCs.match <- match(DCs,rownames(total.integrated.harmony@assays$RNA))
DCs <- rownames(total.integrated.harmony@assays$RNA)[DCs.match]
DCs <- DCs[!is.na(DCs)]
for(i in seq_along(DCs)){
  p.feature.DCs[[i]] <- FeaturePlot(total.integrated.harmony,features = DCs[i],cols = mycols,
                                    reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.DCs <- cowplot::plot_grid(plotlist = p.feature.DCs, ncol = 5, nrow = 3, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.DCs.pdf"),p.feature.comb.DCs,width = 14,height = 9)
# EMT_like_cell
p.feature.EMT_like_cell <- list()
EMT_like.match <- match(EMT_like_cell,rownames(total.integrated.harmony@assays$RNA))
EMT_like_cell <- rownames(total.integrated.harmony@assays$RNA)[EMT_like.match]
EMT_like_cell <- EMT_like_cell[!is.na(EMT_like_cell)]
for(i in seq_along(EMT_like_cell)){
  p.feature.EMT_like_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = EMT_like_cell[i],
                                              cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.EMT_like <- cowplot::plot_grid(plotlist = p.feature.EMT_like_cell, ncol = 3, nrow = 1, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.EMT_like.pdf"),p.feature.comb.EMT_like,width = 8,height = 3)
# Myeloid_cell
p.feature.Myeloid_cell <- list()
Myeloid.match <- match(Myeloid_cell,rownames(total.integrated.harmony@assays$RNA))
Myeloid_cell <- rownames(total.integrated.harmony@assays$RNA)[Myeloid.match]
Myeloid_cell <- Myeloid_cell[!is.na(Myeloid_cell)]
for(i in seq_along(Myeloid_cell)){
  p.feature.Myeloid_cell[[i]] <- FeaturePlot(total.integrated.harmony,features = Myeloid_cell[i],
                                             cols = mycols,reduction = "tsne",pt.size = .7) +
    NoLegend() + mytheme
}
p.feature.comb.myeloid <- cowplot::plot_grid(plotlist = p.feature.Myeloid_cell, ncol = 4, nrow = 3, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.myeloid.pdf"),p.feature.comb.myeloid,width = 11,height = 9)
# NK
p.feature.NK <- list()
NK.match <- match(NK,rownames(total.integrated.harmony@assays$RNA))
NK <- rownames(total.integrated.harmony@assays$RNA)[NK.match]
NK <- NK[!is.na(NK)]
for(i in seq_along(NK)){
  p.feature.NK[[i]] <- FeaturePlot(total.integrated.harmony,features = NK[i],cols = mycols,reduction = "tsne",pt.size = .1) +
    NoLegend() + mytheme
}
p.feature.comb.NK <- cowplot::plot_grid(plotlist = p.feature.NK, ncol = 3, nrow = 2, align = "h")
ggsave(paste0("figures/celltypes/markers/","p.feature.NK.pdf"),p.feature.comb.NK,width = 8,height = 6)

## ----------------- DotPlot ----------------------------------

mytheme <- theme(panel.border = element_rect(size = 1.5,colour = "grey75"),
                 axis.text = element_text(size = 12,face = "bold",colour = "grey45"),
                 axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
                 axis.title = element_text(size = 14,face = "bold",colour = "grey30"),
                 legend.position = "none")
mycols <- c("grey88","DarkCyan")
# iCAF
p.dot.iCAF <- DotPlot(total.integrated.harmony,features=iCAF,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.iCAF.pdf"),p.dot.iCAF,width = 4,height = 10)
# myCAF
p.dot.myCAF <- DotPlot(total.integrated.harmony,features=myCAF,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.myCAF.pdf"),p.dot.myCAF,width = 3,height = 10)
# ductal_cell
p.dot.ductal <- DotPlot(total.integrated.harmony,features=ductal_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.ductal.pdf"),p.dot.ductal,width = 3.5,height = 10)
# Acinar_cell
p.dot.acinar <- DotPlot(total.integrated.harmony,features=Acinar_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.acinar.pdf"),p.dot.acinar,width = 3.5,height = 10)
# Endocrine_cell
p.dot.endocrine <- DotPlot(total.integrated.harmony,features = Endocrine_cell,cols = mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.endocrine.pdf"),p.dot.endocrine,width = 2,height = 10)
# Stellate_cell
p.dot.stellate <- DotPlot(total.integrated.harmony,features=Stellate_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.stellate.pdf"),p.dot.stellate,width = 2,height = 10)
# Fibroblast_cell
p.dot.fibroblast <- DotPlot(total.integrated.harmony,features=Fibroblast_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.fibroblast.pdf"),p.dot.fibroblast,width = 8,height = 10)
# Endothelial_cell
p.dot.endothelial <- DotPlot(total.integrated.harmony,features=Endothelial_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.endothelial.pdf"),p.dot.endothelial,width = 2.5,height = 10)
# Mast_cell
p.dot.mast <- DotPlot(total.integrated.harmony,features=Mast_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.mast.pdf"),p.dot.mast,width = 2,height = 10)
# Macrophage_cell
p.dot.macrophage <- DotPlot(total.integrated.harmony,features=Macrophage_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.macrophage.pdf"),p.dot.macrophage,width = 2.8,height = 10)
# T_cells
p.dot.T <- DotPlot(total.integrated.harmony,features=T_cells,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.T.pdf"),p.dot.T,width = 3.5,height = 10)
# B_cell
p.dot.B <- DotPlot(total.integrated.harmony,features=B_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.B.pdf"),p.dot.B,width = 4.4,height = 10)
# Epithelial_cell
p.dot.epithelial <- DotPlot(total.integrated.harmony,features=Epithelial_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.epithelial.pdf"),p.dot.epithelial,width = 4.8,height = 10)
# Lymphocytes
p.dot.lymphocytes <- DotPlot(total.integrated.harmony,features=Lymphocytes,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.lymphocytes.pdf"),p.dot.lymphocytes,width = 2,height = 10)
# DCs
p.dot.DCs <- DotPlot(total.integrated.harmony,features=DCs,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.DCs.pdf"),p.dot.DCs,width = 4.8,height = 10)

# EMT_like_cell
p.dot.EMT_like <- DotPlot(total.integrated.harmony,features=EMT_like_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.EMT_like.pdf"),p.dot.EMT_like,width = 2,height = 10)
# Myeloid_cell
p.dot.myeloid <- DotPlot(total.integrated.harmony,features=Myeloid_cell,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.myeloid.pdf"),p.dot.myeloid,width = 4.5,height = 10)
# NK
p.dot.NK <- DotPlot(total.integrated.harmony,features=NK,cols=mycols) + 
  labs(x = "Markers", y = "Clusters") + theme_bw() + mytheme
ggsave(paste0("figures/celltypes/markers/","p.dot.NK.pdf"),p.dot.NK,width = 3,height = 10)

figures <- list(
  p.feature.comb.iCAF = p.feature.comb.iCAF,
  p.feature.comb.myCAF = p.feature.comb.myCAF,
  p.feature.comb.ductal = p.feature.comb.ductal,
  p.feature.comb.acinar = p.feature.comb.acinar,
  p.feature.comb.endocrine = p.feature.comb.endocrine,
  p.feature.comb.stellate = p.feature.comb.stellate,
  p.feature.comb.fibroblast = p.feature.comb.fibroblast,
  p.feature.comb.endothelial = p.feature.comb.endothelial,
  p.feature.comb.mast = p.feature.comb.mast,
  p.feature.comb.macrophage = p.feature.comb.macrophage,
  p.feature.comb.T = p.feature.comb.T,
  p.feature.comb.B = p.feature.comb.B,
  p.feature.comb.epithelial = p.feature.comb.epithelial,
  p.feature.comb.lymphocytes = p.feature.comb.lymphocytes,
  p.feature.comb.DCs = p.feature.comb.DCs,
  p.feature.comb.EMT_like = p.feature.comb.EMT_like,
  p.feature.comb.myeloid = p.feature.comb.myeloid,
  p.feature.comb.NK = p.feature.comb.NK,
  p.dot.iCAF = p.dot.iCAF,
  p.dot.myCAF = p.dot.myCAF,
  p.dot.ductal = p.dot.ductal,
  p.dot.acinar = p.dot.acinar,
  p.dot.endocrine = p.dot.endocrine,
  p.dot.endocrine = p.dot.endocrine,
  p.dot.fibroblast = p.dot.fibroblast,
  p.dot.endothelial = p.dot.endothelial,
  p.dot.mast = p.dot.mast,
  p.dot.macrophage = p.dot.macrophage,
  p.dot.T = p.dot.T,
  p.dot.B = p.dot.B,
  p.dot.epithelial = p.dot.epithelial,
  p.dot.lymphocytes = p.dot.lymphocytes,
  p.dot.DCs = p.dot.DCs,
  p.dot.EMT_like = p.dot.EMT_like,
  p.dot.myeloid = p.dot.myeloid,
  p.dot.NK = p.dot.NK
)
save(figures, file = "figures/figures.rda")

## ----------------- Cell type annotation ----------------------------------

# Label the clusters for cell populations
current.clust.ids <- c(0:39)
new.clust.ids <- c("T cells","Myeloid","T cells","Epithelial","Myeloid","T cells","NK cells","Myeloid",
                   "Myeloid","T cells","Myeloid","T cells","Fibroblast","Epithelial","B cells","Acinar cells",
                   "Myeloid","Myeloid","Mast cells","Myeloid","B cells","Epithelial","Myeloid","Stellate cells",
                   "Myeloid","Epithelial","Endothelial","Endothelial","T cells","B cells","T cells","Acinar cells",
                   "Unknown","B cells","Myeloid","Acinar cells","Epithelial","T cells","B cells","Epithelial")
total.integrated.harmony$CellTypes <- plyr::mapvalues(total.integrated.harmony$seurat_clusters, from = current.clust.ids, to = new.clust.ids)

mypal <- c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
           "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699","#E78AC3")
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "bottom",
                 plot.title = element_blank(),
                 legend.text = element_text(size = 14))
p.ct <- DimPlot(total.integrated.harmony,reduction = "tsne",group.by = "CellTypes",cols = mypal) + mytheme
ggsave(paste0("figures/celltypes/","celltypes.pdf"),p.ct, width = 6,height = 6.8)

p.raw <- DimPlot(total.integrated.harmony,reduction = "tsne",group.by = "seurat_clusters",label = T) + 
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/","raw.pdf"),p.raw, width = 6,height = 6)

## ---- cell type highlight ----
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 plot.title = element_blank(),
                 legend.position = "none")

p.acinar <- DimPlot(total.integrated.harmony, reduction = "tsne", cols.highlight = "#B8793C",
                    cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Acinar cells"))) + mytheme
ggsave(paste0("figures/celltypes/acinar.pdf"),p.acinar,width = 4,height = 4)
p.B <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#43D9FD",
               cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "B cells"))) + mytheme
ggsave(paste0("figures/celltypes/B.pdf"),p.B,width = 4,height = 4)
p.NK <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#FD6741",
               cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "NK cells"))) + mytheme
ggsave(paste0("figures/celltypes/NK.pdf"),p.NK,width = 4,height = 4)
p.T <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#437AFD",
                cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "T cells"))) + mytheme
ggsave(paste0("figures/celltypes/T.pdf"),p.T,width = 4,height = 4)
p.epi <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#43FD69",
               cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Epithelial"))) + mytheme
ggsave(paste0("figures/celltypes/Epithelial.pdf"),p.epi,width = 4,height = 4)
p.mast <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#669865",
                   cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Mast cells"))) + mytheme
ggsave(paste0("figures/celltypes/Mast.pdf"),p.mast,width = 4,height = 4)
p.mye <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#FEDA88",
                  cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Myeloid"))) + mytheme
ggsave(paste0("figures/celltypes/Myeloid.pdf"),p.mye,width = 4,height = 4)
p.endo <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#7E6699",
                 cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Endothelial"))) + mytheme
ggsave(paste0("figures/celltypes/Endothelial.pdf"),p.endo,width = 4,height = 4)
p.ste <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#993333",
                  cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Stellate cells"))) + mytheme
ggsave(paste0("figures/celltypes/Stellate.pdf"),p.ste,width = 4,height = 4)
p.fibro <- DimPlot(total.integrated.harmony,reduction = "tsne",cols.highlight = "#C543FD",
                 cells.highlight = rownames(subset(total.integrated.harmony@meta.data, CellTypes == "Fibroblast"))) + mytheme
ggsave(paste0("figures/celltypes/Fibroblast.pdf"),p.fibro,width = 4,height = 4)

## DE genes heatmap
Idents(total.integrated.harmony) <- total.integrated.harmony$CellTypes
all.markers <- FindAllMarkers(object = total.integrated.harmony, min.pct = 0.25, only.pos = TRUE)
save(all.markers,file = "output/all.markers.celltypes.rda")
top20 <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top20.exp <- AverageExpression(total.integrated.harmony,assays = "RNA",features = unique(top20$gene))
tmp = top20.exp$RNA
tmp <- log2(tmp+1)
tmp[tmp>5] = 5
# gaps <- c(20,40,60,77,97,117,137,156,168,187)
colorRampPalette(c("#4393C3", "white", "#D6604D"))(100)
pheatmap(tmp,cluster_rows = F,cluster_cols = F,show_rownames = F,
         fontsize = 14,angle_col = 45, border = F,
         color = rev(RColorBrewer::brewer.pal(10,"RdYlBu")),
         filename = paste0("figures/celltypes/","diff.gene.expr.pdf"),
         width = 8,height = 12)

## ------------ TME --------------------

Idents(total.integrated.harmony) <- total.integrated.harmony$CellTypes
mypal <- c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
           "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699","#E78AC3")
mypal2 <- c("#c14089","#6f5553","#E5C494","#738f4c","#bb6240","#66C2A5","#2dfd29","#0c0fdc")

# split by sampleName
p3.1.1 <- DimPlot(total.integrated.harmony, reduction = "tsne", label = T, pt.size = 0.7, split.by = "SampleName", ncol = 5, cols = mypal)
p3.1.2 <- DimPlot(total.integrated.harmony,reduction = "tsne",label = F,pt.size = 0.2,group.by = "SampleName")
ggsave(paste0("figures/TME/","tsne.splitBy.SampleName.pdf"),p3.1.1,width = 21,height = 30)
ggsave(paste0("figures/TME/","tsne.groupBy.SampleName.pdf"),p3.1.2,width = 9,height = 5)
# split/group by clinical stage
p3.2.1 <- DimPlot(total.integrated.harmony,reduction = "tsne", label = T, pt.size = 0.2, split.by = "ClinicalStage", ncol = 3, cols = mypal)
p3.2.2 <- DimPlot(total.integrated.harmony,reduction = "tsne", label = T, pt.size = 0.2, group.by = "ClinicalStage", cols = mypal2[1:5])
ggsave(paste0("figures/TME/","tsne.splitBy.ClinicalStage.pdf"),p3.2.1, width = 12,height = 8)
ggsave(paste0("figures/TME/","tsne.groupBy.ClinicalStage.pdf"),p3.2.2, width = 7, height = 6)
# split by class (primary, metastatic, paracancerous)
p3.3.1 <- DimPlot(total.integrated.harmony,reduction = "tsne",label = T,pt.size = 0.2, split.by = "Class", ncol = 2, cols = mypal)
p3.3.2 <- DimPlot(total.integrated.harmony,reduction = "tsne",label = T,pt.size = 0.2, group.by = "Class", cols = mypal2[1:4])
ggsave(paste0("figures/TME/","tsne.splitBy.Class.pdf"),p3.3.1, width = 9,height = 8)
ggsave(paste0("figures/TME/","tsne.groupBy.Class.pdf"),p3.3.2, width = 7,height = 6)
# split by tissue (tissue vs. pbmc)
p3.4.1 <- DimPlot(total.integrated.harmony,reduction = "tsne",label = T,pt.size = 0.2,split.by = "Tissue",ncol = 2, cols = mypal)
p3.4.2 <- DimPlot(total.integrated.harmony,reduction = "tsne",label = T,pt.size = 0.2, group.by = "Tissue", cols = mypal2[1:2])
ggsave(paste0("figures/TME/","tsne.splitBy.Tissue.pdf"),p3.4.1,width = 9,height = 4)
ggsave(paste0("figures/TME/","tsne.groupBy.Tissue.pdf"),p3.4.2, width = 7, height = 6)
# split by patient (PDAC vs Normal)
p3.5.1 <- DimPlot(total.integrated.harmony, reduction = "tsne", label = T, pt.size = 0.2, split.by = "Patient", ncol = 2, cols = mypal)
p3.5.2 <- DimPlot(total.integrated.harmony, reduction = "tsne", label = T, pt.size = 0.2, group.by = "Patient", cols = mypal2[1:2])
ggsave(paste0("figures/TME/","tsne.splitBy.Patient.pdf"),p3.5.1,width = 9,height = 4)
ggsave(paste0("figures/TME/","tsne.groupBy.Patient.pdf"),p3.5.2,width = 7,height = 6)

## --- Stacked Barplot (cell proportion) ----
# SampleName vs. Celltypes
p.sn <- ggplot(total.integrated.harmony@meta.data, aes(x = SampleName, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.direction = "vertical", legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)) +
  labs(x = "Samples", y = "Cell Fraction") +
  guides(fill = guide_legend(ncol = 1))
ggsave(paste0("figures/TME/","cell.prop_vs_Samples.pdf"),p.sn,width = 8,height = 12)

# Patient vs. Celltypes
p.pt <- ggplot(total.integrated.harmony@meta.data, aes(x = Patient, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 14),legend.title = element_text(size = 15)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Patient", y = "Cell Fraction")
ggsave(paste0("figures/TME/","cell.prop_vs_Patient.pdf"),p.pt, width = 14,height = 3)

# ClinicalStage vs. Celltypes
p.cs <- ggplot(total.integrated.harmony@meta.data, aes(x = ClinicalStage, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 14),legend.title = element_text(size = 15)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "ClinicalStage", y = "Cell Fraction")
ggsave(paste0("figures/TME/","cell.prop_vs_ClinicalStage.pdf"),p.cs, width = 14,height = 3)

# Tissue vs. Celltypes
p.ts <- ggplot(total.integrated.harmony@meta.data, aes(x = Tissue, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 14),legend.title = element_text(size = 15)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Tissue", y = "Cell Fraction")
ggsave(paste0("figures/TME/","cell.prop_vs_Tissue.pdf"),p.ts, width = 14,height = 3)

# Tissue vs. Celltypes
p.cls <- ggplot(total.integrated.harmony@meta.data, aes(x = Class, fill = CellTypes)) +
  geom_bar(position = "fill") + theme_bw() + coord_flip() +
  scale_fill_manual(values = mypal) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,face = "bold",colour = "grey50"),
        axis.title = element_text(size = 14,face = "bold",colour = "grey25"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 14),legend.title = element_text(size = 15)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Class", y = "Cell Fraction")
ggsave(paste0("figures/TME/","cell.prop_vs_Class2.pdf"),p.cls, width = 18,height = 3)


## ---------------- TCGA check -------------------

clust.ct <- as.matrix(table(total.integrated.harmony$seurat_clusters,total.integrated.harmony$CellTypes))
cls.name <- NULL
for(i in 1:nrow(clust.ct)){
  for(j in 1:ncol(clust.ct)){
    if(clust.ct[i,j]>0) ct.name <- paste0("C",rownames(clust.ct)[i],"(",colnames(clust.ct)[j],")")
  }
  cls.name <- c(cls.name,ct.name)
}

load("output/all.markers.rda")

all.markers.1 <- all.markers
all.markers.1$clusters <- plyr::mapvalues(all.markers.1$cluster, from = seq(0,39,1),to = cls.name)
paad_mrna <- xlsx::read.xlsx("input/tcga/PAAD_mrna.xlsx",sheetIndex = 1)
all.markers.1$gene <- gsub("IL8","CXCL8",all.markers.1$gene)
all.markers.1$gene <- gsub("EMR2","ADGRE2",all.markers.1$gene)
all.markers.1$gene <- gsub("GPR110","ADGRF1",all.markers.1$gene)
all.markers.1$gene <- gsub("C19orf10","MYDGF",all.markers.1$gene)
all.markers.1$gene <- gsub("GPR56","ADGRG1",all.markers.1$gene)
all.markers.1$gene <- gsub("C19orf59","MCEMP1",all.markers.1$gene)
all.markers.1$gene <- gsub("MRC1L1","MRC1",all.markers.1$gene)
all.markers.1$gene <- gsub("KIAA1598","SHTN1",all.markers.1$gene)
all.markers.1$gene <- gsub("PPAP2B","PLPP3",all.markers.1$gene)
all.markers.1$gene <- gsub("ERO1LB","ERO1B",all.markers.1$gene)
all.markers.1$gene <- gsub("PPAP2A","PLPP1",all.markers.1$gene)
all.markers.1$gene <- gsub("GPR116","ADGRF5",all.markers.1$gene)
all.markers.1$gene <- gsub("C8orf47","ERICH5",all.markers.1$gene)
all.markers.1$gene <- gsub("DARC","ACKR1",all.markers.1$gene)
all.markers.1$gene <- gsub("ELTD1","ADGRL4",all.markers.1$gene)
all.markers.1$gene <- gsub("C9orf89","CARD19",all.markers.1$gene)
all.markers.1$gene <- gsub("PPAPDC1B","PLPP5",all.markers.1$gene)

all.markers.1 <- all.markers.1[all.markers.1$gene %in% paad_mrna$Updated.Name,]
paad_mrna.1 <- paad_mrna[which(paad_mrna$Updated.Name %in% all.markers.1$gene),][,c(4:9)]
colnames(paad_mrna.1)[1] <- "gene"
dat <- dplyr::right_join(x = all.markers.1, y = paad_mrna.1, by = "gene")
dat$clusters <- as.factor(dat$clusters)

mypal <- c("#437AFD","#FEDA88","#437AFD","#43FD69","#FEDA88","#437AFD","#FD6741","#FEDA88",
           "#FEDA88","#437AFD","#FEDA88","#437AFD","#C543FD","#43FD69","#43D9FD","#B8793C",
           "#FEDA88","#FEDA88","#669865","#FEDA88","#43D9FD","#43FD69","#FEDA88","#993333",
           "#FEDA88","#43FD69","#7E6699","#7E6699","#437AFD","#43D9FD","#437AFD","#B8793C",
           "#E685C1","#43D9FD","#FEDA88","#B8793C","#43FD69","#437AFD","#43D9FD","#43FD69")
# visualize cox coefficient
p.cox <- ggplot(data = dat[,c(8,9)],aes(x = clusters,y=Cox.coefficient,fill=clusters)) + 
  geom_hline(yintercept = 0,size=1,color="grey10") + geom_violindot(dots_size = 0.1) + 
  scale_fill_manual(values = mypal)+ theme_bw() +
  ggtitle("TCGA-PAAD (n=183)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "grey10",size = 1.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,hjust = 1,vjust = 1,colour = "grey10"),
        axis.title = element_text(size = 16, face = "bold",colour = "grey10"),
        axis.text.y = element_text(size = 14,colour = "grey10"),
        plot.title = element_text(size = 14,colour = "grey10"))
ggsave(paste0("figures/celltypes/","TCGA.Cox.pdf"),p.cox,width = 16,height = 3.5)

sum.len <- 0
ave.cox <- NULL
for(i in c(0:39)){
  len <- length(dat$cluster[dat$cluster==i])
  ave.cox.tmp <- sum(dat$Cox.coefficient[(sum.len+1):(sum.len+len)])/len
  ave.cox.rep <- rep(ave.cox.tmp,len)
  ave.cox <- c(ave.cox,ave.cor.rep)
  sum.len <- sum.len + len
  
}
dat <- cbind(dat,ave.cor)

dat.1 <- dat[,c("clusters","ave.cox")]
dat.1 <- unique(dat.1)
# visualize average cox coefficient
p.ave.cox <- ggpubr::ggbarplot(data = dat.1, x = "clusters",y = "ave.cox",color = mypal,fill = mypal) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 14)) +
  geom_hline(yintercept = 0,size=0.7,color="grey10")
ggsave(paste0("figures/celltypes/","TCGA.ave.Cox.pdf"),p.ave.cox,width = 18,height = 3.5)

save(dat, file = "output/all.markers.with.cox.rda")

# ---- Check gene ----

FeaturePlot(total.integrated.harmony,features = Epithelial_cell[i],
            cols = mycols,reduction = "tsne",pt.size = .1)
mycols <- c("grey88","DarkCyan")
mytheme <- theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank())

## ---- B cell ----
# BCL11A, DNASE1L3, CD79A, BANK1
# p.BCL11A <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "BCL11A",
#                         cols = c("grey88","DarkCyan"), pt.size = 1) +
#   theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(paste0("figures/celltypes/features/","BCL11A.pdf"),p.BCL11A,width = 4,height = 4.2)
# # DNASE1L3
# p.DNASE1L3 <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "DNASE1L3",
#                         cols = c("grey88","DarkCyan"), pt.size = 1) +
#   theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(paste0("figures/celltypes/features/","DNASE1L3.pdf"),p.DNASE1L3,width = 4,height = 4.2)
# # CD79A
# p.CD79A <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "CD79A",
#                           cols = c("grey88","DarkCyan"), pt.size = 1) +
#   theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(paste0("figures/celltypes/features/","CD79A.pdf"),p.CD79A,width = 4,height = 4.2)
# 
# p.BANK1 <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "BANK1",
#                        cols = c("grey88","DarkCyan"), pt.size = 1) +
#   theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(paste0("figures/celltypes/features/","BANK1.pdf"),p.BANK1,width = 4,height = 4.2)

# validation the expression of IL35,BTK and HIF1A that were reported in 
# B cell that were associated with poor PDAC diagnosis
# p.IL35 <- VlnPlot(total.integrated.harmony,group.by = "CellTypes",features = "EBI3")+NoLegend()+
#   theme(axis.title.x = element_blank())
# ggsave(paste0("figures/celltypes/features/","IL35.pdf"),p.IL35,width = 6,height = 2.5)
# 
# p.BTK <- VlnPlot(total.integrated.harmony,group.by = "CellTypes",features = "BTK",pt.size = 0)+NoLegend()+
#   theme(axis.title.x = element_blank())
# ggsave(paste0("figures/celltypes/features/","BTK.pdf"),p.BTK,width = 6,height = 2.5)
# 
# p.HIF1A <- VlnPlot(total.integrated.harmony,group.by = "CellTypes",features = "HIF1A",pt.size = 0)+NoLegend()+
#   theme(axis.title.x = element_blank())
# ggsave(paste0("figures/celltypes/features/","HIF1A.pdf"),p.HIF1A,width = 6,height = 2.5)



## ---- Myeloid cell -----

# p.PTX3 <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "PTX3",
#                       cols = c("grey88","DarkCyan"), pt.size = 1) +
#   theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# ggsave(paste0("figures/celltypes/features/","PTX3.pdf"),p.PTX3,width = 4,height = 4.2)

p.SPP1 <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "SPP1",
                      cols = c("grey88","DarkCyan"), pt.size = 1) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/features/","SPP1.pdf"),p.SPP1,width = 4,height = 4.2)

p.MARCO <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "MARCO",
                      cols = c("grey88","DarkCyan"), pt.size = 1) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/features/","MARCO.pdf"),p.MARCO,width = 4,height = 4.2)

p.APOE <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "APOE",
                       cols = c("grey88","DarkCyan"), pt.size = 1) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/features/","APOE.pdf"),p.APOE,width = 4,height = 4.2)

p.CD68 <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "CD68",
            cols = c("grey88","DarkCyan"), pt.size = 1) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/features/","CD68.pdf"),p.CD68,width = 4,height = 4.2)

p.SIRPA <- FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "SIRPA",
                      cols = c("grey88","DarkCyan"), pt.size = 1) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
ggsave(paste0("figures/celltypes/features/","SIRPA.pdf"),p.SIRPA,width = 4,height = 4.2)


## cell stemness calculation

stemness <- function(exprs,stem.sig.file){
  stem.sig <- read.delim(stem.sig.file,header = F,row.names = 1)
  comm.genes <- intersect(rownames(stem.sig), rownames(exprs))
  exprs <- exprs[comm.genes,]
  stem.sig <- stem.sig[comm.genes,]
  s <- apply(exprs,2,function(z) {cor(z, stem.sig, method = "sp", use = "complete.obs")})
  names(s) <- colnames(exprs)
  s <- s - min(s)
  s <- s/max(s)
  return(s)
}
DefaultAssay(total.integrated.harmony) <- "RNA"
total.integrated.harmony <- ScaleData(total.integrated.harmony)
stemness.scores <- stemness(exprs = GetAssayData(total.integrated.harmony,assay = "RNA",slot = "scale.data"),stem.sig.file = "input/pcbc-stemsig.tsv")
total.integrated.harmony$stemness.scores <- stemness.scores

p.stemness <-FeaturePlot(total.integrated.harmony,reduction = "tsne",features = "stemness.scores",
            cols = c("grey88","DarkCyan"), pt.size = 0.3) +
  theme(panel.border = element_rect(size = 1.5, colour = "grey75"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 14))
ggsave(paste0("figures/celltypes/","stemness.scores.pdf"),p.stemness,width = 6,height = 5.5)

## Immune checpoint 
load("input/signature_collection.rda")
IC.genes <- signature_collection$Immune_Checkpoint
IC.genes[33] <- "C10orf54"

mytheme <- theme(panel.border = element_rect(size = 1.5,colour = "grey75"),
                 axis.text = element_text(size = 12,face = "bold",colour = "grey45"),
                 axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
                 axis.title = element_text(size = 14,face = "bold",colour = "grey30"),
                 legend.position = "none")
mycols <- c("grey88","DarkCyan")
Idents(total.integrated.harmony) <- total.integrated.harmony$CellTypes
p.dot.iCAF <- DotPlot(total.integrated.harmony,features=IC.genes,cols=mycols) + 
  labs(x = "Immune Checkpoint Ligands/Receptors", y = "Celltypes") + theme_bw() + mytheme

IC.exp <- AverageExpression(total.integrated.harmony,assays = "RNA",features = IC.genes)
IC.exp <- log2(IC.exp$RNA+1)
IC.exp[IC.exp>3] = 3
tmp = t(IC.exp)
pheatmap(tmp,cluster_rows = F,border = NA,
         color = rev(RColorBrewer::brewer.pal(10,"RdYlBu")),
         angle_col = 45,fontsize = 14,
         cellheight = 30,cellwidth = 30,
         filename = paste0("figures/celltypes/ImmuneCheckpoint.pdf"),
         width = 16,height = 6)


