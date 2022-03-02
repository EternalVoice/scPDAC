**The pipeline will implement the following analysis contents: **

****

##### Figure 1 (Dissection of the tumor microenvironment in PDAC)

<img src="./Figures/Figure 1.jpg" style="zoom:50%;" />

**A.** Samples collected from [GSE155698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698); **B.** Visualization of single-cell RNA-seq data of 124,575 cells by tSNE; **C.** Single-cell resolution heatmap of top expressed genes for each cell type; **D.** Proportion of 11 cell types among 16 PDAC patients with tissue, PBMC, and paratumor, and 3 healthy PBMC samples.

****

##### Figure 2(Myeloid cells exert immune-suppressive potentials)

<img src="./Figures/Figure 2.jpg" style="zoom: 25%;" />

**A.** Genes exclusively expressed in C3 (epithelial/tumor cells) (ave.cor = 0.182) are associated with poor prognosis, where ave.cor presents the positively with favorable prognosis of PDAC, suggesting a immune-suppressive function of C1 (Myeloid, ave.cor = 0.023) cells in the microenvironment of PDAC.  **B.** Expression of immunosuppressive markers consists of SPP1, MARCO, APOE, CD56 and SIRPA. **C.** Myeloid cells have a relatively high stemness scores, algorithm to calculate stemness score was adopted from (scCancer, Guo et al., 2021). **D.** Immune checkpoint inhibitors, CD86, HAVCR2, CD48, C10orf54 (VSIR), LGALS9, TNFRSF14, LGALS3 have a relatively higher expression values in myeloid cells than in others.

****

##### Figure 3(Characterization of major myeloid cell lineages)

<img src="./Figures/Figure 3.jpg" style="zoom: 25%;" />

**A.** t-SNE plot showing 6 myeloid clusters of sample from PDAC patient. Proportion of each myeloid cell lineage from primary tumor, metastatic tumor, paracancerous tissue (**B**) and tissue and peripheral blood (**C**). **D.** Bubble plot showing selected cell type-specific markers across all clusters. Size of dots represents the fraction of cell expressing a particular marker, and intensity of color indicates the level of average mean expression. **E.** t-SNE plots showing expression of some immunosuppressive markers in the myeloid subclusters. **F-I.** Kaplan-Meier survival analysis of some myeloid subclusters, indicating macrophage, mast cells, CD14Mono, CD16Mono, and DCs.

****

##### Figure 4(Differentially expressed gene analysis (PDAC.tumor_vs_paratumor))

<img src="./Figures/Figure 4.jpg" style="zoom: 25%;" />

**A.** Tissue prevalence estimated by Ro/e score between tumor and normal samples. **B.** The differentially expressed genes for each cell type between tumor and paratumor samples from PDAC patients. Top 10 DE genes were shown, and the points dotted in red indicate significant genes. **C-G.** Differential pathway enriched in tumor and paratumor for each cell type.

****

##### Figure 5(Overlapped Differentially expressed genes associated with necroptosis)

<img src="./Figures/Figure 5.jpg" style="zoom: 25%;" />

Analysis of overlapped necroptosis-associated genes upregulated (**A**) or downregulated (**B**) in tumor for each cell type. **C.** The overlapped genes and their roles in necroptosis pathway. Genes in red are upregulated in tumor and in blue and downregulated in tumor. Arrows represent activation, and stacked line represents repression.

****

##### Figure 6(A novel immunological RTM population in paratumor)

<img src="./Figures/Figure 6.jpg" style="zoom: 25%;" />

**A.** Schema showing the procedures to identify the subclusters of macrophages. **B.** Bubble plot showing selected cell type-specific markers across all clusters. Size of dots represent the fraction of cells expressing a particular marker, and intensity of color indicates the level of average mean expression. **C.** t-SNE plots showing expression of GLUL and SQSTM1. **D.** t-SNE plots showing the subclusters of macrophages between tumor and paratumor samples. **E.** Violin plots showing a RTM (tissue resident macrophage, which is GLUL and SQSTM1 double negative) cluster-specific markers. **F.** Correlation between the GLUL-SQSTM1-RTM cluster-specific genes in CD8A. **G.** TCGA validation of the genes illustrated in **F.**

****

##### Figure 7(HSP90AA1+HSP90AB1+Mast cells are pro-immune)

<img src="./Figures/Figure 7.jpg" style="zoom: 25%;" />

**A.** Schema showing the procedures to distinguish mast cells by HSP90AA1 and HSP90AB1. **B.** Scatterplot showing the correlation between the HSP90AA1/HSP90AB1 expression and CD8A in this study. Spearman correlation. **C.** Violin plots showing the expression of HSP90AA1+HSP90AB1+Mast cluster-specific genes. **D.** Scatterplot showing the correlation between top 6 HSP90AA1+HSP90AB1+Mast cluster-specific genes and CD8A in this study.

****

##### Figure 8 (JAK3+TLR4+CD16 Mono cells are anti-immune)

<img src="./Figures/Figure 8.jpg" style="zoom: 25%;" />

**A.** Schema showing the procedures to distinguish mast cells by JAK3 and TLR4. **B.** Scatterplot showing the correlation between the JAK3/TLR4 expression and CD8A in this study. Spearman correlation. **C.** Violin plots showing the expression of JAK3+TLR4+Mast cluster-specific genes. **D.** Scatterplot showing the correlation between JAK3+TLR4+Mast cluster-specific genes and CD8A in this study.

****

##### Figure 9 (Detection of the cellular cluster in PDAC)

<img src="./Figures/Figure 9.jpg" style="zoom: 25%;" />

tSNE plots combing GLUL+SQSTM1+RTM (n = 998, in purple), HSP90AA1+HSP90AB1+Mast cells(n = 2,151, in blue), and JAK3+TLR4+CD16 Mono (n = 1,234, in orange) from PDAC (Peng et al., 2019)

****

##### Extended Figure 1 (Heterogeneity of TME in PDAC)

<img src="./Figures/Extended Figure 1.jpg" style="zoom: 25%;" />

**A.** t-SNE plots of cells from 39 samples profiled in this study, with each cell color coded to indicate the associated cell types. Proportion of 11 celltypes among 16 PDAC patients and healthy donors with primary tumor, metastatic tumor, paracancerous, and healthy (**B**), clinical stage I-IV (**C**), PDAC patient and healthy (**D**), tissue and PBMC (**E**).

****

##### Extended Figure 2 (B cell may play a tumor-suppressive role in PDAC)

<img src="./Figures/Extended Figure 2.jpg" style="zoom: 25%;" />

**A.** Top100 genes were used to calculate cox coefficients for each cluster using TCGA (PAAD, n = 183) data. **B.** t-SNE plots showing the expression *BCL11A* and *DNASE1L3*.**C.** Scatterplots showing the correlation between *BCL11A*, *DNASE1L3* and *CD8A* using TCGA (PAAD, n = 183) data. **D.** The bodymaps showing the expression of *BCL11A* and *DNASE1L3* between tumor and normal samples via GEPIA 2 (http://gepia2.cancer-pku.cn/#index). **E.** Kaplan-Meier overall survival analysis of high and low group of *BCL11A* (top) and *DNASE1L3* (bottom).**F.** RFS (Disease Free Survival) analysis of high and low group of *BCL11A* (top) and *DNASE1L3* (bottom).Hazards ratio was calculated based on Cox PH Model, and 95% CI (Confidence Interval) was applied.

****

##### Extended Figure 3(DE gene and pathway analysis for sample from PBMC)

<img src="./Figures/Extended Figure 3.jpg" style="zoom: 25%;" />

**A.** t-SNE plots of myeloid lineages of PBMC from healthy donors and PDAC patients. **B.** Proportion of each myeloid cell lineage from PDAC patients and healthy donors. **C.** Tissue prevalence estimated by Ro/e score between primary tumor, metastatic tumor, and paratumor (left), tissue and PBMC (bottom) from PDAC patients samples, PDAC patient and healthy donors of PBMCs (right). **D.** The differentially expressed genes for each celltype between PDAC patient and normal samples of PBMC. Top 10 DE genes were shown, and the points dotted in red indicate significant genes. **E-I.** Differential pathway enriched in PDAC and Normal of PBMC for each celltype.

****

##### Extended Figure 4 (Differential genes(PDAC.tissue tumor_vs_paratumor) and their associated pathway)

<img src="./Figures/Extended Figure 4.jpg" style="zoom:25%;" />

**A.** The differentially expressed genes for each celltype between tumor patient and paratumor from tissue of PDAC patients. Top 10 DE genes were shown, and the points dotted in red indicate significant genes. **B-F.** Differential pathway enriched in tumor and paratumor from PDAC tissue for each celltype.

****

##### Extended Figure 5 (Differential genes(PDAC. tissue_vs_pbmc) and their associated pathway)

<img src="./Figures/Extended Figure 5.jpg" style="zoom: 25%;" />

**A.** The differentially expressed genes for each celltype between tissue and PBMC from PDAC patients. Top 10 DE genes were shown, and the points dotted in red indicate significant genes. **B-F.** Differential pathway enriched in tissue and PBMC from PDAC patients for each celltype.

****

##### Extended Figure 6 (Visualization of overlapped necroptosis-associated DE genes)

<img src="./Figures/Extended Figure 6.jpg" style="zoom: 25%;" />

Violin plots showing overlapped necroptosis-associated genes (denoted in Figure 5) upregulated (A) and downregulated (B) in tumor. In each violin plot, paratumor, normal, and PBMC are in left, and tumor, PDAC, tissue are in right. *, p < 0.05, **, p < 0.01, ***, p < 0.001.

****

##### Extended Figure 7 (**The** **novel** immunological RTM has no direct correction with epithelial cells)

<img src="./Figures/Extended Figure 7.jpg" style="zoom: 25%;" />

**A.** Scatterplots showing the correlation between the expression of GLUL<sup>-</sup>SQSTM1<sup>-</sup> cluster-specific genes (denoted in Figure 6) and EPCAM. **B.** TCGA validation of the correlations. Data are from TCGA-PAAD (n = 178).

****

##### Extended Figure 8 (HSP90AA1<sup>-</sup>HSP90AB1<sup>-</sup>Mast cells are anti-immune)

<img src="./Figures/Extended Figure 8.jpg" style="zoom: 25%;" />

**A.** Scatterplot showing the correlation between the rest of HSP90AA1+HSP90AB1+Mast cluster-specific genes (showing in Figure 7) and CD8A in this study. **B.** Violin plots showing the expression of HSP90AA1-HSP90AB1-Mast cluster-specific genes. **C.** Scatterplots showing the correlation between HSP90AA1-HSP90AB1-Mast cluster-specific genes and CD8A in this study.

****

##### Extended Figure 9 (JAK3<sup>-</sup>TLR4<sup>-</sup>CD16 Mono. cells are pro-immune)

<img src="./Figures/Extended Figure 9.jpg" style="zoom: 25%;" />

**A.** Violin plots showing the expression of JAK3-TLR4-Mast cluster-specific high expressed genes  **B.** Scatterplot showing the correlation between the JAK3/TLR4 expression and CD8A in this study. Spearman correlation.

****

##### Extended Figure 10 (Several myeloid subsets may act as pro/anti-immune regulators in a non-necroptosis way)

<img src="./Figures/Extended Figure 10.jpg" style="zoom: 25%;" />

Schema showing the procedures to distinguish CD14 Mono. (A), CD16 Mono. (B), and DCs (C) by BIRC3, PPIA, and CHMP1B, respectively.
