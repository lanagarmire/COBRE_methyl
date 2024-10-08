# Multi-omics Analysis of Umbilical Cord Hematopoietic Stem Cells from a Multi-ethnic Cohort of Hawaii Reveals the Transgenerational Effect of Maternal Pre-Pregnancy Obesity <br>

**Background** <br> Maternal obesity is a health concern that may predispose newborns to a high risk of medical problems later in life. To understand the transgenerational effect of maternal obesity, we conducted a multi-omics study, using DNA methylation and gene expression in the CD34+/CD38-/Lin-umbilical cord blood hematopoietic stem cells (uHSCs) and metabolomics of the cord blood, all from a multi-ethnic cohort (n=72) from Kapiolani Medical Center for Women and Children in Honolulu, Hawaii (collected between 2016 and 2018).

**Citation** <br>
Du, Yuheng, Paula A. Benny, Yuchen Shao, Ryan J. Schlueter, Alexandra Gurary, Annette Lum-Jones, Cameron B. Lassiter, et al. n.d. “Multi-Omics Analysis of Umbilical Cord Hematopoietic Stem Cells from a Multi-Ethnic Cohort of Hawaii Reveals the Transgenerational Effect of Maternal Pre-Pregnancy Obesity.” *medRxiv*. [https://doi.org/10.1101/2024.07.27.24310936.](https://www.medrxiv.org/content/10.1101/2024.07.27.24310936v1) <br>

## Figure 1. Analysis Pipeline
![Analysis Pipeline](https://github.com/yhdu36/COBRE_methyl/blob/main/Figure_01/Figure_01_COBRE%20analysis%20pipeline.png?raw=true)


## Figure 2. Demographic statistics visualization

    # Input: clinical information pd
    cobre_mom_baby = read.csv('clinical_pd.csv')
    # Run Figure_02/cobre mother newborn statistics visualization.R
	ggplot(...)

## Figure 3. Differential methylation analysis
### Step 1: Data Acquisition
#### Option 1:  Preprocessing Raw Illumina 450K .idat files
    # Input: Processed .idat with ChAMP package in R
    Rscript Figure_03/1_load_cobre.r
#### Option 2: Use preprocessed methylation matrix
    library(GEOquery)
    library(limma)
    library(umap)
    # load series and platform data from GEO
    gset <- getGEO("GSE273075", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    myCombat <- exprs(gset)
	mypd <- pData(gset)
	
### Step 2: Source of Variance (SOV) Analysis 
#### Input
> **myCombat**: Preprocessed and normalized methylation matrix<br>
> **my_pd**:  clinical data
#### **Output**
> **Ftab**: Confounder names <br>
> **Fmean**: Avergaed F-statistics <br>
> **finalvars**: Significant counfounders selected to be adjusted with F-stat < threshold (**Default F-stat<1**).

    Rscript Figure_03/2_sov.R

### Step 3: Differential Methylation Analysis
#### Input
> **myCombat**: Preprocessed and normalized methylation matrix<br>
> **cobre_pd**:  clinical data<br>
> **cobre_sov_res**: Source of variance analysis results 
#### **Output**
>  **Plot**: Differential volcano plot<br>
>  **cobre_limma_sig**: Significant differential CpG sites <br>
>  **diff_cpg_gene**: Differential CpG sites with gene info annotation
 
    Rscript Figure_03/3_limma volcano.R
    
### Step 4: CpG region distribution

    Rscript Figure_03/gene pie distribution plot.R


## Figure 4. Functional enrichment analysis
### Step 1: KEGG enrcihment
#### Input
> **diff_cpg_gene**: Differential CpG sites with gene info annotation<br>
> **missmethyl_kegg**:  Subset of KEGG pathways with supergroups from 'Cellular Processes','Environmental Information Processing', 'Genetic Information Processing', 'Metabolism','Organismal Systems'
#### **Output**
> **gst.kegg.hypo.modified**: KEGG pathways enriched with hypomethylated sites <br>
> **gst.kegg.hyper.modified**: KEGG pathways enriched with hypermethylated sites

    Rscript Figure_04/1_KEGG_enrichment_missmethyl.r

### Step 2: KEGG immune and protein scores
#### Input
> **cobre_beta**: Methylation beta matrix<br>
> **cobre_pd**: Clinical data<br>
> **anno**: Illumina450K annotation file <br>
> **missmethyl_kegg**:  Subset of KEGG pathways with supergroups from 'Cellular Processes','Environmental Information Processing', 'Genetic Information Processing', 'Metabolism','Organismal Systems'
#### **Output**
> **Plot**: Immune and protein pathway scores comparison

    Figure_04/2_check_immune_protein_scores.ipynb

### Step 3: Stemness score
#### Input
> **cobre_beta**: Methylation beta matrix<br>
> **cobre_pd**: Clinical data<br>
> **anno**: Illumina450K annotation file
#### **Output**
> **SC.entropyscaled.Plasticity**: Sample-wise scaled stemness score calculated with Shannon entropy

    Figure_04/4_stemness.ipynb

### Step 4: Weighted Correlation Network Analysis (WGCNA)
#### Input
> **cpg_mean_promotor_tss_limma_adj**: Methylation confounder adjusted beta matrix, averaged on gene promoter regions (See: Figure_03/7_gsea_adjustbeta_methyl_to_gene.r)  <br>
> **cobre_pd**: Clinical data<br>
> **anno**: Illumina450K annotation file
#### **Output**
> **TOM**: Co-expressed networks with module annotation  <br>
> **modTraitCor**: Module and pheno trait correlation matrix

    Rscript Figure_04/5_WGCNA.R

### Step 4: Protein-Protein Interaction (PPI) Network
#### Input
> **diff_cpg_gene**: Differential CpG sites with gene info annotation
#### **Output**
> **KP_tss200_tss1500**: KEGG pathways enriched with significant protein-protein pairs  <br>
> **OrderAll_tss200_tss1500_sub**: Node information for Cytoscape

    Rscript Figure_04/6_PPI.R

## Figure 5. Multi-omics integration analysis
### Step 1: RNA-seq preprocessing

#### Option 1:  Preprocessing Raw fastq files
    # Obtain .bam files from STAR alignment
    # Use featureCount from subread to generate the gene count matrix
    Figure_05/RNA_seq/1_fc_cobre.sbat
    
#### Option 2: Use preprocessed methylation matrix
    # load series and platform data from GEO
    gset <- getGEO("GSE273075", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL20301", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    expr <- exprs(gset)

### Step 2: Metabolomics data preprocessing

 1. **Batch removal**

	Figure_05/01_3batch_preproc_vsn_combat.ipynb
	
 2. **Impute missing values**

	Figure_05/02_missing_impute_RF.R

### Step 3: DIABLO Multiomics Integration
#### Input
> **combat_edata3_vsn**: Metabolomics data after imputation and variance stabilizing normalization  <br>
> **dds**: RNA-seq counts for DESeq2  <br>
> **cobre_pd**: clinical information  <br>
> **cobre_methyl_sig**: Significant CpG sites  <br>
> **anno**: Illumina450K annotation file
#### **Output**
> **result.diablo.tcga**: DIABLO integration result

    Figure_05/03_diablo_3omics_3batch_vsn.ipynb

