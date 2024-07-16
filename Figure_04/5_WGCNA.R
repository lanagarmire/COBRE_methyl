library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
allowWGCNAThreads()

type = "unsigned"
corType = "pearson"
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

#saveRDS(cpg_mean_promotor_tss_limma_sub,'cpg_mean_promotor_tss_limma_sub.RDS')

cobre_pd = readRDS('/Users/Brainfactory/cobre_pd.rds')
#cobre_beta = readRDS('/Users/Brainfactory/cobre_beta.rds')
#cobre_beta_adjall = readRDS('/Users/Brainfactory/cobre_final_adj_betas_allparityhemo.rds')
cpg_mean_promotor_tss_limma_adj = readRDS('/Users/Brainfactory/cpg_mean_promotor_tss_limma_adj.RDS')

cobre_pd_control = cobre_pd[cobre_pd$Sample_Group=='control',]
cobre_pd_obese = cobre_pd[cobre_pd$Sample_Group=='obese',]
cobre_beta_control = cobre_beta[,colnames(cobre_beta)%in%cobre_pd_control$Sample_Name]
cobre_beta_obese = cobre_beta[,colnames(cobre_beta)%in%cobre_pd_obese$Sample_Name]

#cpg_mean_promotor_tss_limma_sub = readRDS('cpg_mean_promotor_tss_limma_sub.RDS')
#dataExpr = t(cobre_beta)
#dataExpr = t(cpg_geomean_promotor_full_beta)
#dataExpr = t(cpg_mean_promotor_tss_limma_sub) #unadj

library(lumi)
dataExpr = t(beta2m(cpg_mean_promotor_tss_limma_adj)) #adj

#data<-readRDS('/Users/Brainfactory/Downloads/Bing\'s WGCNA code/quant.metabolitemeasurements.rds')
#dataExpr<-t(data[,-1])
#m.mad <- apply(dataExpr,1,mad)
#dataExprVar <- dataExpr[which(m.mad >  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]


#dataExpr <- as.data.frame(t(dataExprVar)) #row=samples, column=cpg/gene/metabolites
#gsg = goodSamplesGenes(dataExpr, verbose = 3)

#dataExpr = t(beta2m(cpg_mean_promotor_tss_limma_sub))
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,networkType=type,verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power = sft$powerEstimate #final_power = 6, adjusted power = 8

#power = 6  #unadjusted mean cpg promotor
#power = 8 adjusted confounding, tss200+1500 
power = 7 #adjusted, mvalue

net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes, #maxBlockSize=nGenes
                       TOMType = type, #minModuleSize = 20,
                       reassignThreshold = 1e-6, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("cobre_4046gene_72sample_adjusted_200_1500",".tom"),
                       verbose = 3)


table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("Module.", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
#traitData_yile <- readRDS('/Users/Brainfactory/Downloads/Bing\'s WGCNA code/clinic.rds')
traitData = cobre_pd
rownames(traitData) = traitData$Sample_Name
traitData = traitData[,-c(1,2,4,5)]

traitData<-as.data.frame(traitData)
levels(traitData$Sample_Group) = c(0,1) #"control" "obese"

traitData$Mat_Ethnicity = as.factor(traitData$Mat_Ethnicity)
levels(traitData$Mat_Ethnicity) = c(0,1,2) #"A"    "C"    "NHPI"

traitData$Gravidity = as.factor(traitData$Gravidity)
levels(traitData$Gravidity) = c(1,2,3,4,5,6) #"1"    "2"    "3"    "4"    "5"    "More"

traitData$Parity = as.factor(traitData$Parity)
levels(traitData$Parity) = c(0,1,2,3) #"0"    "1"    "2"    "More"

traitData$Pat_Ethnicity = as.factor(traitData$Pat_Ethnicity)
levels(traitData$Pat_Ethnicity) = c(0,1,2) #"A"    "C"    "NHPI"

traitData$Sex = as.factor(traitData$Sex)
levels(traitData$Sex) = c(0,1) #"F" "M"

obese<-traitData$Sample_Group
MEs_col = orderMEs(MEs_col)
MEs_colpheno = orderMEs(cbind(MEs_col, obese))
pdf("COBRE_yile_Module-trait-heatmap_adjusted200_1500_merge0.15_mval.pdf", width=8, height=8)
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(6,7,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
#plotTOM = dissTOM^7
#diag(plotTOM) = NA
#TOMplot(plotTOM, net$dendrograms, moduleColors, main = "Network heatmap plot, all genes")
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
quantile(TOM)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("COBRE_yile_WGCNA_0.03_adj_200_1500_merge0.15", ".edges.txt", sep=""),
                               nodeFile = paste("COBRE_yile_WGCNA_0.03_adj_200_1500_merge0.15", ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.03,
                               nodeNames = probes, nodeAttr = moduleColors)


if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf("COBRE_yile_WGCNA_0.03_adj_200_1500_merge0.15_mval.pdf", width=8, height=8)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

