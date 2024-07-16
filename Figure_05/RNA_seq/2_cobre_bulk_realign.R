cobre_fc = read.table('/Users/Brainfactory/cobre_counts_name.txt',sep='\t',header = T)
cobre_fc_expr = cobre_fc[,-c(2,3,4,5,6)]
rownames(cobre_fc_expr) = cobre_fc_expr$Geneid
#cobre_fc_expr = cobre_fc_expr[,-c(1)]
colnames(cobre_fc_expr) = gsub('X.home.yhdu.cobre_rna.cobre_align_final.','',colnames(cobre_fc_expr))
colnames(cobre_fc_expr) = gsub('.Aligned.sortedByCoord.out.bam','',colnames(cobre_fc_expr))
colnames(cobre_fc_expr) = gsub('Sample_','MT',colnames(cobre_fc_expr))
cobre_bulk_pd = cobre_pd
rownames(cobre_bulk_pd) = cobre_bulk_pd$Sample_Name
#cobre_bulk_pd = cobre_bulk_pd[,-c(1)]
save(cobre_bulk_pd,cobre_bulk_pd,file='cobre_bulk_fc_pd.rdata')

cobre_bulk_pd_sub = cobre_bulk_pd[rownames(cobre_bulk_pd)%in%colnames(cobre_fc_expr),]
cobre_fc_expr_sub = cobre_fc_expr[,colnames(cobre_fc_expr)%in%rownames(cobre_bulk_pd)]
save(cobre_bulk_pd_sub,cobre_fc_expr_sub,file='cobre_bulk_fc_pd_sub.rdata')

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = cobre_fc_expr_sub,
                              colData = cobre_bulk_pd_sub,
                              design = ~ Sample_Group)
dds <- DESeq(dds)
res <- results(dds)
res = res[complete.cases(res),]
res[res$padj<0.05,]


#limma-voom
library(edgeR)

labels = dds$Sample_Group
mm <- model.matrix(~0 + labels)
y <- voom(assay(dds), mm, plot = T)

fit <- lmFit(y, mm)
contr <- makeContrasts( labelsobese - labelscontrol , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))


design = model.matrix(~ Sample_Group + Sex + Net_Weight_Gain + 
                        Mat_Age + Mat_Ethnicity + Pat_Ethnicity + Gravidity + Gestational_Age, 
                      data = cobre_bulk_pd_sub)
y <- voom(assay(dds), design, plot = T)
fit = lmFit(y, design)
fit = eBayes(fit)
allg.limma <- topTable(fit, coef=2, n=dim(fit)[1])
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)


library("pheatmap")
library("RColorBrewer")

barplot(round( colSums(assay(dds)), 1 ),col=cobre_bulk_pd_sub$Sample_Group)

sampleDists <- dist(t(assay(dds)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( cobre_bulk_pd_sub$Sample_Group,rownames(cobre_bulk_pd_sub), sep = " - " ) 
colnames(sampleDistMatrix) <- rownames(cobre_bulk_pd_sub)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix,cluster_rows=FALSE, cluster_cols=FALSE,
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, col = colors,fontsize_col=7)





library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

cobre_cpg_mean_promotor = readRDS("/Users/Brainfactory/Downloads/GSE72094_RAW/cpg_mean_promotor.rds")

sum(rownames(cobre_cpg_mean_promotor)%in%rownames(assay(dds))) #15535 genes in common

common_id = rownames(cobre_cpg_mean_promotor)[rownames(cobre_cpg_mean_promotor)%in%rownames(assay(dds))]


#rownames(rld2) = bulk_id[match(names(bulk_id), rownames(rld2))]

common_bulk = assay(dds)[rownames(assay(dds))%in%common_id,]
common_methyl = cobre_cpg_mean_promotor[rownames(cobre_cpg_mean_promotor)%in%common_id,]


#mean_common_bulk = aggregate(common_bulk, list(row.names(common_bulk)), mean)
#rownames(mean_common_bulk) = mean_common_bulk$Group.1
#mean_common_bulk = mean_common_bulk[,-c(1)]

common_methyl <- common_methyl[rownames(common_bulk),]
common_methyl <- common_methyl[,colnames(common_bulk)]

library(lumi)
gene_pearson = diag(cor(t(common_bulk),t(beta2m(common_methyl))))
hist(gene_pearson,main='Gene level pearson correlation distribution')


diff_cobre_gene_list = unique(diff_cpg_gene$Gene)

diff_bulk_df = common_bulk[rownames(common_bulk) %in%diff_cobre_gene_list,]

write.table(diff_bulk_df,file='diff_bulk_df.txt',quote=FALSE)
