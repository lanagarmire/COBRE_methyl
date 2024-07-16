#differential
library(DESeq2)
library(edgeR)
library(limma)
library(lumi)
library(ggplot2)
#cobre_bulk_pd_sub,cobre_fc_expr_sub,file='cobre_bulk_fc_pd_sub.rdata
load('cobre_bulk_fc_pd_sub.rdata',verbose=TRUE)
counts = as.matrix(cobre_fc_expr_sub)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

mm = model.matrix(~ Sample_Group + Sex + Net_Weight_Gain + 
                             Mat_Age + Mat_Ethnicity + Pat_Ethnicity + Gravidity + Gestational_Age, 
                           data = cobre_bulk_pd_sub)

y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
#ontr <- makeContrasts( labelsObese - labelsNormal , levels = colnames(coef(fit)))
#tmp <- contrasts.fit(fit, contr)
#tmp <- eBayes(tmp)
tmp <- eBayes(fit)
#top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table = topTable(tmp, coef=2, n=dim(fit)[1])
cobre_bulk_diff_gene_toptable = top.table

saveRDS(cobre_bulk_diff_gene_toptable,'cobre_bulk_diff_gene_toptable.RDS')
cobre_bulk_diff_gene_toptable = readRDS('/Users/Brainfactory/prelim_data/cobre_bulk_diff_gene_toptable.RDS')
top.table = cobre_bulk_diff_gene_toptable
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))



cobre_bulk_diff_gene_toptable
common_diff_gene_hyper = rownames(cobre_bulk_diff_gene_toptable)[rownames(cobre_bulk_diff_gene_toptable)%in%hyper_island200_1500_sort_highestfc_gene$Gene]
common_diff_gene_hypo = rownames(cobre_bulk_diff_gene_toptable)[rownames(cobre_bulk_diff_gene_toptable)%in%hypo_island200_1500_sort_highestfc_gene$Gene]

#868
diff_hyper = cobre_bulk_diff_gene_toptable[rownames(cobre_bulk_diff_gene_toptable)%in%hyper_island200_1500_sort_highestfc_gene$Gene,]
diff_hyper
diff_hyper_cpg = hyper_island200_1500_sort_highestfc_gene[hyper_island200_1500_sort_highestfc_gene$Gene%in%common_diff_gene_hyper,]
rownames(diff_hyper_cpg) = diff_hyper_cpg$Gene
diff_hyper_cpg = diff_hyper_cpg[rownames(diff_hyper),]

#plot(diff_hyper_cpg$logFC,diff_hyper$logFC)
df_hyper_logfc =  as.data.frame(cbind(diff_hyper_cpg$logFC,diff_hyper$logFC))
colnames(df_hyper_logfc) = c('CpG','Gene')
rownames(df_hyper_logfc) = rownames(diff_hyper_cpg)
ggplot(df_hyper_logfc, aes(x=CpG, y=Gene))+
  geom_point( )+
  geom_smooth(method=lm)+
  ggtitle('logFC comparison between hypermethylated sites and genes')+
  xlab('Hypermethylatde CpG sites logFC') +
  ylab('Corresponding genes logFC')

#261
diff_hypo = cobre_bulk_diff_gene_toptable[rownames(cobre_bulk_diff_gene_toptable)%in%hypo_island200_1500_sort_highestfc_gene$Gene,]
diff_hypo_cpg = hypo_island200_1500_sort_highestfc_gene[hypo_island200_1500_sort_highestfc_gene$Gene%in%common_diff_gene_hypo,]
rownames(diff_hypo_cpg) = diff_hypo_cpg$Gene
diff_hypo_cpg = diff_hypo_cpg[rownames(diff_hypo),]
#plot(diff_hypo_cpg$logFC,diff_hypo$logFC)

df_hypo_logfc =  as.data.frame(cbind(diff_hypo_cpg$logFC,diff_hypo$logFC))
colnames(df_hypo_logfc) = c('CpG','Gene')
rownames(df_hypo_logfc) = rownames(diff_hypo_cpg)
ggplot(df_hypo_logfc, aes(x=CpG, y=Gene))+
  geom_point( )+
  geom_smooth(method=lm)+
  ggtitle('logFC comparison between hypomethylated sites and genes')+
  xlab('Hypomethylatde CpG sites logFC') +
  ylab('Corresponding genes logFC')




cobre_pd = readRDS('/Users/Brainfactory/cobre_pd.rds')
cobre_beta = readRDS('/Users/Brainfactory/cobre_beta.rds')

cobre_beta_sig_gene_sub = cobre_beta[rownames(cobre_beta)%in%diff_gene_tss200_tss1500$CpG,colnames(cobre_beta)%in%rownames(cobre_bulk_pd_sub)]
cobre_beta_sig_gene_sub_mval = beta2m(cobre_beta_sig_gene_sub)

load('/Users/Brainfactory/cobre_bulk_fc_pd_sub.rdata') #cobre_bulk_pd_sub,cobre_fc_expr_sub,
cobre_bulk_pd_sub
cobre_fc_expr_sub = cobre_fc_expr_sub[rowSums(cobre_fc_expr_sub)!=0,]

hyper_gene = diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$Type=='Hyper'&diff_gene_tss200_tss1500$Island=='Island',]
hypo_gene = diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$Type=='Hypo'&diff_gene_tss200_tss1500$Island=='Island',]

cobre_bulk_pd_sub_sorted = cobre_bulk_pd_sub[order(cobre_bulk_pd_sub$Sample_Group),]
hyper_sorted = hyper_gene[order(-hyper_gene$logFC),]
hypo_sorted = hypo_gene[order(hypo_gene$logFC),]

hyper_island200_1500_sort_highestfc_gene = hyper_sorted[!duplicated(hyper_sorted$Gene),]
hypo_island200_1500_sort_highestfc_gene = hypo_sorted[!duplicated(hypo_sorted$Gene),]

cobre_hyper_sub = cobre_beta_sig_gene_sub_mval[rownames(cobre_beta_sig_gene_sub_mval)%in%hyper_island200_1500_sort_highestfc_gene$CpG,]
cobre_hypo_sub = cobre_beta_sig_gene_sub_mval[rownames(cobre_beta_sig_gene_sub_mval)%in%hypo_island200_1500_sort_highestfc_gene$CpG,]

#959 common hyper genes, final 951 without 0 genes, 1443 from 200+1500
common_hyper_gene = rownames(cobre_fc_expr_sub)[rownames(cobre_fc_expr_sub)%in%hyper_island200_1500_sort_highestfc_gene$Gene]
common_hyper_gene_cpg = hyper_island200_1500_sort_highestfc_gene[hyper_island200_1500_sort_highestfc_gene$Gene%in%common_hyper_gene,]
rownames(common_hyper_gene_cpg) = common_hyper_gene_cpg$CpG

#321 common hypo genes, final 312 without 0 genes, 572 from 200+1500
common_hypo_gene = rownames(cobre_fc_expr_sub)[rownames(cobre_fc_expr_sub)%in%hypo_island200_1500_sort_highestfc_gene$Gene]
common_hypo_gene_cpg = hypo_island200_1500_sort_highestfc_gene[hypo_island200_1500_sort_highestfc_gene$Gene%in%common_hypo_gene,]
rownames(common_hypo_gene_cpg) = common_hypo_gene_cpg$CpG

cobre_fc_expr_sub_hyper = cobre_fc_expr_sub[rownames(cobre_fc_expr_sub)%in%hyper_island200_1500_sort_highestfc_gene$Gene,]
cobre_fc_expr_sub_hypo = cobre_fc_expr_sub[rownames(cobre_fc_expr_sub)%in%hypo_island200_1500_sort_highestfc_gene$Gene,]

cobre_fc_expr_sub_hyper_log = log2(cobre_fc_expr_sub_hyper+1)
cobre_fc_expr_sub_hypo_log = log2(cobre_fc_expr_sub_hypo+1)

cobre_hyper_sub_m = cobre_beta_sig_gene_sub_mval[rownames(cobre_beta_sig_gene_sub_mval)%in%common_hyper_gene_cpg$CpG,]
cobre_hypo_sub_m = cobre_beta_sig_gene_sub_mval[rownames(cobre_beta_sig_gene_sub_mval)%in%common_hypo_gene_cpg$CpG,]

testaa = merge(common_hyper_gene_cpg,cobre_hyper_sub_m,by=0, all=TRUE) 
testbb = merge(common_hypo_gene_cpg,cobre_hypo_sub_m,by=0, all=TRUE) 

cobre_hyper_sub_m_cpg_genename = testaa
rownames(cobre_hyper_sub_m_cpg_genename) = cobre_hyper_sub_m_cpg_genename$Gene
#cobre_hyper_sub_m_cpg_genename = cobre_hyper_sub_m_cpg_genename[,-c(1:12)]

cobre_hypo_sub_m_cpg_genename = testbb
rownames(cobre_hypo_sub_m_cpg_genename) = cobre_hypo_sub_m_cpg_genename$Gene
#cobre_hypo_sub_m_cpg_genename = cobre_hypo_sub_m_cpg_genename[,-c(1:12)]

#cobre_hyper_sub_beta = m2beta(cobre_hyper_sub)
#cobre_hypo_sub_beta = m2beta(cobre_hypo_sub)

#For correlation:
library(pheatmap)
cobre_hyper_sub_m_cpg_genename = cobre_hyper_sub_m_cpg_genename[rownames(cobre_fc_expr_sub_hyper_log),rownames(cobre_bulk_pd_sub_sorted)]
cobre_fc_expr_sub_hyper_log = cobre_fc_expr_sub_hyper_log[,rownames(cobre_bulk_pd_sub_sorted)]
cor_hyper = cor(t(cobre_hyper_sub_m_cpg_genename),t(cobre_fc_expr_sub_hyper_log))
#cor_hyper_complete = cor_hyper[complete.cases(cor_hyper),]
quantile(diag(cor_hyper)[complete.cases(diag(cor_hyper))])

diag_cor_hyper = diag(cor_hyper)
names(diag_cor_hyper) = rownames(cor_hyper)

diag_cor_hyper_high = diag_cor_hyper[diag_cor_hyper<(-0.2)] #123 genes, 196 genes from 200+1500
#plot(diff_hyper_cpg$logFC,diff_hyper$logFC)
df_hyper_logfc_high =  merge(diff_hyper_cpg[diff_hyper_cpg$Gene%in%names(diag_cor_hyper_high),],
                             diff_hyper[rownames(diff_hyper)%in%names(diag_cor_hyper_high),],
                             by=0,all=TRUE)
rownames(df_hyper_logfc_high) = df_hyper_logfc_high$Row.names
df_hyper_logfc_high = df_hyper_logfc_high[,c('logFC.x','logFC.y')]
colnames(df_hyper_logfc_high) = c('CpG','Gene')
#rownames(df_hyper_logfc) = rownames(diff_hyper_cpg)
ggplot(df_hyper_logfc_high[!rownames(df_hyper_logfc_high)%in%c('LRCH3','SAFB2','BCAS3'),], aes(x=CpG, y=Gene))+
  geom_point( )+
  geom_smooth(method=lm)+
  ggtitle('logFC comparison between top negatively correlated \nhypermethylated sites and genes')+
  xlab('Hypermethylatde CpG sites logFC') +
  ylab('Corresponding genes logFC')



diff_hyper_cpg[diff_hyper_cpg$Gene%in%names(diag_cor_hyper_high),]$logFC
diff_hyper[rownames(diff_hyper)%in%names(diag_cor_hyper_high),]$logFC

pdf('cor_hyper2.pdf')
pheatmap(cor_hyper,cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()
#colnames(cor_hyper)[colSums(is.na(cor_hyper)) > 0]
#NA: "PAX7"     "RIPPLY2"  "UNCX"     "MIR2110"  "HOXC11"   "HOXC8"    "ANKRD34C" "MIR1181" 

cobre_hyper_sub_m_cpg_genename['PAX7',]
cobre_fc_expr_sub_hyper_log['PAX7',]

#hyper_anno = cobre_bulk_pd_sub_sorted[,c('Sample_Group'), drop=FALSE]
#colnames(hyper_anno) = c('Group')
#pheatmap(cor_hyper, annotation_row = hyper_anno, cluster_rows = F,
#         annotation_col = hyper_anno, cluster_cols = F)

cobre_hypo_sub_m_cpg_genename = cobre_hypo_sub_m_cpg_genename[rownames(cobre_fc_expr_sub_hypo_log),rownames(cobre_bulk_pd_sub_sorted)]
#cobre_fc_expr_sub_hypo_log
cor_hypo = cor(t(cobre_hypo_sub_m_cpg_genename),t(cobre_fc_expr_sub_hypo_log))
quantile(diag(cor_hypo)[complete.cases(diag(cor_hypo))])
diag(cor_hypo)[diag(cor_hypo)<(-0.1)]

hypo_anno = cobre_bulk_pd_sub_sorted[,c('Sample_Group'), drop=FALSE]
colnames(hypo_anno) = c('Group')
pheatmap(cor_hypo, annotation_row = hypo_anno, cluster_rows = F,
         annotation_col = hypo_anno, cluster_cols = F)

#highest anti-correlated cpg genes to select gene and sites for pathway enrichment
#diag(cor_hyper)[complete.cases(diag(cor_hyper))]


#new final hypo part
cobre_hypo_sub_m_cpg_genename = cobre_hypo_sub_m_cpg_genename[rownames(cobre_fc_expr_sub_hypo_log),rownames(cobre_bulk_pd_sub_sorted)]
cobre_fc_expr_sub_hypo_log = cobre_fc_expr_sub_hypo_log[,rownames(cobre_bulk_pd_sub_sorted)]
cor_hypo = cor(t(cobre_hypo_sub_m_cpg_genename),t(cobre_fc_expr_sub_hypo_log))
#cor_hypo_complete = cor_hypo[complete.cases(cor_hypo),]
quantile(diag(cor_hypo)[complete.cases(diag(cor_hypo))])

diag_cor_hypo = diag(cor_hypo)
names(diag_cor_hypo) = rownames(cor_hypo)

diag_cor_hypo_high = diag_cor_hypo[diag_cor_hypo<(-0.1)] #123 genes, 118 200+1500
#plot(diff_hypo_cpg$logFC,diff_hypo$logFC)
df_hypo_logfc_high =  merge(diff_hypo_cpg[diff_hypo_cpg$Gene%in%names(diag_cor_hypo_high),],
                            diff_hypo[rownames(diff_hypo)%in%names(diag_cor_hypo_high),],
                            by=0,all=TRUE)
rownames(df_hypo_logfc_high) = df_hypo_logfc_high$Row.names
df_hypo_logfc_high = df_hypo_logfc_high[,c('logFC.x','logFC.y')]
colnames(df_hypo_logfc_high) = c('CpG','Gene')
#rownames(df_hypo_logfc) = rownames(diff_hypo_cpg)
ggplot(df_hypo_logfc_high[rownames(df_hypo_logfc_high),], aes(x=CpG, y=Gene))+
  geom_point( )+
  geom_smooth(method=lm)+
  ggtitle('logFC comparison between top negatively correlated \nhypomethylated sites and genes')+
  xlab('Hypermethylatde CpG sites logFC') +
  ylab('Corresponding genes logFC')

pdf('cor_hypo.pdf')
pheatmap(cor_hypo,cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()




diag_cor_hyper_high = diag_cor_hyper[diag_cor_hyper<(-0.2)] 

cobre_hyper_sub_m_cpg_genename[1:5,1:5]
diag_cor_hyper_high_cpg_logfc = common_hyper_gene_cpg[common_hyper_gene_cpg$Gene%in%names(diag_cor_hyper_high),]
rownames(diag_cor_hyper_high_cpg_logfc) = diag_cor_hyper_high_cpg_logfc$Gene
hyper_logfc_cor_merged = merge(diag_cor_hyper_high,diag_cor_hyper_high_cpg_logfc,by=0)


ggplot(hyper_logfc_cor_merged, aes(x=2^logFC, y=x))+
  geom_point( )+
  ggtitle('Reproduced logFC comparison between top negatively correlated \nhypermethylation logFC vs. correlation')+
  xlab('Hypermethylated CpG sites logFC') +
  ylab('Correlation')

hyper_logfc_cor_merged_highcorr = hyper_logfc_cor_merged[hyper_logfc_cor_merged$x< (-0.399),]


diag_cor_hypo_high = diag_cor_hypo[diag_cor_hypo<(-0.2)] 

cobre_hypo_sub_m_cpg_genename[1:5,1:5]
diag_cor_hypo_high_cpg_logfc = common_hypo_gene_cpg[common_hypo_gene_cpg$Gene%in%names(diag_cor_hypo_high),]
rownames(diag_cor_hypo_high_cpg_logfc) = diag_cor_hypo_high_cpg_logfc$Gene
hypo_logfc_cor_merged = merge(diag_cor_hypo_high,diag_cor_hypo_high_cpg_logfc,by=0)


ggplot(hypo_logfc_cor_merged, aes(x=-2^(-logFC), y=x))+
  geom_point( )+
  ggtitle('Reproduced logFC comparison between top negatively correlated \nhypomethylation logFC vs. correlation')+
  xlab('Hypomethylated CpG sites logFC') +
  ylab('Correlation')


hyper_logfc_cor_merged$fc = 2^(hyper_logfc_cor_merged$logFC)
hypo_logfc_cor_merged$fc = -2^(-hypo_logfc_cor_merged$logFC)


hyper_gene_relationship = hyper_logfc_cor_merged[,c('Row.names','x','fc')]
hypo_gene_relationship = hypo_logfc_cor_merged[,c('Row.names','x','fc')]
