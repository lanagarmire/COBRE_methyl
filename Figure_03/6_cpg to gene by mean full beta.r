library(limma)

cobre_pd = readRDS('cobre_pd.rds')
cobre_beta = readRDS('cobre_beta.rds')


anno_450k = readRDS('FullAnnot_450k.rds') 
colnames(anno_450k) = c('CpG','Island','Gene','Group')
anno_450k = anno_450k[complete.cases(anno_450k),]
#anno_450k_promotor = subset(anno_450k,Group%in%c('TSS1500','TSS200'))

cobre_beta_sub = cobre_beta[rownames(cobre_beta)[rownames(cobre_beta)%in%rownames(anno_450k)],]
cobre_beta_sub_mval = beta2m(cobre_beta_sub)
#cpg_mean_promotor_full_beta = NULL
#cpg_mean_promotor_full_beta_mval = NULL
cpg_geomean_promotor_full_beta = NULL

for(gene_name in unique(anno_450k$Gene)){
  cpg = anno_450k[anno_450k$Gene==gene_name,]$CpG
  cpg_df = NULL
  if(sum(rownames(cobre_beta_sub)%in%cpg)!=0){
    cpg_sub = cobre_beta_sub[rownames(cobre_beta_sub)%in%cpg,]
    if(length(cpg_sub)==72){
      cpg_df = as.data.frame(t(cpg_sub))
      rownames(cpg_df) = c(gene_name)
    }else{
      #print(gene_name)
      cpg_df = as.data.frame(t(apply(cpg_sub, 2, function(x) exp(mean(log(x)))))) #geometric mean
      #cpg_df = as.data.frame(t( colMeans(cpg_sub)))
      rownames(cpg_df) = c(gene_name)
    }
  }
  cpg_geomean_promotor_full_beta = rbind(cpg_geomean_promotor_full_beta,cpg_df)
}

saveRDS(cpg_geomean_promotor_full_beta,'cpg_geomean_promotor_full_beta.RDS')

saveRDS(cpg_mean_promotor_full_beta_mval,'cpg_mean_promotor_full_beta_mval.RDS')

saveRDS(cpg_mean_promotor_full_beta,'cpg_mean_promotor_full_beta.RDS')

cpg_mean_promotor_full_beta = readRDS('cpg_mean_promotor_full_beta.RDS')

design = model.matrix(~ Sample_Group + Sex + Net_Weight_Gain + 
                        Mat_Age + Mat_Ethnicity + Pat_Ethnicity + Gravidity + Gestational_Age, 
                      data = cobre_pd)
options(connectionObserver = NULL)
library(limma)
library(lumi)
fit = lmFit(beta2m(cpg_mean_promotor_full_beta), design)
fit = eBayes(fit)

allg.limma <- topTable(fit, coef=2, n=dim(fit)[1])
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05) #344 genes
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05) #18759

both_pos <- sigg.limma[c('logFC', 'adj.P.Val')]

write.csv(sigg.limma,'sigg.limma_gene_fullbeta_aggregated.csv')

anno_450k_promotor = subset(anno_450k,Group%in%c('TSS1500','TSS200'))
sigg.limma.gene_full_promotor = sigg.limma[rownames(sigg.limma)[rownames(sigg.limma)%in%anno_450k_promotor$Gene],] #323

sigg.limma.gene_full_promotor_order = sigg.limma.gene_full_promotor[order(sigg.limma.gene_full_promotor$logFC,decreasing=TRUE),]
write.csv(sigg.limma.gene_full_promotor_order,'sigg.limma.gene_full_promotor_order.csv')

#both_pos$logFC <- -both_pos$logFC
logfccut=0
both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)

both_pos$Type <- 'NonSig'
both_pos$Type[both_pos$logFC > logfccut] <- 'Hyper'
both_pos$Type[both_pos$logFC < -logfccut] <- 'Hypo'

hypernum <- nrow(subset(both_pos, Type == 'Hyper'))
hyponum <- nrow(subset(both_pos, Type == 'Hypo'))

both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
both_pos <- subset(both_pos, Type != 'NonSig')

library(ggplot2)
library(RColorBrewer)

if(nrow(both_pos) > 50){
  myColor <- densCols(both_pos$logFC, -log10(both_pos$adj.P.Val), 
                      colramp = colorRampPalette(hcl.colors(10, palette = "PuOr",rev = F)))
                      #colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
}else{
  myColor <- rep('blue', nrow(both_pos))
}


both_pos$denscolor <- myColor
rgbmat <- t(col2rgb(myColor))
rgbmat <- as.data.frame(rgbmat, stringsAsFactors = FALSE)
both_pos <- cbind(both_pos, rgbmat)
both_pos <- both_pos[order(-both_pos$blue, both_pos$red, both_pos$green),]
both_pos1 <- subset(both_pos, blue >= red)
both_pos2 <- subset(both_pos, blue < red)
both_pos2 <- both_pos2[order(-both_pos2$blue, both_pos2$red, -both_pos2$green),]
both_pos <- rbind(both_pos1, both_pos2)

both_nonsig <- nonsigg.limma[c('logFC', 'adj.P.Val')]
both_nonsig$Type <- 'NonSig'
both_nonsig <- rbind(both_nonsig, both_pos_nonsig)
both_nonsig$denscolor <- '#C0C0C0'
both_nonsig$red <- 192
both_nonsig$green <- 192
both_nonsig$blue <- 192
nonsignum <- nrow(both_nonsig)

both_pos <- rbind(both_nonsig, both_pos)
colorpart <- subset(both_pos, Type != 'NonSig')
graypart <- subset(both_pos, Type == 'NonSig')
both_pos <- rbind(graypart, colorpart)

hypernum <- sum(both_pos$Type == 'Hyper')
hyponum <- sum(both_pos$Type == 'Hypo')
nonsignum <- sum(both_pos$Type == 'NonSig')

p <- ggplot(both_pos, aes(x=logFC, y=-log10(adj.P.Val)))

xtitle <- 'log2FC'
if(useM == TRUE){
  xtitle <- 'log2FC'
}else{
  xtitle <- 'Difference'
}

p + geom_point(color = both_pos$denscolor, position='jitter') + 
  xlab(xtitle) + ylab('-log10(adj.P.Val)') + 
  ggtitle(paste0('Aggregated Differential Genes (Control vs Obesed)'), 
          subtitle = paste0('COBRE', ' (Hyper sites = ', hypernum, ', Hypo sites= ', hyponum, 
                            ', NonSig = ', nonsignum, ')')) + 
  geom_hline(yintercept = -log10(0.05), color = 'red') + 
  geom_vline(xintercept = logfccut, color = 'blue') + 
  geom_vline(xintercept = -logfccut, color = 'blue') + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  theme_bw() + 
  theme(panel.grid = element_blank())


both_pos$gene = rownames(both_pos)
gseaprobes <- both_pos[c('gene', 'adj.P.Val', 'logFC', 'Type')]
manhattantable <- makeManhattantable()
plotManhattan()
