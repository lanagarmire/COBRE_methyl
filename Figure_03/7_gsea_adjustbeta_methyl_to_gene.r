#calculate cell-type adjusted beta
library(limma)

cobre_pd = readRDS('cobre_pd.rds')
myCombat = readRDS('cobre_combat_final.rds')
myCombat = readRDS('cobre_beta.rds')

cobre_mycombat_cap = myCombat
cobre_mycombat_cap[cobre_mycombat_cap < 0.001]=0.001
cobre_mycombat_cap[cobre_mycombat_cap > 0.999]=0.999

cobre_mycombat_cap_mvalue = beta2m(cobre_mycombat_cap)


design2 = model.matrix(~ Sex + Net_Weight_Gain + 
                         Mat_Age + Mat_Ethnicity + Pat_Ethnicity + Gravidity + Gestational_Age+Parity+Hemoglobin, 
                       data = cobre_pd)

fit_nogroup = lmFit(cobre_mycombat_cap_mvalue, design2)
fit_nogroup = eBayes(fit_nogroup)
residuals = residuals(fit_nogroup,cobre_mycombat_cap_mvalue)

adj.m <-residuals+matrix(apply(cobre_mycombat_cap_mvalue, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

adj.betas <- m2beta(adj.m)
saveRDS(adj.betas,'cobre_final_adj_betas_allparityhemo.rds')
saveRDS(cobre_pd,'cobre_pd_final.rds')
#design3 = model.matrix(~ Sample_Group, data = cobre_pd)
#fit_onlygroup = lmFit(adj.m, design3)
#fit_onlygroup = eBayes(fit_onlygroup)
#allg.limma.adj <- topTable(fit_onlygroup, coef=2, n=dim(fit1)[1])
#sigg.limma.adj <- subset(allg.limma.adj, P.Value < 0.05)



anno_450k = readRDS('/Users/Brainfactory/FullAnnot_450k.rds') 
colnames(anno_450k) = c('CpG','Island','Gene','Group')
anno_450k = anno_450k[complete.cases(anno_450k),]
anno_450k_promotor = subset(anno_450k,Group%in%c('TSS1500','TSS200'))

adj.betas_with_anno = adj.betas[rownames(adj.betas)%in%anno_450k_promotor$CpG,]

cpg_mean_promotor = NULL

for(gene_name in unique(anno_450k_promotor$Gene)){
  cpg = anno_450k_promotor[anno_450k_promotor$Gene==gene_name,]$CpG
  cpg_df = NULL
  if(sum(rownames(adj.betas_with_anno)%in%cpg)!=0){
    cpg_sub = adj.betas_with_anno[rownames(adj.betas_with_anno)%in%cpg,]
    if(length(cpg_sub)==72){
      cpg_df = as.data.frame(t(cpg_sub))
      rownames(cpg_df) = c(gene_name)
    }else{
      cpg_df = as.data.frame(t(colMeans(cpg_sub)))
      rownames(cpg_df) = c(gene_name)
    }
  }
  cpg_mean_promotor = rbind(cpg_mean_promotor,cpg_df)
}

saveRDS(cpg_mean_promotor,'cpg_mean_promotor2_usingcobrebetards.rds')
saveRDS(cpg_mean_promotor22,'cpg_mean_promotor.rds')

Description = rep(NA,dim(cpg_mean_promotor)[1])
cpg_mean_promotor_full = cbind(rownames(cpg_mean_promotor),Description,(cpg_mean_promotor))
colnames(cpg_mean_promotor_full)[1] = c('Name')
colnames(cpg_mean_promotor_full)[2] = c('Description')
write.table(cpg_mean_promotor_full,file='cpg_mean_promotor.txt',quote = FALSE,na='NA',sep='\t',row.names = FALSE)


cpg_mean_promotor_full_mval = cbind(rownames(cpg_mean_promotor),Description,beta2m(cpg_mean_promotor))
colnames(cpg_mean_promotor_full_mval)[1] = c('Name')
colnames(cpg_mean_promotor_full_mval)[2] = c('Description')
write.table(cpg_mean_promotor_full_mval,file='cpg_mean_promotor_mval.txt',quote = FALSE,na='NA',sep='\t',row.names = FALSE)
