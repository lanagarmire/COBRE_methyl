options(connectionObserver = NULL)
library(limma)
library(lumi)
library(ggplot2)
library(RColorBrewer)

myCombat = readRDS('cobre_combat_final.rds')
cobre_pd = readRDS('cobre_pd.rds')
cobre_sov_res = load('cobre_final_res.Rdata')


cobre_pd$Sample_Group = relevel(factor(cobre_pd$Sample_Group), ref="control")
cobre_pd$Mat_Ethnicity = relevel(factor(cobre_pd$Mat_Ethnicity), ref="A")
cobre_pd$Pat_Ethnicity = relevel(factor(cobre_pd$Pat_Ethnicity), ref="A")

design = model.matrix(~ Sample_Group + Sex + Net_Weight_Gain + 
                        Mat_Age + Mat_Ethnicity + Pat_Ethnicity + Gravidity + Gestational_Age, 
                      data = cobre_pd)

fit = lmFit(beta2m(myCombat), design)
fit = eBayes(fit)

allg.limma <- topTable(fit, coef=2, n=dim(fit)[1])
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)

promoteronly <- TRUE
TSS200only <- TRUE
useM <- TRUE
type <- 'All'
logfccut <- 0
tag <- 'Obesity to Control'
gestwkdmr <- read.table('geo_bumpres_v5_1000_pval0.05.txt', 
                        sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#Draw probe valcano
both_pos <- sigg.limma[c('logFC', 'adj.P.Val')]
both_pos <- as.data.frame(both_pos, stringsAsFactors = FALSE)
both_pos$Type <- 'NonSig'
both_pos$Type[both_pos$logFC > logfccut] <- 'Hyper'
both_pos$Type[both_pos$logFC < -logfccut] <- 'Hypo'

hypernum <- nrow(subset(both_pos, Type == 'Hyper'))
hyponum <- nrow(subset(both_pos, Type == 'Hypo'))

both_pos_nonsig <- subset(both_pos, Type == 'NonSig')
both_pos <- subset(both_pos, Type != 'NonSig')


if(nrow(both_pos) > 50){
  myColor <- densCols(both_pos$logFC, -log10(both_pos$adj.P.Val), 
                      colramp = colorRampPalette(rev(rainbow(10, end=4/6))))
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



gestwkproberatio <- function(platform = '450k', totalprobes = both_pos, gestwkreg = gestwkdmr, 
                             pvalgradient = TRUE, pvalend = 0.1, pvalcutoff = 0.05, 
                             adjusttype = type){

  if(platform == 'EPIC'){
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    chrinfo <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
  }else if(platform == '450k'){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    chrinfo <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  }else if(platform == '27k'){
    library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
    chrinfo <- IlluminaHumanMethylation27kanno.ilmn12.hg19::Locations
  }
  
  gestwkdmrir <- IRanges(start = gestwkreg$start + 1, end = gestwkreg$end)
  gestwkdmrgr <- GRanges(seqnames = gestwkreg$seqnames, ranges = gestwkdmrir, strand = gestwkreg$strand)
  gestwkdmrgr <- reduce(gestwkdmrgr)
  
  
  overnums <- c()
  
  if(pvalgradient == TRUE){
    
    pvalseq <- seq(0, pvalend, 0.001)
    
  }else{
    pvalseq <- pvalcutoff
  }
  
  for(i in pvalseq){
    
    diffprobes <- subset(totalprobes, adj.P.Val < i)
    
    diffprobes <- row.names(diffprobes)
    diffprobes <- intersect(diffprobes, row.names(chrinfo))
    diffprobes <- chrinfo[diffprobes,]
    
    diffprobesir <- IRanges(start = diffprobes$pos, width = 1)
    diffprobesgr <- GRanges(seqnames = diffprobes$chr, ranges = diffprobesir, strand = diffprobes$strand)
    
    
    
    intprobes <- findOverlaps(query = diffprobesgr, subject = gestwkdmrgr, type = c('within'))
    
    overprobes <- diffprobes[intprobes@from,]
    overdmrs <- gestwkreg[intprobes@to,]
    
    unioverprobes <- unique(overprobes)
    unioverprobes <- totalprobes[row.names(unioverprobes),]
    unioverdmrs <- unique(overdmrs)
    
    overnum <- nrow(unioverprobes)
    
    if(i == pvalcutoff){
      
      platprobes <- IRanges(start = chrinfo$pos, width = 1)
      platprobes <- GRanges(seqnames = chrinfo$chr, ranges = platprobes, strand = chrinfo$strand)
      gestwkprobes <- findOverlaps(query = platprobes, subject = gestwkdmrgr, type = c('within'))
      
      res <- list()
      
      res$gestprobes <- unioverprobes
      res$diffprobes <- diffprobes
      
      
    }
    
    overnums <- c(overnums, overnum)
    
    print(i)
    
    
  }
  
  gradientres <- data.frame(logfc = pvalseq, absnums = overnums, stringsAsFactors = FALSE)
  
  res$gradientres <- gradientres
  
  return(res)
  
}

gestwkprobe450k <- gestwkproberatio(platform = '450k',totalprobes = both_pos)  #20 probes 15648


if(nrow(gestwkprobe450k$gestprobes) > 0){
  both_pos[row.names(gestwkprobe450k$gestprobes),]$Type <- 'NonSig'
  both_pos[row.names(gestwkprobe450k$gestprobes),]$denscolor <- '#C0C0C0'
  both_pos[row.names(gestwkprobe450k$gestprobes),]$red <- 192
  both_pos[row.names(gestwkprobe450k$gestprobes),]$green <- 192
  both_pos[row.names(gestwkprobe450k$gestprobes),]$blue <- 192
}

colorpart <- subset(both_pos, Type != 'NonSig')
graypart <- subset(both_pos, Type == 'NonSig')
both_pos <- rbind(graypart, colorpart)

hypernum <- sum(both_pos$Type == 'Hyper')
hyponum <- sum(both_pos$Type == 'Hypo')
nonsignum <- sum(both_pos$Type == 'NonSig')

p <- ggplot(both_pos, aes(x=logFC, y=-log10(adj.P.Val)))

if(useM == TRUE){
  xtitle <- 'log2FC'
}else{
  xtitle <- 'Difference'
}

p + geom_point(color = both_pos$denscolor, position='jitter') + 
  xlab(xtitle) + ylab('-log10(adj.P.Val)') + 
  ggtitle(paste0('Differential Methylation Probes (Control vs Obesed)'), 
          subtitle = paste0('COBRE adj GestWk', ' (Hyper = ', hypernum, ', Hypo = ', hyponum, 
                            ', NonSig = ', nonsignum, ')')) + 
  geom_hline(yintercept = -log10(0.05), color = 'red') + 
  geom_vline(xintercept = logfccut, color = 'blue') + 
  geom_vline(xintercept = -logfccut, color = 'blue') + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  theme_bw() + 
  theme(panel.grid = element_blank())



cobre_limma_sig = sigg.limma[!rownames(sigg.limma)%in%rownames(gestwkprobe450k$gestprobes),]

saveRDS(cobre_limma_sig,'cobre_limma_sig_1w5.RDS')

saveRDS(sigg.limma,file='sigg.limma_1w5.rds')

table(cobre_limma_sig[cobre_limma_sig$logFC< (-0.1),]$Type)

cobre_limma_sig = sigg.limma
diffs_cpg = rownames(cobre_limma_sig)
rownames(cobre_limma_sig)%in%rownames(sigg.limma)


cobre_limma_sig$Type = ifelse(cobre_limma_sig$logFC>0, 'Hyper', 'Hypo')

anno_450k = readRDS('anno_450k.rds')
anno_450k = as.data.frame(anno_450k)
colnames(anno_450k) = c('CpG','Gene','Group')

anno_450k = readRDS('FullAnnot_450k.rds') 
colnames(anno_450k) = c('CpG','Island','Gene','Group')
anno_450k_promotor_TSS200 = anno_450k[anno_450k$Group %in% c('TSS200'),]
diffs_cpg_promotor_TSS200 = diffs_cpg[diffs_cpg%in%anno_450k_promotor_TSS200$CpG] #1w5: 2347 9k:1489

anno_450k_promotor_TSS1500 = anno_450k[anno_450k$Group %in% c('TSS1500'),]
diffs_cpg_promotor_TSS1500  = diffs_cpg[diffs_cpg%in%anno_450k_promotor_TSS1500$CpG]

diffs_cpg_df = cbind(rownames(cobre_limma_sig),cobre_limma_sig) 
colnames(diffs_cpg_df)[1] = c('CpG')
diffs_cpg = diffs_cpg_df$CpG
diff_gene = anno_450k[anno_450k$CpG%in%diffs_cpg,] #15648

diff_cpg_gene = merge(anno_450k,diffs_cpg_df,by='CpG')
diff_cpg_gene = diff_cpg_gene[complete.cases(diff_cpg_gene),]

dim(diff_cpg_gene[diff_cpg_gene$Type=='Hypo'&diff_cpg_gene$Group=='TSS200',])
length(unique(diff_cpg_gene[diff_cpg_gene$Type=='Hypo'&diff_cpg_gene$Group=='TSS200',]$Gene))
nrow(diff_cpg_gene[diff_cpg_gene$Group=='TSS200',])
diff_cpg_gene = diff_cpg_gene[complete.cases(diff_cpg_gene),]

table(diff_cpg_gene[complete.cases(diff_cpg_gene),]$Group)

saveRDS(diff_cpg_gene,'COBRE_diff_cpg_gene1w5.rds')

BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)


diff_cpg_gene = readRDS('COBRE_diff_cpg_gene1w5.rds')

diff_cpg_gene_hyper = diff_cpg_gene[diff_cpg_gene$Type=='Hyper',]
diff_cpg_gene_hyper[order(diff_cpg_gene_hyper$logFC,decreasing=TRUE),][1:20,]
diff_cpg_gene_hypo = diff_cpg_gene[diff_cpg_gene$Type=='Hypo',]
diff_cpg_gene_hypo[order(diff_cpg_gene_hypo$logFC),][1:20,]

diff_gene = diff_cpg_gene[complete.cases(diff_cpg_gene),]
diff_gene_tss200 = diff_cpg_gene[diff_cpg_gene$Group=='TSS200',]
diff_gene_tss200_tss1500 = diff_cpg_gene[diff_cpg_gene$Group%in%c('TSS200','TSS1500'),]

dim(unique(diff_gene_tss200[diff_gene_tss200$Type=='Hypo',])) #289
length(unique(diff_gene_tss200[diff_gene_tss200$Type=='Hypo',]$Gene)) #275

dim(unique(diff_gene_tss200[diff_gene_tss200$Type=='Hyper',])) #1200
length(unique(diff_gene_tss200[diff_gene_tss200$Type=='Hyper',]$Gene)) #1112

cobre_limma_sig_tss200 = cobre_limma_sig[rownames(cobre_limma_sig)%in%diff_gene_tss200$CpG,c('logFC','adj.P.Val')]
cobre_limma_sig_tss200_tss1500 = cobre_limma_sig[rownames(cobre_limma_sig)%in%diff_gene_tss200_tss1500$CpG,c('logFC','adj.P.Val')]

