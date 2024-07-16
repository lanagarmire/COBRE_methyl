#missMethyl version 1.28.0
BiocManager::install("missMethyl")

options(connectionObserver = NULL)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

diff_cpg_gene = readRDS('/Users/Brainfactory/COBRE_diff_cpg_gene1w5.rds')
hypo_genes = diff_cpg_gene[diff_cpg_gene$Type=='Hypo'&diff_cpg_gene$Group%in%c('TSS200','TSS1500'),c('CpG','Gene')]
hyper_genes = diff_cpg_gene[diff_cpg_gene$Type=='Hyper'&diff_cpg_gene$Group%in%c('TSS200','TSS1500'),c('CpG','Gene')]


#replace the internal kegg database in missMethyl
missmethyl_kegg = readRDS('/Users/Brainfactory/Downloads/missmethyl_kegg.rds')
########### NOTE!!!!############
## To use the subset of KEGG pathways with supergroups from:
## 'Cellular Processes','Environmental Information Processing',
## 'Genetic Information Processing', 'Metabolism','Organismal Systems'
##
## After trace(gometh,edit=T)
## change line 29 to:
## kegg<-readRDS('missmethyl_kegg.rds')
################################
trace(gometh,edit=T)


gst.kegg.hypo.modified <- gometh(sig.cpg=diff_cpg_gene[diff_cpg_gene$Type=='Hypo',]$CpG, 
                                 array.type = "450K",
                                 collection="KEGG",sig.genes =TRUE)

gst.kegg.hyper.modified <- gometh(sig.cpg=diff_cpg_gene[diff_cpg_gene$Type=='Hyper',]$CpG, 
                                  array.type = "450K",
                                  collection="KEGG", genomic.features = c("TSS200","TSS1500"),
                                  prior.prob =FALSE,sig.genes =TRUE)
View(gst.kegg.hypo.modified)
View(gst.kegg.hyper.modified)
write.csv(gst.kegg.hyper.modified,'gst.kegg.hyper.modified.csv')
write.csv(gst.kegg.hypo.modified,'gst.kegg.hypo.modified.csv')
