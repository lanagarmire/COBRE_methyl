#library(pandaR)
#data(pandaToyData)
#pandaResult <- panda(pandaToyData$motif, pandaToyData$expression, pandaToyData$ppi)

library(STRINGdb)
library(PANDA)
library('EnsDb.Hsapiens.v79')
data('dfPPI')
data("GENE2GOtopLite")
data("GENE2KEGG")
data("KEGGID2NAME")

#cobre_limma_sig = readRDS('/Users/Brainfactory/cobre_limma_sig.RDS')
diff_cpg_gene = readRDS('/Users/Brainfactory/COBRE_diff_cpg_gene1w5.rds')

diff_gene_tss200_tss1500 = diff_cpg_gene[diff_cpg_gene$Group%in%c('TSS200','TSS1500'),]

#cobre_limma_sig_tss200 = cobre_limma_sig[rownames(cobre_limma_sig)%in%diff_gene_tss200$CpG,c('logFC','adj.P.Val')]
#cobre_limma_sig_tss200_tss1500 = cobre_limma_sig[rownames(cobre_limma_sig)%in%diff_gene_tss200_tss1500$CpG,c('logFC','adj.P.Val')]
#cobre_limma_sig_tss200_tss1500 = cobre_limma_sig[rownames(cobre_limma_sig)%in%diff_gene_tss200_tss1500$CpG,c('logFC','adj.P.Val')]


diff_gene = diff_cpg_gene[complete.cases(diff_cpg_gene),]
#diff_gene_tss200 = diff_cpg_gene[diff_cpg_gene$Group=='TSS200',]
diff_gene_tss200_tss1500 = diff_cpg_gene[diff_cpg_gene$Group%in%c('TSS200','TSS1500'),]
hyper_gene = diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$logFC>0,]$Gene
hypo_gene = diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$logFC<0,]$Gene

geneSymbols <-  unique(diff_gene_tss200_tss1500$Gene)  #new ethnic 1489; unique 1363
geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
geneIDs2 <- geneIDs2[complete.cases(geneIDs2),] #8366 #new ethnic 6104 #TSS200 1413: one symbol map to multiple id
geneIDs2$Dir <- 'Hyper'
geneIDs2$Dir[geneIDs2$SYMBOL %in% hypo_gene] <- 'Hypo'


gene_id_strindb = diff_gene_tss200_tss1500[,c('Gene','Type')]
colnames(gene_id_strindb) = c('SYMBOL','Type')
string_db <- STRINGdb$new(version="10", species=9606, input_directory="", score_threshold = 500)
#stringids <- string_db$map(geneIDs2, "SYMBOL", removeUnmappedRows = TRUE)
stringids <- string_db$map(gene_id_strindb, "SYMBOL", removeUnmappedRows = TRUE)
ids <- unique(stringids$STRING_id) #TSS200: 1306 protein ids #TSSBOTH: 2579 #newboth1w5 3827
ppis <- string_db$get_interactions(ids)
tmp <- ppis
ppis <- merge(ppis, stringids, by.x = 'from', by.y = 'STRING_id')
ppinames <- c('from', 'to', 'combined_score', 'SYMBOL', 'Type')
ppis <- ppis[ppinames]
names(ppis) <- c('from', 'to', 'combined_score', 'from_symbol', 'from_dir')
ppis <- merge(ppis, stringids, by.x = 'to', by.y = 'STRING_id')
names(ppis)[c(6, 7)] <- c('to_symbol', 'to_dir')

pandappi <- ppis[c('from_symbol', 'to_symbol')]

ppi_tss200<- ppis[c('from_symbol', 'to_symbol')]
ppi_tss200_tss1500 <- ppis[c('from_symbol', 'to_symbol')]


getppi_string <- function(gene_data){
  gene_id_strindb = gene_data[,c('Gene','Type')]
  colnames(gene_id_strindb) = c('SYMBOL','Type')
  string_db <- STRINGdb$new(version="10", species=9606, input_directory="",score_threshold=500)
  stringids <- string_db$map(gene_id_strindb, "SYMBOL", removeUnmappedRows = TRUE)
  ids <- unique(stringids$STRING_id)
  ppis <- string_db$get_interactions(ids)
  tmp <- ppis
  ppis <- merge(ppis, stringids, by.x = 'from', by.y = 'STRING_id')
  ppinames <- c('from', 'to', 'combined_score', 'SYMBOL', 'Type')
  ppis <- ppis[ppinames]
  names(ppis) <- c('from', 'to', 'combined_score', 'from_symbol', 'from_dir')
  ppis <- merge(ppis, stringids, by.x = 'to', by.y = 'STRING_id') 
  names(ppis)[c(6, 7)] <- c('to_symbol', 'to_dir')
  return(unique(ppis[c('from_symbol', 'to_symbol')]))
}


#OrderAll_tss200 <- SignificantPairs(PPIdb = pandappi) #c('from_symbol', 'to_symbol')
OrderAll_tss200_tss1500 <- SignificantPairs(PPIdb = ppi_tss200_tss1500)
#GP = GOpredict(Pfile=OrderAll_tss200_tss1500, PPIdb=dfPPI, Gene2Annotation=GENE2GOtopLite, p_value=0.001)
KP_tss200_tss1500 = KEGGpredict(Pfile=OrderAll_tss200_tss1500, PPIdb=dfPPI, Gene2Annotation=GENE2KEGG, p_value=0.005, IDtoNAME=KEGGID2NAME)

sum(OrderAll_tss200_tss1500$Sym_A%in%names(GENE2KEGG))
sum(OrderAll_tss200_tss1500$Sym_B%in%names(GENE2KEGG))

GP_type = diff_gene_tss200[diff_gene_tss200$Gene%in%GP$Symbol,c('Gene','Type')]
GP_tss200 = merge(GP,GP_type,by.x='Symbol',by.y='Gene')
write.csv(GP_tss200,'GP_tss200_GO_PANDA.csv',row.names = F)

OrderAll_tss200_sub = subset(OrderAll_tss200,Sym_A%in%GP_tss200$Symbol&Sym_B%in%GP_tss200$Symbol)
OrderAll_tss200_tss1500_sub = subset(OrderAll_tss200_tss1500,Sym_A%in%GP_tss200$Symbol&Sym_B%in%GP_tss200$Symbol)

OrderAll_toy=SignificantPairs(PPIdb=dfPPI)


gene_kegg = diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$Gene%in%KP_tss200_tss1500$Symbol,]
gene_kegg_remove_dup  = c()
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
for (i in unique(gene_kegg$Gene)) {
  res_type=''
  temp = gene_kegg[gene_kegg$Gene==i,]
  if(dim(temp)[1]>2){
    res_type = Mode(temp$Type)
  }
  if(dim(temp)[1]==2){
    res_type = temp[temp$logFC==max(temp$logFC),]$Type
  }
  if(dim(temp)[1]==1){
    res_type = temp$Type
  }
  gene_kegg_remove_dup = c(gene_kegg_remove_dup,res_type)
}
gene_kegg_remove_dup = cbind(unique(gene_kegg$Gene),gene_kegg_remove_dup)
colnames(gene_kegg_remove_dup) = c('Gene','Type')

#KP_type = unique(diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$Gene%in%KP_tss200_tss1500$Symbol,c('Gene','Type')])
#KP_type = KP_type[!KP_type$Gene=='DIMT1L',]
#KP_tss200_tss1500 = KP_tss200_tss1500[!KP_tss200_tss1500$Symbol=='DIMT1L',]
#KP_tss200_tss1500 = merge(KP_tss200_tss1500,KP_type,by.x='Symbol',by.y='Gene')

KP_tss200_tss1500 = merge(KP_tss200_tss1500,gene_kegg_remove_dup,by.x='Symbol',by.y='Gene')
OrderAll_tss200_tss1500_sub = subset(OrderAll_tss200_tss1500,Sym_A%in%KP_tss200_tss1500$Symbol&Sym_B%in%KP_tss200_tss1500$Symbol)

OrderAll_tss200_tss1500_sub$source_anno = ''
OrderAll_tss200_tss1500_sub$end_anno = ''

for (i in 1:dim(OrderAll_tss200_tss1500_sub)[1]) {
  temp = OrderAll_tss200_tss1500_sub[i,]
  source = temp$Sym_A
  end = temp$Sym_B
  temp$source_anno = KP_tss200_tss1500[KP_tss200_tss1500$Symbol==source,]$Type
  temp$end_anno = KP_tss200_tss1500[KP_tss200_tss1500$Symbol==end,]$Type
  OrderAll_tss200_tss1500_sub[i,] = temp
}


write.csv(KP_tss200_tss1500,'KEGG_tss200_tss1500_PANDA_anno_050322.csv',quote = F,row.names = F)
write.csv(OrderAll_tss200_tss1500_sub,'KEGG_tss200_tss1500_PANDA_nodes_050322.csv',quote = F,row.names = F)

write.csv(KP_tss200_tss1500,'KEGG_tss200_tss1500_0.005_PANDA_anno_051722.csv',quote = F,row.names = F)
write.csv(OrderAll_tss200_tss1500_sub,'KEGG_tss200_tss1500_0.005_PANDA_nodes_051722.csv',quote = F,row.names = F)


KP_tss200_tss1500 = read.csv('/Users/Brainfactory/KEGG_tss200_tss1500_PANDA_anno_050322.csv')
OrderAll_tss200_tss1500_sub = read.csv('/Users/Brainfactory/KEGG_tss200_tss1500_PANDA_nodes_050322.csv')

only_one_path = table(KP_tss200_tss1500$PathName)[table(KP_tss200_tss1500$PathName)<5]
gene_only_one_path = KP_tss200_tss1500[KP_tss200_tss1500$PathName%in%names(only_one_path),]
gene_multi_path = KP_tss200_tss1500[!KP_tss200_tss1500$PathName%in%names(only_one_path),]

OrderAll_tss200_tss1500_sub_multi_path = subset(OrderAll_tss200_tss1500_sub,(!Sym_A%in%gene_only_one_path$Symbol)&(!Sym_B%in%gene_only_one_path$Symbol))

write.csv(gene_multi_path,'/Users/Brainfactory/KEGG_tss200_tss1500_PANDA_anno_multi_path_more2_050322.csv',quote = F,row.names = F)
write.csv(OrderAll_tss200_tss1500_sub_multi_path,'/Users/Brainfactory/KEGG_tss200_tss1500_PANDA_nodes_multi_path_more2_050322.csv',quote = F,row.names = F)

write.csv(gene_multi_path,'KEGG_tss200_tss1500_0.005more5_PANDA_anno_051722.csv',quote = F,row.names = F)
write.csv(OrderAll_tss200_tss1500_sub_multi_path,'KEGG_tss200_tss1500_0.005more5_PANDA_nodes_051722.csv',quote = F,row.names = F)

table(gene_multi_path$PathName)

OrderAll_tss200_tss1500_sub_multi_path

KP_type = unique(diff_gene_tss200_tss1500[diff_gene_tss200_tss1500$Gene%in%KP$Symbol,c('Gene','Type')])

#TSS200: 
#double = c('DIMT1L')
checkdouble = table(KP_type$Gene)
KP_type = KP_type[!KP_type$Gene=='DIMT1L',]

#TSS1500:
doubles = names(checkdouble[checkdouble>1])
KP_type = KP_type[!KP_type$Gene%in%doubles,]
KP = KP[!KP$Symbol%in%doubles,]
KP_tss200_1500 = merge(KP,KP_type,by.x='Symbol',by.y='Gene')
OrderAll_tss200_tss1500_sub = subset(OrderAll_tss200_tss1500,Sym_A%in%KP_tss200_1500$Symbol&Sym_B%in%KP_tss200_1500$Symbol)
write.csv(KP_tss200_1500,'KEGG_tss200_tss1500_PANDA_anno.csv',quote = F,row.names = F)
write.csv(OrderAll_tss200_tss1500_sub,'KEGG_tss200_tss1500_PANDA_nodes.csv',quote = F,row.names = F)



library(EnsDb.Hsapiens.v86)


diff_gene_tss200_tss1500 = diff_cpg_gene[diff_cpg_gene$Group%in%c('TSS200','TSS1500'),]




