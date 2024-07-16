#> table(anno_450k$Relation_to_Island)
#Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
#150254   24844   62870  176047   22300   49197 
#> table(anno_450k$UCSC_RefGene_Group)
#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200 
#22737   17494   42685  161677   68984   52283 
#shelf: 47144
#shore: 112067
#cpgisland: 150254
#opensea: 176047

library(ChAMPdata)
data(probe.features)

table(probe.features$feature)
#1stExon   3'UTR   5'UTR    Body     IGR TSS1500  TSS200 
#22737   17494   42685  161677  119717   68984   52283 

table(probe.features$cgi)
#island opensea   shelf   shore 

isle_percentage_norm = data.frame(matrix(0, ncol = 3, nrow = 4))
colnames(isle_percentage_norm) = c('All DM Sites','Hypo Sites','Hyper Sites')
rownames(isle_percentage_norm) = c('Shelf','Shore','CpG Island','Opensea')

isle_percentage_norm[,1]=c(dim(dm_shelf)[1]/47144/(dim(dm_shelf)[1]/47144+dim(dm_shore)[1]/112067+dim(dm_cpgisland)[1]/150254+dim(dm_opensea)[1]/176047),
                           dim(dm_shore)[1]/112067/(dim(dm_shelf)[1]/47144+dim(dm_shore)[1]/112067+dim(dm_cpgisland)[1]/150254+dim(dm_opensea)[1]/176047),
                           dim(dm_cpgisland)[1]/150254/(dim(dm_shelf)[1]/47144+dim(dm_shore)[1]/112067+dim(dm_cpgisland)[1]/150254+dim(dm_opensea)[1]/176047),
                           dim(dm_opensea)[1]/176047/(dim(dm_shelf)[1]/47144+dim(dm_shore)[1]/112067+dim(dm_cpgisland)[1]/150254+dim(dm_opensea)[1]/176047))

isle_percentage_norm[,2]=c(dim(hypo_shelf)[1]/47144/(dim(hypo_shelf)[1]/47144+dim(hypo_shore)[1]/112067+dim(hypo_cpgisland)[1]/150254+dim(hypo_opensea)[1]/176047),
                           dim(hypo_shore)[1]/112067/(dim(hypo_shelf)[1]/47144+dim(hypo_shore)[1]/112067+dim(hypo_cpgisland)[1]/150254+dim(hypo_opensea)[1]/176047),
                           dim(hypo_cpgisland)[1]/150254/(dim(hypo_shelf)[1]/47144+dim(hypo_shore)[1]/112067+dim(hypo_cpgisland)[1]/150254+dim(hypo_opensea)[1]/176047),
                           dim(hypo_opensea)[1]/176047/(dim(hypo_shelf)[1]/47144+dim(hypo_shore)[1]/112067+dim(hypo_cpgisland)[1]/150254+dim(hypo_opensea)[1]/176047))

isle_percentage_norm[,3]=c(dim(hyper_shelf)[1]/47144/(dim(hyper_shelf)[1]/47144+dim(hyper_shore)[1]/112067+dim(hyper_cpgisland)[1]/150254+dim(hyper_opensea)[1]/176047),
                           dim(hyper_shore)[1]/112067/(dim(hyper_shelf)[1]/47144+dim(hyper_shore)[1]/112067+dim(hyper_cpgisland)[1]/150254+dim(hyper_opensea)[1]/176047),
                           dim(hyper_cpgisland)[1]/150254/(dim(hyper_shelf)[1]/47144+dim(hyper_shore)[1]/112067+dim(hyper_cpgisland)[1]/150254+dim(hyper_opensea)[1]/176047),
                           dim(hyper_opensea)[1]/176047/(dim(hyper_shelf)[1]/47144+dim(hyper_shore)[1]/112067+dim(hyper_cpgisland)[1]/150254+dim(hyper_opensea)[1]/176047))

mdata_norm <- expand.grid(colnames(isle_percentage_norm),rownames(isle_percentage_norm))
colnames(mdata_norm) = c('df','region')
mdata_norm$value = rep(0,12)
for (i in 1:12){
  mdata_norm[i,3]=isle_percentage_norm[mdata_norm[i,2],mdata_norm[i,1]]
}

library(wesanderson)
ggplot(mdata_norm, aes(x = df, y = value*100, fill = region)) +
  geom_bar(position="dodge",stat="identity", width = 0.7) +
  theme_minimal(base_size = 12)+
  ggtitle("Isle Region") +
  labs(x = "", y = "Normalized Percentage of Sites", fill = "Isle Region") +
  theme(axis.text.x = element_text(angle = 0))+scale_fill_manual(values = wes_palette("Moonrise2", n = 4))


all_isle_norm = mdata_norm[mdata_norm$df=='All DM Sites',]
p1_norm=ggplot(all_isle_norm, aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black')) +
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('All DM sites')+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = value/4 + c(0, cumsum(value)[-length(value)]),
#              label = percent(value)), size=5)

p2_norm=ggplot(mdata_norm[mdata_norm$df=='Hypo Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','white','white')) +
  scale_fill_brewer(palette="Blues") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

p3_norm=ggplot(mdata_norm[mdata_norm$df=='Hyper Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','white','white')) +
  scale_fill_brewer(palette='Reds') + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))


library(gridExtra)

p_norm <- ggplot()

grid.arrange(grobs = list(p1_norm,p2_norm,p3_norm), 
             top = "CpG Sites Region Distribution",ncol=3, byrow=TRUE)



library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
anno450 = getAnnotation('IlluminaHumanMethylation450kanno.ilmn12.hg19')

#> table(anno_450k$UCSC_RefGene_Group)
#1stExon   3'UTR   5'UTR    Body    TSS1500  TSS200 
#22737      17494   42685   161677   68984   52283 
#gene region
#all sites
all_body = all_df_loc[grepl('Body',all_df_loc$feat.cgi, fixed=TRUE),]
all_igr = all_df_loc[grepl('IGR',all_df_loc$feat.cgi, fixed=TRUE),]
all_200 = all_df_loc[grepl('TSS200',all_df_loc$feat.cgi, fixed=TRUE),]
all_1500 = all_df_loc[grepl('TSS1500',all_df_loc$feat.cgi, fixed=TRUE),]
all_5 = all_df_loc[grepl('5\'UTR',all_df_loc$feat.cgi, fixed=TRUE),]
all_3 = all_df_loc[grepl('3\'UTR',all_df_loc$feat.cgi, fixed=TRUE),]
all_1st = all_df_loc[grepl('1stExon',all_df_loc$feat.cgi, fixed=TRUE),]

#dm sites
dm_body = all_body
dm_igrS = all_igr
dm_T200 = all_200
dm_1500 = all_1500
dm_5UTR = all_5
dm_3UTR = all_3
dm_1ste = all_1st


#hypo sites 
hypo_body = dm_body[dm_body$logFC<0,]
hypo_igrS = dm_igrS[dm_igrS$logFC<0,]
hypo_T200 = dm_T200[dm_T200$logFC<0,]
hypo_1500 = dm_1500[dm_1500$logFC<0,]
hypo_5UTR = dm_5UTR[dm_5UTR$logFC<0,]
hypo_3UTR = dm_3UTR[dm_3UTR$logFC<0,]
hypo_1ste = dm_1ste[dm_1ste$logFC<0,]

#hyper sites
hyper_body = dm_body[dm_body$logFC>0,]
hyper_igrS = dm_igrS[dm_igrS$logFC>0,]
hyper_T200 = dm_T200[dm_T200$logFC>0,]
hyper_1500 = dm_1500[dm_1500$logFC>0,]
hyper_5UTR = dm_5UTR[dm_5UTR$logFC>0,]
hyper_3UTR = dm_3UTR[dm_3UTR$logFC>0,]
hyper_1ste = dm_1ste[dm_1ste$logFC>0,]

gene_percentage_norm = data.frame(matrix(0, ncol = 3, nrow = 7))
colnames(gene_percentage_norm) = c('All DM Sites','Hypo Sites','Hyper Sites')
rownames(gene_percentage_norm) = c('Body','IGR','TSS200','TSS1500',
                              '5\'UTR','3\'UTR','1stExon')


gene_percentage_norm[,1]=c(dim(dm_body)[1]/161677/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_igrS)[1]/119717/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_T200)[1]/52283/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_1500)[1]/68984/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_5UTR)[1]/42685/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_3UTR)[1]/17494/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737),
                           dim(dm_1ste)[1]/22737/(dim(dm_body)[1]/161677+dim(dm_igrS)[1]/119717+dim(dm_T200)[1]/52283+dim(dm_1500)[1]/68984+dim(dm_5UTR)[1]/42685+dim(dm_3UTR)[1]/17494+dim(dm_1ste)[1]/22737))

gene_percentage_norm[,2]=c(dim(hypo_body)[1]/161677/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_igrS)[1]/119717/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_T200)[1]/52283/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_1500)[1]/68984/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_5UTR)[1]/42685/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_3UTR)[1]/17494/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737),
                           dim(hypo_1ste)[1]/22737/(dim(hypo_body)[1]/161677+dim(hypo_igrS)[1]/119717+dim(hypo_T200)[1]/52283+dim(hypo_1500)[1]/68984+dim(hypo_5UTR)[1]/42685+dim(hypo_3UTR)[1]/17494+dim(hypo_1ste)[1]/22737))

gene_percentage_norm[,3]=c(dim(hyper_body)[1]/161677/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_igrS)[1]/119717/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_T200)[1]/52283/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_1500)[1]/68984/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_5UTR)[1]/42685/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_3UTR)[1]/17494/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737),
                           dim(hyper_1ste)[1]/22737/(dim(hyper_body)[1]/161677+dim(hyper_igrS)[1]/119717+dim(hyper_T200)[1]/52283+dim(hyper_1500)[1]/68984+dim(hyper_5UTR)[1]/42685+dim(hyper_3UTR)[1]/17494+dim(hyper_1ste)[1]/22737))

gene_mdata_norm <- expand.grid(colnames(gene_percentage_norm),rownames(gene_percentage_norm))
colnames(gene_mdata_norm) = c('df','region')
gene_mdata_norm$value = rep(0,21)
for (i in 1:21){
  gene_mdata_norm[i,3]=gene_percentage_norm[gene_mdata_norm[i,2],gene_mdata_norm[i,1]]
}

library(wesanderson)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(gene_mdata_norm, aes(x = df, y = value*100, fill = region)) +
  geom_bar(position="dodge",stat="identity", width = 0.7) +
  theme_minimal(base_size = 12)+
  ggtitle("Gene Region") +
  labs(x = "", y = "Percentage of Sites", fill = "Gene Region") +
  theme(axis.text.x = element_text(angle = 0))+
  #scale_fill_brewer(palette = "YlOrRd")
  scale_fill_manual(values =cbPalette)


p4=ggplot(gene_mdata_norm[gene_mdata_norm$df=='All DM Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','white','white','white','white')) +
  scale_fill_brewer(palette="Greens") + theme_void()+
  ggtitle('All DM sites')+
  theme(plot.title = element_text(hjust = 0.5))

p5=ggplot(gene_mdata_norm[gene_mdata_norm$df=='Hypo Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','white','white','white','white')) +
  scale_fill_brewer(palette="Blues") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

p6=ggplot(gene_mdata_norm[gene_mdata_norm$df=='Hyper Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")), size=3.5, position = position_stack(vjust = 0.45),
            color=c('black','black','black','white','white','white','white')) +
  scale_fill_brewer(palette='Reds') + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))


library(gridExtra)

p <- ggplot()

grid.arrange(grobs = list(p4,p5,p6), top = "Normalized CpG Sites Distribution across Gene Regions",ncol=3, byrow=TRUE)

