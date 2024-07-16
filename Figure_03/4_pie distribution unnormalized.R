sigg.limma = readRDS('cobre_limma_sig_1w5.RDS')
all_df = sigg.limma
all_df$CpG = rownames(all_df)
load('CpG_Gene_loc_convert.RData')
colnames(df_cpg_gene_loc) <- c('CpG','CHR','gene','feat.cgi')
all_df_loc = merge(all_df,df_cpg_gene_loc,by ='CpG')


#isle region
#all sites
all_shelf = all_df_loc[grepl('shelf',all_df_loc$feat.cgi, fixed=TRUE),]
all_shore = all_df_loc[grepl('shore',all_df_loc$feat.cgi, fixed=TRUE),]
all_cpgisland = all_df_loc[grepl('island',all_df_loc$feat.cgi, fixed=TRUE),]
all_opensea = all_df_loc[grepl('opensea',all_df_loc$feat.cgi, fixed=TRUE),]
#dm sites
dm_shelf = all_shelf[grepl('shelf',all_shelf$feat.cgi, fixed=TRUE),]
dm_shore = all_shore[grepl('shore',all_shore$feat.cgi, fixed=TRUE),]
dm_cpgisland = all_cpgisland[grepl('island',all_cpgisland$feat.cgi, fixed=TRUE),]
dm_opensea = all_opensea[grepl('opensea',all_opensea$feat.cgi, fixed=TRUE),]



#hypo sites 
hypo_shelf = dm_shelf[dm_shelf$logFC<0,]
hypo_shore = dm_shore[dm_shore$logFC<0,]
hypo_cpgisland = dm_cpgisland[dm_cpgisland$logFC<0,]
hypo_opensea = dm_opensea[dm_opensea$logFC<0,]
#hyper sites
hyper_shelf = dm_shelf[dm_shelf$logFC>0,]
hyper_shore = dm_shore[dm_shore$logFC>0,]
hyper_cpgisland = dm_cpgisland[dm_cpgisland$logFC>0,]
hyper_opensea = dm_opensea[dm_opensea$logFC>0,]

isle_percentage = data.frame(matrix(0, ncol = 3, nrow = 4))
colnames(isle_percentage) = c('All DM Sites','Hypo Sites','Hyper Sites')
rownames(isle_percentage) = c('Shelf','Shore','CpG Island','Opensea')

#isle_percentage[,1]=c(dim(all_shelf)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_shore)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_cpgisland)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_opensea)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]))

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

isle_percentage[,1]=c(dim(dm_shelf)[1]/(dim(dm_shelf)[1]+dim(dm_shore)[1]+dim(dm_cpgisland)[1]+dim(dm_opensea)[1]),
                      dim(dm_shore)[1]/(dim(dm_shelf)[1]+dim(dm_shore)[1]+dim(dm_cpgisland)[1]+dim(dm_opensea)[1]),
                      dim(dm_cpgisland)[1]/(dim(dm_shelf)[1]+dim(dm_shore)[1]+dim(dm_cpgisland)[1]+dim(dm_opensea)[1]),
                      dim(dm_opensea)[1]/(dim(dm_shelf)[1]+dim(dm_shore)[1]+dim(dm_cpgisland)[1]+dim(dm_opensea)[1]))

isle_percentage[,2]=c(dim(hypo_shelf)[1]/(dim(hypo_shelf)[1]+dim(hypo_shore)[1]+dim(hypo_cpgisland)[1]+dim(hypo_opensea)[1]),
                      dim(hypo_shore)[1]/(dim(hypo_shelf)[1]+dim(hypo_shore)[1]+dim(hypo_cpgisland)[1]+dim(hypo_opensea)[1]),
                      dim(hypo_cpgisland)[1]/(dim(hypo_shelf)[1]+dim(hypo_shore)[1]+dim(hypo_cpgisland)[1]+dim(hypo_opensea)[1]),
                      dim(hypo_opensea)[1]/(dim(hypo_shelf)[1]+dim(hypo_shore)[1]+dim(hypo_cpgisland)[1]+dim(hypo_opensea)[1]))

isle_percentage[,3]=c(dim(hyper_shelf)[1]/(dim(hyper_shelf)[1]+dim(hyper_shore)[1]+dim(hyper_cpgisland)[1]+dim(hyper_opensea)[1]),
                      dim(hyper_shore)[1]/(dim(hyper_shelf)[1]+dim(hyper_shore)[1]+dim(hyper_cpgisland)[1]+dim(hyper_opensea)[1]),
                      dim(hyper_cpgisland)[1]/(dim(hyper_shelf)[1]+dim(hyper_shore)[1]+dim(hyper_cpgisland)[1]+dim(hyper_opensea)[1]),
                      dim(hyper_opensea)[1]/(dim(hyper_shelf)[1]+dim(hyper_shore)[1]+dim(hyper_cpgisland)[1]+dim(hyper_opensea)[1]))


mdata <- expand.grid(colnames(isle_percentage),rownames(isle_percentage))
colnames(mdata) = c('df','region')
mdata$value = rep(0,12)
for (i in 1:12){
  mdata[i,3]=isle_percentage[mdata[i,2],mdata[i,1]]
}

library(wesanderson)
ggplot(mdata, aes(x = df, y = value*100, fill = region)) +
  geom_bar(position="dodge",stat="identity", width = 0.7) +
  theme_minimal(base_size = 12)+
  ggtitle("Isle Region") +
  labs(x = "", y = "Percentage of Sites", fill = "Isle Region") +
  theme(axis.text.x = element_text(angle = 0))+
  scale_fill_manual(values = wes_palette("Moonrise2", n = 4))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

all_isle = mdata[mdata$df=='All DM Sites',]
p1=ggplot(all_isle, aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('All DM sites')+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(aes(y = value/4 + c(0, cumsum(value)[-length(value)]),
#              label = percent(value)), size=5)

p2=ggplot(mdata[mdata$df=='Hypo Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

p3=ggplot(mdata[mdata$df=='Hyper Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette='Reds') + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))


library(gridExtra)

p <- ggplot()

grid.arrange(grobs = list(p1,p2,p3), top = "CpG Sites Region Distribution",ncol=3, byrow=TRUE)




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

gene_percentage = data.frame(matrix(0, ncol = 3, nrow = 7))
colnames(gene_percentage) = c('All DM Sites','Hypo Sites','Hyper Sites')
rownames(gene_percentage) = c('Body','IGR','TSS200','TSS1500',
                              '5\'UTR','3\'UTR','1stExon')

#isle_percentage[,1]=c(dim(all_shelf)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_shore)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_cpgisland)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]),
#                      dim(all_opensea)[1]/(dim(all_shelf)[1]+dim(all_shore)[1]+dim(all_cpgisland)[1]+dim(all_opensea)[1]))


gene_percentage[,1]=c(dim(dm_body)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_igrS)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_T200)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_1500)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_5UTR)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_3UTR)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]),
                      dim(dm_1ste)[1]/(dim(dm_body)[1]+dim(dm_igrS)[1]+dim(dm_T200)[1]+dim(dm_1500)[1]+dim(dm_5UTR)[1]+dim(dm_3UTR)[1]+dim(dm_1ste)[1]))

gene_percentage[,2]=c(dim(hypo_body)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_igrS)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_T200)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_1500)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_5UTR)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_3UTR)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]),
                      dim(hypo_1ste)[1]/(dim(hypo_body)[1]+dim(hypo_igrS)[1]+dim(hypo_T200)[1]+dim(hypo_1500)[1]+dim(hypo_5UTR)[1]+dim(hypo_3UTR)[1]+dim(hypo_1ste)[1]))

gene_percentage[,3]=c(dim(hyper_body)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_igrS)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_T200)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_1500)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_5UTR)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_3UTR)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]),
                      dim(hyper_1ste)[1]/(dim(hyper_body)[1]+dim(hyper_igrS)[1]+dim(hyper_T200)[1]+dim(hyper_1500)[1]+dim(hyper_5UTR)[1]+dim(hyper_3UTR)[1]+dim(hyper_1ste)[1]))


gene_mdata <- expand.grid(colnames(gene_percentage),rownames(gene_percentage))
colnames(gene_mdata) = c('df','region')
gene_mdata$value = rep(0,21)
for (i in 1:21){
  gene_mdata[i,3]=gene_percentage[gene_mdata[i,2],gene_mdata[i,1]]
}

library(wesanderson)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(gene_mdata, aes(x = df, y = value*100, fill = region)) +
  geom_bar(position="dodge",stat="identity", width = 0.7) +
  theme_minimal(base_size = 12)+
  ggtitle("Gene Region") +
  labs(x = "", y = "Percentage of Sites", fill = "Gene Region") +
  theme(axis.text.x = element_text(angle = 0))+
  #scale_fill_brewer(palette = "YlOrRd")
  scale_fill_manual(values =cbPalette)


p4=  
  ggplot(gene_mdata[gene_mdata$df=='All DM Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Greens") + theme_void()+
  ggtitle('All DM sites')+
  theme(plot.title = element_text(hjust = 0.5))

p5=ggplot(gene_mdata[gene_mdata$df=='Hypo Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Blues") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

p6=ggplot(gene_mdata[gene_mdata$df=='Hyper Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette='Reds') + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))


library(gridExtra)

p <- ggplot()

grid.arrange(grobs = list(p4,p5,p6), top = "CpG Sites Distribution across Gene Regions",ncol=3, byrow=TRUE)
