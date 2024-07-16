isle_df = read.csv("C:\\Users\\yhdu\\OneDrive - Michigan Medicine\\Desktop\\cobre cpg probes location pie chart\\gene_isle_unnormalized_data_piechart.csv")
region_df = read.csv("C:\\Users\\yhdu\\OneDrive - Michigan Medicine\\Desktop\\cobre cpg probes location pie chart\\gene_region_unnormalized_data_piechart.csv")


isle_p1 = ggplot(isle_df[isle_df$df=='Illumina450K',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black')) +
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('All Illumina450K')+
  theme(plot.title = element_text(hjust = 0.5))


isle_p2 = ggplot(isle_df[isle_df$df=='All DM Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black')) +
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('All DM Sites')+
  theme(plot.title = element_text(hjust = 0.5))

isle_p3 = ggplot(isle_df[isle_df$df=='Hypermethylation',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black')) +
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))

isle_p4 = ggplot(isle_df[isle_df$df=='Hypomethylation',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%")),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black')) +
  scale_fill_brewer(palette="Pastel1") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

library(ggpubr)
ggarrange(isle_p1,isle_p2,isle_p3,isle_p4,ncol=4,common.legend = TRUE, legend="right")


###########
###########
# regions
###########
###########

region_p1 = ggplot(region_df[region_df$df=='Illumina450K',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%"), x = 1.2),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black','black','black','black')) +
  scale_fill_brewer(palette="Pastel2") + theme_void()+
  ggtitle('All Illumina450K')+
  theme(plot.title = element_text(hjust = 0.5))


region_p2 = ggplot(region_df[region_df$df=='All DM Sites',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%"), x = 1.2),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black','black','black','black')) +
  scale_fill_brewer(palette="Pastel2") + theme_void()+
  ggtitle('All DM Sites')+
  theme(plot.title = element_text(hjust = 0.5))

region_p3 = ggplot(region_df[region_df$df=='Hypermethylation',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%"), x = 1.2),size=3.5,position = position_stack(vjust = 0.45),
            color=c('black','black','black','black','black','black','black')) +
  scale_fill_brewer(palette="Pastel2") + theme_void()+
  ggtitle('Hyper Sites')+
  theme(plot.title = element_text(hjust = 0.5))

region_p4 = ggplot(region_df[region_df$df=='Hypomethylation',], aes(x = '', y = value, fill = region))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(value*100,1),"%"), x = 1.2),size=3.5,position = position_stack(vjust = 0.5),
            color=c('black','black','black','black','black','black','black')) +
  scale_fill_brewer(palette="Pastel2") + theme_void()+
  ggtitle('Hypo Sites')+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(ggarrange(isle_p1,isle_p2,isle_p3,isle_p4,ncol=4,common.legend = TRUE, legend="right"),
          ggarrange(region_p1,region_p2,region_p3,region_p4,ncol=4,common.legend = TRUE, legend="right"),
          nrow = 2)


ggarrange(ggarrange(isle_p1,isle_p2,isle_p3,isle_p4,ncol=4,common.legend = TRUE, legend="right"),
          ggarrange(region_p1,region_p2,region_p3,region_p4,ncol=4,common.legend = TRUE, legend="right"),
          nrow = 2)
isle_count = read.csv("/Users/Brainfactory/gene_isle_count_data_piechart.csv")
region_count = read.csv("/Users/Brainfactory/gene_region_count_data_piechart.csv")

library(reshape2)

isle_count_df = melt(isle_count,value.name = 'Isle')
colnames(isle_count_df) = c('Isle','df','Count')

observed1 = cbind(isle_count_df[isle_count_df$df=='Illumina450K','Count'],isle_count_df[isle_count_df$df=='DM','Count'])
chisq.test(observed1)
observed2 = cbind(isle_count_df[isle_count_df$df=='Illumina450K','Count'],isle_count_df[isle_count_df$df=='Hypermethylation','Count'])
chisq.test(observed2)
observed3 = cbind(isle_count_df[isle_count_df$df=='Illumina450K','Count'],isle_count_df[isle_count_df$df=='Hypomethylation','Count'])
chisq.test(observed3)
observed4 = cbind(isle_count_df[isle_count_df$df=='Hypermethylation','Count'],isle_count_df[isle_count_df$df=='Hypomethylation','Count'])
chisq.test(observed4)


