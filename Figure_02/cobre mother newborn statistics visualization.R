

p1 = ggplot(cobre_mom_baby, aes(x = `Weight.Kg.`, y = Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle('P=0.0485 (*)')+theme_minimal()+
  xlab('Baby Weight (kg)')+ylab('Sample Group')+
  scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9))

p2 = ggplot(cobre_mom_baby, aes(x = `HeadCirc.cm.`, y = Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+
  ggtitle("P=0.1139 (n.s.)") +
  theme_minimal()+xlab('Head Circumference (cm)')+ylab('Sample Group')+
  scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())


p3 = ggplot(cobre_mom_baby, aes(x = `Length.cm.`, y = Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle("P=0.5753 (n.s.)") +
  theme_minimal()+xlab('Baby Length (cm)')+ylab('Sample Group')+
  scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())


p4 = ggplot(cobre_mom_baby, aes(x = `APGAR.5`, y = Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle("P=0.1856 (n.s.)") +
  theme_minimal()+xlab('APGAR Score 5 min')+ylab('Sample Group')+
  scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())


p5 = ggplot(cobre_pd[,c('Sample_Group','Mat_Age')], aes(x=Mat_Age, y=Sample_Group, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle(paste0("P = ",round(t.test(cobre_pd[cobre_pd$Sample_Group=='control','Mat_Age'],
                                     cobre_pd[cobre_pd$Sample_Group=='obese','Mat_Age'])$p.value,4), " (n.s.)")) +
  xlab("Maternal Age (Year)")+
  ylab("Sample Group")+scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9))


p6 = ggplot(cobre_pd[,c('Sample_Group','Gestational_Age')], aes(x=Gestational_Age, y=Sample_Group, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle(paste0("P = ",round(t.test(cobre_pd[cobre_pd$Sample_Group=='control','Gestational_Age'],
                                     cobre_pd[cobre_pd$Sample_Group=='obese','Gestational_Age'])$p.value,4), " (n.s.)")) +
  xlab("Gestational Age (Week)")+
  ylab("Sample Group")+scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())

p7 = ggplot(cobre_pd[,c('Sample_Group','Net_Weight_Gain')], aes(x=Net_Weight_Gain, y=Sample_Group, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle(paste0("P = ",round(t.test(cobre_pd[cobre_pd$Sample_Group=='control','Net_Weight_Gain'],
                                     cobre_pd[cobre_pd$Sample_Group=='obese','Net_Weight_Gain'])$p.value,4), " (n.s.)")) +
  xlab("Net Weight Gain (lb)")+
  ylab("Sample Group")+scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())


p8 = ggplot(cobre_pd[,c('Sample_Group','Hemoglobin')], aes(x=Hemoglobin, y=Sample_Group, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_density_ridges(aes(fill = Sample_Group), alpha=0.2, scale=0.8,color='transparent') +
  geom_boxplot(aes(fill = Sample_Group), width = 0.06)+theme_minimal()+
  ggtitle(paste0("P = ",round(t.test(cobre_pd[cobre_pd$Sample_Group=='control','Hemoglobin'],
                                     cobre_pd[cobre_pd$Sample_Group=='obese','Hemoglobin'])$p.value,4), " (n.s.)"))+
  xlab("Hemoglobin (g/dl)")+
  ylab("Sample Group")+scale_y_discrete(expand = c(0,0.2))+ 
  theme(plot.title = element_text(size=9),axis.title.y = element_blank())



p9 = ggplot(cobre_pd[,c('Sample_Group','Mat_Ethnicity')],aes(x=Mat_Ethnicity, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_bar(stat="count", position=position_dodge(), alpha=0.7,width = 0.5)+
  theme_minimal()+
  ggtitle(paste0("P = ",round(chisq.test(table(cobre_pd[,c('Sample_Group','Mat_Ethnicity')]), correct = F)$p.value,4), " (***)"))+
  ylab("Counts")+ylim(c(0,25))+
  xlab("Maternal Ethnicity")+ #coord_flip()+ 
  theme(plot.title = element_text(size=9))#+scale_x_discrete(labels= c("NHPI","Caucasian","Asian"))

p10 = ggplot(cobre_pd[,c('Sample_Group','Pat_Ethnicity')],aes(x=Pat_Ethnicity, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_bar(stat="count", position=position_dodge(), alpha=0.7,width = 0.6)+
  theme_minimal()+
  ggtitle(paste0("P = ",round(chisq.test(table(cobre_pd[,c('Sample_Group','Pat_Ethnicity')]), correct = F)$p.value,4)," (***)"))+
  ylab("Counts")+ylim(c(0,25))+
  xlab("Paternal Ethnicity")+ #coord_flip()+ 
  theme(plot.title = element_text(size=9))#+scale_x_discrete(labels= c("NHPI","Caucasian","Asian"))

p11= ggplot(cobre_pd[,c('Sample_Group','Sex')],aes(x=Sex, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_bar(stat="count", position=position_dodge(), alpha=0.7,width = 0.6)+
  theme_minimal()+
  ggtitle(paste0("P = ",round(chisq.test(table(cobre_pd[,c('Sample_Group','Sex')]), correct = F)$p.value,4)," (n.s.)"))+
  ylab("Counts")+ylim(c(0,25))+
  xlab("Baby Sex")+ #coord_flip()+
  scale_x_discrete(labels= c("Male","Female"))+ 
  theme(plot.title = element_text(size=9))

p12 = ggplot(cobre_pd[,c('Sample_Group','Parity')],aes(x=Parity, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_bar(stat="count", position=position_dodge(), alpha=0.7,width = 0.7)+
  theme_minimal()+
  ggtitle(paste0("P = ",round(chisq.test(table(cobre_pd[,c('Sample_Group','Parity')]), correct = F)$p.value,4)," (***)"))+
  ylab("Counts")+ylim(c(0,25))+
  xlab("Parity")+ #coord_flip()+ 
  theme(plot.title = element_text(size=9))#+scale_x_discrete(labels= c("Male","Female"))


p13 = ggplot(as.data.frame(table(cobre_pd[,c('Sample_Group','Gravidity')])),aes(x=Gravidity,y=Freq, fill=Sample_Group)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))+
  geom_bar(stat="identity", alpha=0.7,position=position_dodge())+
  theme_minimal()+
  ggtitle(paste0("P = ",round(chisq.test(table(cobre_pd[,c('Sample_Group','Gravidity')]), correct = F)$p.value,4)," (**)"))+
  ylab("Counts")+ylim(c(0,25))+
  xlab("Gravidity")+ #coord_flip()+ 
  theme(plot.title = element_text(size=9))#+scale_x_discrete(labels= c("Male","Female"))


png('COBRE_FINAL_FINAL_FIGURES\\New Fig 1 baby and mother stats.png', height=3000, width=4000,res=300)
ggarrange(ggarrange(p1, p2, p3,p4, ncol = 4,labels = c("A","B", "C","D"), common.legend = TRUE, legend="right" ), 
          #ggplot() + theme_void()+theme(plot.margin = margin(0,0,0,0, 'cm')),
          ggarrange(p5,p6,p7,p8, ncol = 4,labels = c("E","F", "G","H"), common.legend = TRUE, legend="right" ),   # First row with scatter plot
          #ggplot() + theme_void()+theme(plot.margin = margin(0,0,0,0, 'cm')),
          ggarrange(p11,p9,p10,p12,p13,  ncol = 5, labels = c("I","J", "K","L",'M'), common.legend = TRUE, legend="right" ), # Second row with box and dot plots
          nrow = 3  , 
          heights = c(1,1,0.8),
          common.legend = TRUE   ,legend="right"                               # Labels of the scatter plot
) 
dev.off()


image = ggarrange(ggarrange(p1, p2, p3,p4, ncol = 4,labels = c("A","B", "C","D"), common.legend = TRUE, legend="right" ), 
                  #ggplot() + theme_void()+theme(plot.margin = margin(0,0,0,0, 'cm')),
                  ggarrange(p5,p6,p7,p8, ncol = 4,labels = c("E","F", "G","H"), common.legend = TRUE, legend="right" ),   # First row with scatter plot
                  #ggplot() + theme_void()+theme(plot.margin = margin(0,0,0,0, 'cm')),
                  ggarrange(p11,p9,p10,p12,p13,  ncol = 5, labels = c("I","J", "K","L",'M'), common.legend = TRUE, legend="right" ), # Second row with box and dot plots
                  nrow = 3  , 
                  heights = c(1,1,0.8),
                  common.legend = TRUE   ,legend="right"                               # Labels of the scatter plot
) 
ggsave(file="COBRE_FINAL_FINAL_FIGURES\\New3 Fig 1 svg baby and mother stats.svg", 
       plot=image, width=10, height=8)





png('COBRE_FINAL_FINAL_FIGURES\\Final Fig 1 Rearranged New  baby and mother stats.png', height=3000, width=3500,res=300)
ggarrange(
           
          ggarrange(p11,p9,p10,p12,p13,  ncol = 5, labels = c("A","B", "C","D",'E'), common.legend = TRUE, legend="right" ), 
          ggarrange(p5,p6,p7,p8, ncol = 4,labels = c("F","G", "H","I"), common.legend = TRUE, legend="right" ), 
          ggarrange(p1, p2, p3,p4, ncol = 4,labels = c("J","K", "L","M"), common.legend = TRUE, legend="right" ), 
          nrow = 3  , 
          heights = c(0.8,1,1),
          common.legend = TRUE   ,legend="right"                               # Labels of the scatter plot
) 
dev.off()


