
ggarrange(ggarrange(isle_p1,isle_p2,isle_p3,isle_p4,ncol=4,common.legend = TRUE, legend="right"),
          ggarrange(region_p1,region_p2,region_p3,region_p4,ncol=4,common.legend = TRUE, legend="right"),
          nrow = 2)




sovplot <- function(restab = MSSmean, clustername = 'Maternal Obesity', plottype = 'MSS', 
                    textsize = 20){
  
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername)#, '_Control')
  
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        ggtitle('Source of Variance (Type 3 Anova)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
    
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        #ggtitle('Source of Variance (Before)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') + 
        #ggtitle('Source of Variance (After)') + 
        ylab(ytitle) + 
        xlab('') + 
        scale_fill_discrete(guide = FALSE) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize))
    )
  }
  
}

load("cobre_final_res.Rdata",verbose = TRUE)
Ftab_before = Ftab
Fmean_before = Fmean
finalvars_before = finalvars
sov_before = sovplot(restab = Fmean_before, plottype = 'F', textsize = 10)+
  ggtitle('Source of Variance (Before)')

load("cobre_final_res_adj_allparityhemo_sov2.rdata",verbose = TRUE)
sov_after = sovplot(restab = Fmean, plottype = 'F', textsize = 10)+
  ggtitle('Source of Variance (After)')

volcanoplot = readRDS("volcanoplot.RDS")



pieplots = ggarrange(ggarrange(isle_p1,isle_p2,isle_p3,isle_p4,ncol=4,common.legend = TRUE, legend="right"),
                     ggarrange(region_p1,region_p2,region_p3,region_p4,ncol=4,common.legend = TRUE, legend="right"),
                     nrow = 2)



figure2 = ggarrange(ggarrange(sov_before,sov_after,ncol=2,labels = c("A","B")),
ggarrange(volcanoplot,pieplots,ncol=2,labels = c("C","D")),nrow = 2)



ggsave(file="COBRE_FINAL_FINAL_FIGURES\\New Fig 2.png", 
       plot=figure2, width=12, height=10, dpi = 300, units = "in", device='png')


ggsave(file="New Fig 2 pieplot.svg", 
       plot=pieplots, width=6, height=5)



