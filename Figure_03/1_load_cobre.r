options(connectionObserver = NULL)
library("ChAMP")
library(parallel)
library(ggplot2)
library(car)
library("lumi")
library("dplyr")
library("AnalyzeFMRI")


dataDir = 'cobre/idat'
myLoad <- champ.load(dataDir,arraytype='450K')
saveRDS(myLoad,'myLoad_final_3.6ChAMP_2.14.0.rds')
champ.QC()
myLoad$pd=readRDS('cobre_pd.rds')
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5) #default BMIQ
champ.SVD(beta=myNorm,pd=myLoad$pd)
myCombat = champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c('Slide','Array'))
champ.SVD(beta=myCombat,pd=myLoad$pd)
saveRDS(myCombat,file='cobre_combat_final.rds')
