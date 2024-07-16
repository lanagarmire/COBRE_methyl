library("ChAMP")
library(parallel)
library(car)
library("lumi")
library("dplyr")
library("AnalyzeFMRI")


myCombat = readRDS('cobre/cobre_combat_final.rds')
mypd = readRDS('cobre/cobre_pd.rds')

sovdat = data.frame(t(myCombat))
pd = data.frame(mypd)[,c(3,6:14)]

sovdatlist <- list()
for(i in 1:ncol(sovdat)){
    sovdatlist[[i]] <- sovdat[,i]
}

Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)

for(i in 2:ncol(pd)){
    varname <- names(pd)[i]
    Ftab[varname] <- numeric()
}

calF <- function(probe = probecol){
    newdata <- pd
    pdnames <- names(newdata)
    newdata$beta <- probe
    
    formstr <- paste0(pdnames, collapse = ' + ')
    formstr <- paste0('beta ~ ', formstr)
    formstr <- as.formula(formstr)
    
    fit <- lm(formstr, data = newdata)
    
    aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
    
    
    F <- aovfit$`F value`
    
    F <- F[2:(length(F)-1)]
    names(F) <- pdnames
    F <- as.data.frame(F, stringsAsFactors = FALSE)
    F <- as.data.frame(t(F))
    row.names(F) <- 1
    
    
    Ftab <- rbind(Ftab, F)
    
    return(Ftab)
}

Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 40)


Ftab <- do.call(rbind, Ftab)


Fmean <- colMeans(Ftab)

Fmean <- Fmean[order(-Fmean)]

Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)

finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))


save(Ftab, Fmean, finalvars, file = "cobre_final_res.rds")




