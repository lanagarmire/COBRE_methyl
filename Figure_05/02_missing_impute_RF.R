library(mice)
library(ranger)
load('merged123.rdata')
impute<-mice(merged123,method='rf',seed=1)
impute_val<-complete(impute,5)
save(impute,file='imputed123.rdata')
write.csv(impute_val,'imputed123_val.csv',row.names=TRUE)

