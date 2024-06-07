# This file is a validation exercise for the imputation techniques used in this project
# Comments are that imputing the variable means for the missing values may corrupt the analysis.
# In this file it is explored whether existing data can be imputed more correctly by the technique 
# than by guessing the missing values as means.

# From the complete cases, values of a metabolite will be erased at 4.17% per iteration and 
# predicted by the imputation algorithm. The deviations for the predictions will be compared 
# with the deviations from mean

library(doParallel)
library(parallel)
library(foreach)
library(dplyr)

setwd("C:/Users/jef_m/OneDrive/Bureaublad/thesis MaStat")
df <- xl.read.file(filename="20231218_exportBMI_export.encrypt.xlsx", 
                   password=Sys.getenv("pw"))
df$ObesityClass <- c("Normal weight", "Overweight", "Obese")[1 + (df$BMI>25) + (df$BMI>30)]
df$Age <- df$`Maternal Age`
df$AgeGroup <- as.factor(c("A", "B", "C", "D", "E")[1 + (df$`Maternal Age`>=20) + 
                                                      (df$`Maternal Age`>=30) + 
                                                      (df$`Maternal Age`>=35) + 
                                                      (df$`Maternal Age`>=40)])
df$AgeGroup <- relevel(df$AgeGroup, ref="B")
df$ID <- sapply(1:nrow(df), sprintf, fmt="%04.0f")
df$Race <- as.factor(df$Race)
df$Race <- relevel(df$Race, ref="White")

relMet1 <- c("met_078", "met_010", "met_017", "met_093", "met_047", "met_012")
for (i in 1:(length(relMet1)-1)) {
  for (j in (i+1):length(relMet1)) {
    met1 <- relMet1[i]
    met2 <- relMet1[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- df[[met1]]/df[[met2]]
  }
}
relMet2 <- c("met_038", "met_041", "met_032")
for (i in 1:(length(relMet2)-1)) {
  for (j in (i+1):length(relMet2)) {
    met1 <- relMet2[i]
    met2 <- relMet2[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- df[[met1]]/df[[met2]]
  }
}
relMet3 <- c("met_012", "met_026", "met_038", "met_041", "met_093")
for (i in 1:(length(relMet3)-1)) {
  for (j in (i+1):length(relMet3)) {
    met1 <- relMet3[i]
    met2 <- relMet3[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- df[[met1]]/df[[met2]]
  }
}
for (i in 1:length(relMet3)) {
  met1 <- relMet3[i]
  met2 <- "met_047"
  ratio <- paste0(met1, "/", met2)
  df[[ratio]] <- df[[met1]]/df[[met2]]
}

source("dataExploration.R")

complete_cases <- subset(df, subset = !is.na(met_002)&!is.na(met_010)&!is.na(met_068))

if (FALSE) {
  set.seed(132)
  foldid <- sample(rep(1:24, each=65), size=nrow(complete_cases), replace=FALSE)
  
  # Detect the number of available cores and create cluster
  cl <- parallel::makeCluster(detectCores())
  # Activate cluster for foreach library
  doParallel::registerDoParallel(cl)
  
  imputeErrors <- foreach::foreach(i=1:24, .combine=rbind,
                                   .packages=c("glmnet", "dplyr", "mice")) %dopar% {
    
    w <- which(foldid==i)
    
    
    ### met_002
    
    # erase a part of the complete data
    erase_met002 <- complete_cases
    erase_met002$met_002[w] <- NA
    IDs <- erase_met002$ID[w]
    
    # repopulate the missing data with the imputation algorithm
    repopulated002 <- boxplotMissing(erase_met002, makePlot=FALSE)
    repopulated002 <- data.frame("ID"=repopulated002$ID, "imp_met_002"=repopulated002$met_002)
    erase_met002 <- merge(erase_met002, repopulated002, by="ID")
    
    # calculate the stratified mean concentration which is the alternative
    stratified_means_met002 <- erase_met002 %>% group_by(Race) %>% summarise(mean(log(met_002), na.rm=TRUE))
    erase_met002 <- merge(erase_met002, stratified_means_met002, by="Race")
    
    # collect deviations for evaluation
    erase_met002 <- merge(erase_met002, complete_cases, by="ID")
    v <- which(erase_met002$ID %in% IDs)
    imputeErrors002 <- log(erase_met002$imp_met_002[v]) - log(erase_met002$met_002.y[v])
    meanErrors002 <- erase_met002$`mean(log(met_002), na.rm = TRUE)`[v] - log(erase_met002$met_002.y[v])
    
    
    
    ### met_068
    
    # erase a part of the complete data
    erase_met068 <- complete_cases
    erase_met068$met_068[w] <- NA
    IDs <- erase_met068$ID[w]
    
    # repopulate the missing data with the imputation algorithm
    repopulated068 <- boxplotMissing(erase_met068, makePlot=FALSE)
    repopulated068 <- data.frame("ID"=repopulated068$ID, "imp_met_068"=repopulated068$met_068)
    erase_met068 <- merge(erase_met068, repopulated068, by="ID")
    
    # calculate the stratified mean concentration which is the alternative
    stratified_means_met068 <- erase_met068 %>% group_by(Race) %>% summarise(mean(log(met_068), na.rm=TRUE))
    erase_met068 <- merge(erase_met068, stratified_means_met068, by="Race")
    
    # collect deviations for evaluation
    erase_met068 <- merge(erase_met068, complete_cases, by="ID")
    v <- which(erase_met068$ID %in% IDs)
    imputeErrors068 <- log(erase_met068$imp_met_068[v]) - log(erase_met068$met_068.y[v])
    meanErrors068 <- erase_met068$`mean(log(met_068), na.rm = TRUE)`[v] - log(erase_met068$met_068.y[v])
    
    
    
    cbind(meanErrors002, imputeErrors002, 
          meanErrors068, imputeErrors068)
  }
  
  # Stop cluster to free up resources
  parallel::stopCluster(cl)
  
  saveRDS(as.data.frame(imputeErrors), "imputeErrors.rds")
  
}

imputeErrors <- readRDS("imputeErrors.rds")

imputeErrors %>% summarise(mean(meanErrors002), mean(imputeErrors002))
imputeErrors %>% summarise(mean(meanErrors068), mean(imputeErrors068))

imputeErrors %>% summarise(sd(meanErrors002), sd(imputeErrors002))
imputeErrors %>% summarise(sd(meanErrors068), sd(imputeErrors068))

imputeErrors %>% summarise(cor(meanErrors002, imputeErrors002))
imputeErrors %>% summarise(cor(meanErrors068, imputeErrors068))

plot(imputeErrors002~meanErrors002, data=imputeErrors)
abline(a=0, b=1, lty="dashed")

plot(imputeErrors068~meanErrors068, data=imputeErrors)
abline(a=0, b=1, lty="dashed")

plot(log(complete_cases$met_002), imputeErrors$meanErrors002, col="red", main="met_002")
points(log(complete_cases$met_002), imputeErrors$imputeErrors002, col="blue")
legend("topright", legend=c("mean imp.", "MICE imp."), col=c("red", "blue"), pch=1)

plot(log(complete_cases$met_068), imputeErrors$meanErrors068, col="red", main="met_068")
points(log(complete_cases$met_068), imputeErrors$imputeErrors068, col="blue")
legend("topright", legend=c("mean imp.", "MICE imp."), col=c("red", "blue"), pch=1)


