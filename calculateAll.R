# Running this R-file updates all tables and all figures of this research project
# The settings and randomization of data is prespecified in this file and all the 
# results generated are conditional on these settings.
# The other files contain merely functions to be called generating some output.


## =============================================================================
# Load data and necessary packages

library(excel.link)
library(rstatix)
library(knitr)
library(glmnet)
library(bnlearn)
library(mice)
library(rms)
library(VGAM)
library(DMwR)
library(dplyr)
library(tibble)
setwd("C:/Users/jef_m/OneDrive/Bureaublad/thesis MaStat")
df <- xl.read.file(filename="20231218_exportBMI_export.encrypt.xlsx", 
                   password=Sys.getenv("pw"))
df$ObesityClass <- c("Normal weight", "Overweight", "Obese")[1 + (df$BMI>25) + (df$BMI>30)]
df$ObesityClass <- factor(df$ObesityClass, levels=c("Normal weight", "Overweight", "Obese"))
df$Age <- df$`Maternal Age`
df$AgeGroup <- c("A", "B", "C", "D", "E")[1 + (df$`Maternal Age`>=20) + 
                                              (df$`Maternal Age`>=30) + 
                                              (df$`Maternal Age`>=35) + 
                                              (df$`Maternal Age`>=40)]
df$AgeGroup <- factor(df$AgeGroup, levels=c("A", "B", "C", "D", "E"))
df$ID <- sapply(1:nrow(df), sprintf, fmt="%04.0f")
df$Race <- factor(df$Race, levels=c("White", "Black", "South Asian", "East Asian", "Mixed"))

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

ethnicities <- c("White", "Black", "South Asian", "East Asian", "Mixed")
obesityClasses <- c("Normal weight", "Overweight", "Obese")
metabolites <- names(df)[-which(names(df)%in%c("Maternal Age","BMI", "Race", "Age", 
                                               "AgeGroup","Smoking status", "Control", 
                                               "ID", "ObesityClass", "met_025", "met_040", 
                                               "met_042", "met_082", "met_108"))]
#colors <- c("red","blue","green","violet","black")
#chars <- c("W", "B", "S", "E", "M")

data_model <- as.data.frame(log(makeX(train=df[,c("BMI", metabolites)], na.impute=TRUE)),
                            check.names=FALSE)
colnames(data_model) <- c("BMI", metabolites)
data_model <- data.frame(Race=df$Race, 
                         #Smoking=relevel(as.factor(df$`Smoking status`), ref="FALSE"), 
                         Smoking=as.numeric(df$`Smoking status`), 
                         Age=(df$`Maternal Age`-mean(df$`Maternal Age`))/sd(df$`Maternal Age`),
                         ObesityClass=df$ObesityClass,
                         data_model, check.names=FALSE)

source("convertToTexTable.R")
source("dataExploration.R")
source("linearRegression.R")

## =============================================================================
# Output results of the data exploration

#source("dataExploration.R")

# generate patient counts for every ethnicity and BMI class
biometrics <- tableBiometrics(df)
convertToTexTable(biometrics, "explore_biometrics.tex", 
                  caption="Sample patient count for every ethnicity and obesity class.", 
                  reflabel="sample_biometrics")

distributionList <- logNormalityAndOutliers(df)
convertToTexTable(distributionList$countOutliers, "explore_outliers.tex", 
                  caption="Number of outliers for every metabolite. Number of outliers that are also an outlier in other metabolites.", 
                  reflabel="sample_outliers")
convertToTexTable(distributionList$nonNormals, "explore_nonnormals.tex", 
                  caption="Listing of metabolites failing the test for log-normality. Cut-off p-value for the Shapiro-Wilk test is 10e-4.", 
                  reflabel="sample_nonnormals")

mets <- distributionList$nonNormals[,"met."]
pvals <- distributionList$nonNormals[,"p-val."]
makeQQplot(df, mets=mets, pvals=pvals)

demonstrateOutlierDetection(df, mets=c("met_002", "met_034"), ethnicity=c("White", "Black"))

correlations <- tableCorrelations(df)
convertToTexTable(correlations, "explore_correlations.tex", 
                  caption="Listing of the strongest correlations (>0.8)", 
                  reflabel="sample_correlations")

#boxplotMissing(df)
#df_complete <- boxplotMissing(df)

# the imputation technique was not found to perform better at predicting missing values than mean imputation
df_complete <- df

stratified_means_met002 <- df_complete %>% group_by(Race) %>% summarise(mean(log(met_002), na.rm=TRUE))
df_complete <- merge(df_complete, stratified_means_met002, by="Race")

stratified_means_met010 <- df_complete %>% group_by(Race) %>% summarise(mean(log(met_010), na.rm=TRUE))
df_complete <- merge(df_complete, stratified_means_met010, by="Race")

stratified_means_met068 <- df_complete %>% group_by(Race) %>% summarise(mean(log(met_068), na.rm=TRUE))
df_complete <- merge(df_complete, stratified_means_met068, by="Race")

df_complete$met_002[which(is.na(df_complete$met_002))] <- 
  exp(df_complete$`mean(log(met_002), na.rm = TRUE)`[which(is.na(df_complete$met_002))])
df_complete$met_010[which(is.na(df_complete$met_010))] <- 
  exp(df_complete$`mean(log(met_010), na.rm = TRUE)`[which(is.na(df_complete$met_010))])
df_complete$met_068[which(is.na(df_complete$met_068))] <- 
  exp(df_complete$`mean(log(met_068), na.rm = TRUE)`[which(is.na(df_complete$met_068))])

met10colnames <- names(df)[which(grepl(names(df), pattern="met_010"))]
w <- which(is.na(df$met_010))
for (name in met10colnames) {
  df_complete[[name]][w] <- with(df_complete[w,], eval(parse(text=name)))
}

df_complete$`mean(log(met_002), na.rm = TRUE)` <- NULL
df_complete$`mean(log(met_010), na.rm = TRUE)` <- NULL
df_complete$`mean(log(met_068), na.rm = TRUE)` <- NULL


# transform metabolites to log scale and transform other variables to a usable format
data_model <- as.data.frame(check.names=FALSE, 
                            log(makeX(train=df_complete[,c("BMI", metabolites)], 
                                      na.impute=FALSE)))
data_model <- data.frame(Race=relevel(as.factor(df_complete$Race), ref="White"), 
                         #Smoking=relevel(as.factor(df$`Smoking status`), ref="FALSE"), 
                         Smoking=as.numeric(df$`Smoking status`), 
                         Age=(df$`Maternal Age`-mean(df$`Maternal Age`))/sd(df$`Maternal Age`),
                         ObesityClass=relevel(as.factor(df$ObesityClass), 
                                              ref="Normal weight"),
                         data_model, check.names=FALSE)
colnames(data_model) <- c("Race", "Smoking", "Age", "ObesityClass", "BMI", metabolites)

# filter out outliers of data set
filtered_data <- filterOutliers(data_model, outliers=distributionList$outliers, 
                                mets=c("met_011", "met_015", "met_017", "met_027", 
                                       "met_029", "met_031", "met_034", "met_038",
                                       "met_049", "met_059", "met_060", "met_064", 
                                       "met_066", "met_073", "met_084", "met_085",
                                       "met_088", "met_090", "met_132", "met_134"))
data_regular <- filtered_data$regulars
data_outlier <- filtered_data$outliers
  

## =============================================================================
# Output results of the data exploration

#source("linearRegression.R")

# white ethnicity data

data_white <- subset(data_regular, subset = Race=="White")
set.seed(2)
trainSelectA <- sample(x=1:nrow(data_white), size=round(0.8*nrow(data_white)))
data_white_train <- data_white[trainSelectA,]
data_white_test <- data_white[-trainSelectA,]

data_white_train_balanced <- oversample(data_white_train)


# black ethnicity data

data_black <- subset(data_regular, subset = Race=="Black")
set.seed(2)
trainSelectB <- sample(x=1:nrow(data_black), size=round(0.8*nrow(data_black)))
data_black_train <- data_black[trainSelectB,]
data_black_test <- data_black[-trainSelectB,]


# retrieve or train OLS models or train if not available in workspace

effects1 <- c("main", "interaction", "main", "interaction")
types1 <- rep("OLS", times=4)
transformations1 <- rep(c("Log", "Inv"), each=2)
balancing1 <- rep("", times=4)

effects <- rep(c("main", "interaction"), times=4)
types <- rep("OLS", times=8)
transformations <- rep(c("Log", "Log", "Inv", "Inv"), times=2)
balancing <- rep(c("", "Balanced"), each=4)

filenames1 <- paste(effects, types, transformations, paste0("ModelWhite", balancing, ".rds"), sep="")
filenames2 <- paste(effects1, types1, transformations1, paste0("ModelBlack", ".rds"), sep="")
filenames <- c(filenames1, filenames2)
if (any(!file.exists(filenames))) {trainOLS(filenames[which(!file.exists(filenames))])}

mainOLSLogModelWhite <- readRDS("mainOLSLogModelWhite.rds")
interactionOLSLogModelWhite <- readRDS("interactionOLSLogModelWhite.rds")
mainOLSInvModelWhite <- readRDS("mainOLSInvModelWhite.rds")
interactionOLSInvModelWhite <- readRDS("interactionOLSInvModelWhite.rds")

mainOLSLogModelWhiteBalanced <- readRDS("mainOLSLogModelWhiteBalanced.rds")
interactionOLSLogModelWhiteBalanced <- readRDS("interactionOLSLogModelWhiteBalanced.rds")
mainOLSInvModelWhiteBalanced <- readRDS("mainOLSInvModelWhiteBalanced.rds")
interactionOLSInvModelWhiteBalanced <- readRDS("interactionOLSInvModelWhiteBalanced.rds")

mainOLSLogModelBlack <- readRDS("mainOLSLogModelBlack.rds")
interactionOLSLogModelBlack <- readRDS("interactionOLSLogModelBlack.rds")
mainOLSInvModelBlack <- readRDS("mainOLSInvModelBlack.rds")
interactionOLSInvModelBlack <- readRDS("interactionOLSInvModelBlack.rds")


# train ridge regression models

mainRidgeLogModelWhite <- trainRidge(data_white_train, transformation="Log", 
                                     interactionEffect=FALSE)
interactionRidgeLogModelWhite <- trainRidge(data_white_train, transformation="Log", 
                                            interactionEffect=TRUE)
mainRidgeInvModelWhite <- trainRidge(data_white_train, transformation="Inv", 
                                     interactionEffect=FALSE)
interactionRidgeInvModelWhite <- trainRidge(data_white_train, transformation="Inv", 
                                            interactionEffect=TRUE)

mainRidgeLogModelWhiteBalanced <- trainRidge(data_white_train_balanced, 
                                             transformation="Log", 
                                             interactionEffect=FALSE)
interactionRidgeLogModelWhiteBalanced <- trainRidge(data_white_train_balanced,
                                                    transformation="Log", 
                                                    interactionEffect=TRUE)
mainRidgeInvModelWhiteBalanced <- trainRidge(data_white_train_balanced,
                                             transformation="Inv", 
                                             interactionEffect=FALSE)
interactionRidgeInvModelWhiteBalanced <- trainRidge(data_white_train_balanced,
                                                    transformation="Inv", 
                                                    interactionEffect=TRUE)

mainRidgeLogModelBlack <- trainRidge(data_black_train, transformation="Log", 
                                     interactionEffect=FALSE)
interactionRidgeLogModelBlack <- trainRidge(data_black_train, transformation="Log", 
                                            interactionEffect=TRUE)
mainRidgeInvModelBlack <- trainRidge(data_black_train, transformation="Inv", 
                                     interactionEffect=FALSE)
interactionRidgeInvModelBlack <- trainRidge(data_black_train, transformation="Inv",
                                            interactionEffect=TRUE)


# train LASSO regression models

mainLASSOLogModelWhite <- trainLASSO(data_white_train, transformation="Log", 
                                     interactionEffect=FALSE)
interactionLASSOLogModelWhite <- trainLASSO(data_white_train, transformation="Log", 
                                            interactionEffect=TRUE)
mainLASSOInvModelWhite <- trainLASSO(data_white_train, transformation="Inv", 
                                     interactionEffect=FALSE)
interactionLASSOInvModelWhite <- trainLASSO(data_white_train, transformation="Inv", 
                                            interactionEffect=TRUE)

mainLASSOLogModelWhiteBalanced <- trainLASSO(data_white_train_balanced, 
                                             transformation="Log", 
                                             interactionEffect=FALSE)
interactionLASSOLogModelWhiteBalanced <- trainLASSO(data_white_train_balanced,
                                                    transformation="Log", 
                                                    interactionEffect=TRUE)
mainLASSOInvModelWhiteBalanced <- trainLASSO(data_white_train_balanced,
                                             transformation="Inv", 
                                             interactionEffect=FALSE)
interactionLASSOInvModelWhiteBalanced <- trainLASSO(data_white_train_balanced,
                                                    transformation="Inv", 
                                                    interactionEffect=TRUE)

mainLASSOLogModelBlack <- trainLASSO(data_black_train, transformation="Log", 
                                     interactionEffect=FALSE)
interactionLASSOLogModelBlack <- trainLASSO(data_black_train, transformation="Log", 
                                            interactionEffect=TRUE)
mainLASSOInvModelBlack <- trainLASSO(data_black_train, transformation="Inv", 
                                     interactionEffect=FALSE)
interactionLASSOInvModelBlack <- trainLASSO(data_black_train, transformation="Inv",
                                            interactionEffect=TRUE)



effects1 <- rep(c("main", "interaction", "main", "interaction", "main", "interaction"), times=2)
types1 <- rep(c("OLS", "OLS", "Ridge", "Ridge", "LASSO", "LASSO"), times=2)
transformations1 <- rep(c("Log", "Inv"), each=6)
balancing1 <- rep("", times=12)

effects2 <- rep(c("main", "interaction", "main", "interaction", "main", "interaction"), times=2)
types2 <- rep(c("OLS", "OLS", "Ridge", "Ridge", "LASSO", "LASSO"), times=2)
transformations2 <- rep(c("Log", "Inv"), each=6)
balancing2 <- rep("Balanced", times=12)

effects <- c(effects1, effects2)
types <- c(types1, types2)
transformations <- c(transformations1, transformations2)
balancing <- c(balancing1, balancing2)


# make table of model evaluations

trainWhite <- tabulatePredictionEvaluation(effects, types, transformations, balancing, 
                                           ethnicity="White", new.data=data_white_train)
trainBlack <- tabulatePredictionEvaluation(effects1, types1, transformations1, balancing1,
                                           ethnicity="Black", new.data=data_black_train)
convertToTexTable(trainWhite, "trainWhite.tex", 
                  caption="Prediction evaluation of the models on white ethnicity training data.", 
                  reflabel="trainWhite")
convertToTexTable(trainBlack, "trainBlack.tex", 
                  caption="Prediction evaluation of the models on black ethnicity training data.", 
                  reflabel="trainBlack")


# make posters of model diagnostics

modelDiagnostics(effects, types, transformations, balancing, 
                 ethnicity="White", new.data=data_white_train)
modelDiagnostics(effects1, types1, transformations1, 
                 ethnicity="Black", new.data=data_black_train)


# evaluate models with test set data

testWhite <- tabulatePredictionEvaluation(effects, types, transformations, balancing, 
                                          ethnicity="White", new.data=data_white_test)
testBlack <- tabulatePredictionEvaluation(effects1, types1, transformations1, balancing1,
                                          ethnicity="Black", new.data=data_black_test)
convertToTexTable(testWhite, "testWhite.tex", 
                  caption="Prediction evaluation of the models on white ethnicity test data.", 
                  reflabel="testWhite")
convertToTexTable(testBlack, "testBlack.tex", 
                  caption="Prediction evaluation of the models on black ethnicity test data.", 
                  reflabel="testBlack")

