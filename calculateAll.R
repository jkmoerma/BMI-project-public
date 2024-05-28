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
library(VGAM)
library(DMwR)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(foreach)
library(doParallel)
setwd("C:/Users/jef_m/OneDrive/Bureaublad/thesis MaStat")
df <- xl.read.file(filename="20231218_exportBMI_export.encrypt.xlsx", 
                   password=Sys.getenv("pw"))
df$ObesityClass <- c("Normal weight", "Overweight", "Obese")[1 + (df$BMI>25) + (df$BMI>30)]
df$ObesityClass <- factor(df$ObesityClass, levels=c("Normal weight", "Overweight", "Obese"))
df$Age <- df$`Maternal Age`
df$ID <- sapply(1:nrow(df), sprintf, fmt="%04.0f")
df$Race <- factor(df$Race, levels=c("White", "Black", "South Asian", "East Asian", "Mixed"))

relMet1 <- c("met_078", "met_010", "met_017", "met_093", "met_047", "met_012")
for (i in 1:(length(relMet1)-1)) {
  for (j in (i+1):length(relMet1)) {
    met1 <- relMet1[i]
    met2 <- relMet1[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- (df[[met1]]/df[[met2]])
  }
}
relMet2 <- c("met_038", "met_041", "met_032")
for (i in 1:(length(relMet2)-1)) {
  for (j in (i+1):length(relMet2)) {
    met1 <- relMet2[i]
    met2 <- relMet2[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- (df[[met1]]/df[[met2]])
  }
}
relMet3 <- c("met_012", "met_026", "met_038", "met_041", "met_093")
for (i in 1:(length(relMet3)-1)) {
  for (j in (i+1):length(relMet3)) {
    met1 <- relMet3[i]
    met2 <- relMet3[j]
    ratio <- paste0(met1, "/", met2)
    df[[ratio]] <- (df[[met1]]/df[[met2]])
  }
}
for (i in 1:length(relMet3)) {
  met1 <- relMet3[i]
  met2 <- "met_047"
  ratio <- paste0(met1, "/", met2)
  df[[ratio]] <- (df[[met1]]/df[[met2]])
}


ethnicities <- c("White", "Black", "South Asian", "East Asian", "Mixed")
obesityClasses <- c("Normal weight", "Overweight", "Obese")
metabolites <- names(df)[-which(names(df)%in%c("Maternal Age","BMI", "Race", "Age", 
                                               "Smoking status", "Control", "AgeGroup",
                                               "ID", "ObesityClass", "met_025", "met_040", 
                                               "met_042", "met_082", "met_108"))]

source("convertToTexTable.R")
source("dataExploration.R")
source("linearRegression.R")
source("metabolicAnalysis.R")

## =============================================================================
# Output results of the data exploration

#source("dataExploration.R")

# generate patient counts for every ethnicity and BMI class
biometrics <- tableBiometrics(df)
convertToTexTable(biometrics, "explore_biometrics.tex", rows.named=TRUE,
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
data_model <- df_complete
for (var in c("BMI", metabolites)) {
  data_model[[var]] <- log(df_complete[[var]])
}
data_model <- data_model %>% 
                mutate(Race=relevel(as.factor(Race), ref="White")) %>%
                mutate(Smoking=as.numeric(`Smoking status`)) %>%
                mutate(Age=(`Maternal Age`-mean(`Maternal Age`))/sd(`Maternal Age`)) %>%
                mutate(ObesityClass=relevel(as.factor(ObesityClass), ref="Normal weight"))
data_model$`Maternal Age` <- NULL
data_model$Control <- NULL
data_model$`Smoking status` <- NULL
data_model$met_025 <- NULL
data_model$met_040 <- NULL
data_model$met_042 <- NULL
data_model$met_082 <- NULL
data_model$met_108 <- NULL


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
# Output results of the regression analysis

#source("linearRegression.R")

# white ethnicity data

data_white <- subset(data_regular, subset = Race=="White")
set.seed(2)
testSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
data_white_test <- data_white[testSelectA,]
data_white_use <- data_white[-testSelectA,]
set.seed(3)
trainSelectA <- sample(x=1:nrow(data_white_use), size=round(0.75*nrow(data_white_use)))
data_white_train <- data_white_use[trainSelectA,]
data_white_valid <- data_white_use[-trainSelectA,]


data_white_train_balanced <- oversample(data_white_train)


# black ethnicity data

data_black <- subset(data_regular, subset = Race=="Black")
set.seed(2)
testSelectB <- sample(x=1:nrow(data_black), size=round(0.2*nrow(data_black)))
data_black_test <- data_black[testSelectB,]
data_black_use <- data_black[-testSelectB,]
set.seed(3)
trainSelectB <- sample(x=1:nrow(data_black_use), size=round(0.75*nrow(data_black_use)))
data_black_train <- data_black_use[trainSelectB,]
data_black_valid <- data_black_use[-trainSelectB,]


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

if (!file.exists("mainRidgeLogModelWhite.rds")) {
  mainRidgeLogModelWhite <- trainRidge(data_white_train, transformation="Log", 
                                       interactionEffect=FALSE)
  interactionRidgeLogModelWhite <- trainRidge(data_white_train, transformation="Log", 
                                              interactionEffect=TRUE)
  mainRidgeInvModelWhite <- trainRidge(data_white_train, transformation="Inv", 
                                       interactionEffect=FALSE)
  interactionRidgeInvModelWhite <- trainRidge(data_white_train, transformation="Inv", 
                                              interactionEffect=TRUE)
  
  saveRDS(mainRidgeLogModelWhite, "mainRidgeLogModelWhite.rds")
  saveRDS(interactionRidgeLogModelWhite, "interactionRidgeLogModelWhite.rds")
  saveRDS(mainRidgeInvModelWhite, "mainRidgeInvModelWhite.rds")
  saveRDS(interactionRidgeInvModelWhite, "mainRidgeInvModelWhite.rds")
}
mainRidgeLogModelWhite <- readRDS("mainRidgeLogModelWhite.rds")
interactionRidgeLogModelWhite <- readRDS("interactionRidgeLogModelWhite.rds")
mainRidgeInvModelWhite <- readRDS("mainRidgeInvModelWhite.rds")
interactionRidgeInvModelWhite <- readRDS("mainRidgeInvModelWhite.rds")

if (!file.exists("mainRidgeLogModelWhiteBalanced.rds")) {
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
  
  saveRDS(mainRidgeLogModelWhiteBalanced, "mainRidgeLogModelWhiteBalanced.rds")
  saveRDS(interactionRidgeLogModelWhiteBalanced, "interactionRidgeLogModelWhiteBalanced.rds")
  saveRDS(mainRidgeInvModelWhiteBalanced, "mainRidgeInvModelWhiteBalanced.rds")
  saveRDS(interactionRidgeInvModelWhiteBalanced, "interactionRidgeInvModelWhiteBalanced.rds")
}
mainRidgeLogModelWhiteBalanced <- readRDS("mainRidgeLogModelWhiteBalanced.rds")
interactionRidgeLogModelWhiteBalanced <- readRDS("interactionRidgeLogModelWhiteBalanced.rds")
mainRidgeInvModelWhiteBalanced <- readRDS("mainRidgeInvModelWhiteBalanced.rds")
interactionRidgeInvModelWhiteBalanced <- readRDS("interactionRidgeInvModelWhiteBalanced.rds")

if (!file.exists("mainRidgeLogModelBlack.rds")) {
  mainRidgeLogModelBlack <- trainRidge(data_black_train, transformation="Log", 
                                       interactionEffect=FALSE)
  interactionRidgeLogModelBlack <- trainRidge(data_black_train, transformation="Log", 
                                              interactionEffect=TRUE)
  mainRidgeInvModelBlack <- trainRidge(data_black_train, transformation="Inv", 
                                       interactionEffect=FALSE)
  interactionRidgeInvModelBlack <- trainRidge(data_black_train, transformation="Inv",
                                              interactionEffect=TRUE)
  
  saveRDS(mainRidgeLogModelBlack, "mainRidgeLogModelBlack.rds")
  saveRDS(interactionRidgeLogModelBlack, "interactionRidgeLogModelBlack.rds")
  saveRDS(mainRidgeInvModelBlack, "mainRidgeInvModelBlack.rds")
  saveRDS(interactionRidgeInvModelBlack, "interactionRidgeInvModelBlack.rds")
  
}
mainRidgeLogModelBlack <- readRDS("mainRidgeLogModelBlack.rds")
interactionRidgeLogModelBlack <- readRDS("interactionRidgeLogModelBlack.rds")
mainRidgeInvModelBlack <- readRDS("mainRidgeInvModelBlack.rds")
interactionRidgeInvModelBlack <- readRDS("interactionRidgeInvModelBlack.rds")


# train LASSO regression models

if (!file.exists("mainLASSOLogModelWhite.rds")) {
  mainLASSOLogModelWhite <- trainLASSO(data_white_train, transformation="Log", 
                                       interactionEffect=FALSE)
  interactionLASSOLogModelWhite <- trainLASSO(data_white_train, transformation="Log", 
                                              interactionEffect=TRUE)
  mainLASSOInvModelWhite <- trainLASSO(data_white_train, transformation="Inv", 
                                       interactionEffect=FALSE)
  interactionLASSOInvModelWhite <- trainLASSO(data_white_train, transformation="Inv", 
                                              interactionEffect=TRUE)
  
  saveRDS(mainLASSOLogModelWhite, "mainLASSOLogModelWhite.rds")
  saveRDS(interactionLASSOLogModelWhite, "interactionLASSOLogModelWhite.rds")
  saveRDS(mainLASSOInvModelWhite, "mainLASSOInvModelWhite.rds")
  saveRDS(interactionLASSOInvModelWhite, "interactionLASSOInvModelWhite.rds")
}
mainLASSOLogModelWhite <- readRDS("mainLASSOLogModelWhite.rds")
interactionLASSOLogModelWhite <- readRDS("interactionLASSOLogModelWhite.rds")
mainLASSOInvModelWhite <- readRDS("mainLASSOInvModelWhite.rds")
interactionLASSOInvModelWhite <- readRDS("interactionLASSOInvModelWhite.rds")

if (!file.exists("mainLASSOLogModelWhiteBalanced.rds")) {
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
  
  saveRDS(mainLASSOLogModelWhiteBalanced, "mainLASSOLogModelWhiteBalanced.rds")
  saveRDS(interactionLASSOLogModelWhiteBalanced, "interactionLASSOLogModelWhiteBalanced.rds")
  saveRDS(mainLASSOInvModelWhiteBalanced, "mainLASSOInvModelWhiteBalanced.rds")
  saveRDS(interactionLASSOInvModelWhiteBalanced, "interactionLASSOInvModelWhiteBalanced.rds")
}
mainLASSOLogModelWhiteBalanced <- readRDS("mainLASSOLogModelWhiteBalanced.rds")
interactionLASSOLogModelWhiteBalanced <- readRDS("interactionLASSOLogModelWhiteBalanced.rds")
mainLASSOInvModelWhiteBalanced <- readRDS("mainLASSOInvModelWhiteBalanced.rds")
interactionLASSOInvModelWhiteBalanced <- readRDS("interactionLASSOInvModelWhiteBalanced.rds")

if (!file.exists("mainLASSOLogModelBlack.rds")) {
  mainLASSOLogModelBlack <- trainLASSO(data_black_train, transformation="Log", 
                                       interactionEffect=FALSE)
  interactionLASSOLogModelBlack <- trainLASSO(data_black_train, transformation="Log", 
                                              interactionEffect=TRUE)
  mainLASSOInvModelBlack <- trainLASSO(data_black_train, transformation="Inv", 
                                       interactionEffect=FALSE)
  interactionLASSOInvModelBlack <- trainLASSO(data_black_train, transformation="Inv", 
                                              interactionEffect=TRUE)
  
  saveRDS(mainLASSOLogModelBlack, "mainLASSOLogModelBlack.rds")
  saveRDS(interactionLASSOLogModelBlack, "interactionLASSOLogModelBlack.rds")
  saveRDS(mainLASSOInvModelBlack, "mainLASSOInvModelBlack.rds")
  saveRDS(interactionLASSOInvModelBlack, "interactionLASSOInvModelBlack.rds")
}
mainLASSOLogModelBlack <- readRDS("mainLASSOLogModelBlack.rds")
interactionLASSOLogModelBlack <- readRDS("interactionLASSOLogModelBlack.rds")
mainLASSOInvModelBlack <- readRDS("mainLASSOInvModelBlack.rds")
interactionLASSOInvModelBlack <- readRDS("interactionLASSOInvModelBlack.rds")



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


# choose models with validation set data

validWhite <- tabulatePredictionEvaluation(effects, types, transformations, balancing, 
                                          ethnicity="White", new.data=data_white_valid)
validBlack <- tabulatePredictionEvaluation(effects1, types1, transformations1, balancing1,
                                          ethnicity="Black", new.data=data_black_valid)
convertToTexTable(validWhite, "validWhite.tex", 
                  caption="Prediction evaluation of the models on white ethnicity validation data.", 
                  reflabel="validWhite")
convertToTexTable(validBlack, "validBlack.tex", 
                  caption="Prediction evaluation of the models on black ethnicity validation data.", 
                  reflabel="validBlack")

plotPredictions(effects, types, transformations, balancing, ethnicity="White", 
                chosen=23, new.data=data_white_valid)
plotPredictions(effects1, types1, transformations1, balancing1, ethnicity="Black", 
                chosen=9, new.data=data_black_valid)


# make posters of model diagnostics of chosen models

modelDiagnostics(effects="main", types="LASSO", transformations="Inv", balancing="Balanced", 
                 ethnicity="White", new.data=data_white_train)
modelDiagnostics(effects="main", types="Ridge", transformations="Inv", 
                 ethnicity="Black", new.data=data_black_train)


# evaluate final models with test set data

testWhite <- tabulatePredictionEvaluation(effects="main", types="LASSO", 
                                          transformations="Inv", balancing="Balanced", 
                                          ethnicity="White", new.data=data_white_test)
testBlack <- tabulatePredictionEvaluation(effects="main", types="Ridge", 
                                          transformations="Inv", balancing="",
                                          ethnicity="Black", new.data=data_black_test)
convertToTexTable(testWhite, "testWhite.tex", 
                  caption="Prediction evaluation of the models on white ethnicity test data.", 
                  reflabel="testWhite")
convertToTexTable(testBlack, "testBlack.tex", 
                  caption="Prediction evaluation of the models on black ethnicity test data.", 
                  reflabel="testBlack")

testPredsWhite <- predictRidgeLASSO(data_white_test, mainLASSOInvModelWhiteBalanced, interactionEffect=FALSE)
testPredsClassWhite <- c("predicted: Normal weight", "predicted: Overweight", 
                         "predicted: Obese") [1 + (testPredsWhite<1/25) + (testPredsWhite<1/30)]
testPredsClassWhite <- factor(testPredsClassWhite, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
testObsClassWhite <- data_white_test$ObesityClass
testPredictTableWhite <- table(testPredsClassWhite, testObsClassWhite)
convertToTexTable(testPredictTableWhite, "testPredictTableWhite.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for white ethnicity patients.",
                  reflabel="testPredictTableWhite")

testPredsBlack <- predictRidgeLASSO(data_black_test, mainRidgeInvModelBlack, interactionEffect=FALSE)
testPredsClassBlack <- c("predicted: Normal weight", "predicted: Overweight", 
                         "predicted: Obese") [1 + (testPredsBlack<1/25) + (testPredsBlack<1/30)]
testPredsClassBlack <- factor(testPredsClassBlack, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
testObsClassBlack <- data_black_test$ObesityClass
testPredictTableBlack <- table(testPredsClassBlack, testObsClassBlack)
convertToTexTable(testPredictTableBlack, "testPredictTableBlack.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for black ethnicity patients.",
                  reflabel="testPredictTableBlack")


## =============================================================================
# Output results of metabolic outlier analysis

#source("metabolicAnalysis.R")


# calculate power and sample size of fitted metabolite beta coefficients
finalWhite <- trainLASSO(oversample(data_white), transformation="Inv")
finalBlack <- trainRidge(data_black, transformation="Inv")

whitePredictions <- 1/predict(finalWhite, newx=as.matrix(data_white[,c(metabolites, "Age", "Smoking")]))[,"s0"]
blackPredictions <- 1/predict(finalBlack, newx=as.matrix(data_black[,c(metabolites, "Age", "Smoking")]))[,"s0"]

if (any(!file.exists(c("whiteSimulation.rds", "blackSimulation.rds")))) {
  whiteSimulation <- simulateAlternative(finalWhite, data_white, 1/whitePredictions,
                                         transformation="Inv", type="LASSO", nsim=1000, oversampled=TRUE)
  blackSimulation <- simulateAlternative(finalBlack, data_black, 1/blackPredictions,
                                         transformation="Inv", type="Ridge", nsim=1000)
  
  saveRDS(whiteSimulation, "whiteSimulation.rds")
  saveRDS(blackSimulation, "blackSimulation.rds")
}

whiteSimulation <- readRDS("whiteSimulation.rds")
blackSimulation <- readRDS("blackSimulation.rds")


# calculate power from the standard deviations of the simulated coefficients

whitePower <- plotPower(whiteSimulation$pvals, plottitle="Power of predictors (white ethnicity)")
blackPower <- plotPower(blackSimulation$pvals, plottitle="Power of predictors (black ethnicity)")

ggsave(filename="whitePower.pdf", whitePower)
ggsave(filename="blackPower.pdf", blackPower)


# check normality of beta coefficients from the simulation under the alternative

simBetasWhite <- pivot_longer(whiteSimulation$betas, cols=colnames(whiteSimulation$betas), 
                              names_to="met", values_to="beta")
standardize <- "(beta-mean(beta, na.rm=TRUE))/sd(beta)"
Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")
simBetasWhite <- simBetasWhite %>% group_by(met) %>% 
  mutate(standardized=eval(parse(text=standardize)),
         normQuantile=eval(parse(text=Zscore)))
plotSimBetasWhite <- ggplot(data=simBetasWhite, aes(x=normQuantile, y=standardized, by=met)) +
                       geom_line(aes(color=met)) + 
                       labs(title = "QQ-plots of simulation beta coefficiets (White)",
                            x = "Normal quantile", y = "Standardized beta coefficient")

simBetasBlack <- pivot_longer(blackSimulation$betas, cols=colnames(blackSimulation$betas), 
                              names_to="met", values_to="beta")
standardize <- "(beta-mean(beta, na.rm=TRUE))/sd(beta)"
Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")
simBetasBlack <- simBetasBlack %>% group_by(met) %>% 
  mutate(standardized=eval(parse(text=standardize)),
         normQuantile=eval(parse(text=Zscore)))
plotSimBetasBlack <- ggplot(data=simBetasBlack, aes(x=normQuantile, y=standardized, by=met)) +
                       geom_line(aes(color=met)) + 
                       labs(title = "QQ-plots of simulation beta coefficiets (Black)",
                            x = "Normal quantile", y = "Standardized beta coefficient")

ggsave(filename="simBetasWhiteQQplot.pdf", plotSimBetasWhite)
ggsave(filename="simBetasBlackQQplot.pdf", plotSimBetasBlack)


# use this assumption to perform a sample size calculation



# perform hypothesis test for patients with similar predicted BMI
if (!file.exists("TukeyCorrected.rds")) {
  TukeyCorrected <- matrix(nrow=2, ncol=3)
  colnames(TukeyCorrected) <- c("22-25", "26-29", "30-35")
  rownames(TukeyCorrected) <- c("White", "Black")
  for (range in colnames(TukeyCorrected)) {
    
    splitted_range <- strsplit(range, split="-")[[1]]
    lower <- as.numeric(splitted_range[1])
    upper <- as.numeric(splitted_range[2])
    
    use_data_white <- subsetMetabolicBands(data_white, whitePredictions, 
                                           band_lower=lower, band_upper=upper)
    TukeyCorrected["White", range] <- multipleCorrTest(use_data_white, nsim0=1000)
    
    use_data_black <- subsetMetabolicBands(data_black, blackPredictions, 
                                           band_lower=lower, band_upper=upper)
    TukeyCorrected["Black", range] <- multipleCorrTest(use_data_black, nsim0=1000)
  }
  saveRDS(TukeyCorrected, "TukeyCorrected.rds")
}
TukeyCorrected <- readRDS("TukeyCorrected.rds")
convertToTexTable(TukeyCorrected, "TukeyCorrected.tex", rows.named=TRUE,
                  caption="Tukey Corrected p-values for hypothesizing that patients with the same predicted BMI also have the same metabolite profile.", 
                  reflabel="TukeyCorrected")


# Calculate power for increasing sample sizes
if (!file.exists("TukeyCorrectedPowerWhite.rds")&file.exists("TukeyCorrectedPowerBlack.rds")) {
  
  ranges <- c("22-25", "26-29", "30-35")
  n_values <- c(50, 100, 200, 500)
  
  TukeyCorrectedPowerWhite <-
    data.frame(stratum = rep("White", each=length(ranges)*length(n_values)),
               range = rep(ranges, each=length(n_values)),
               n = rep(n_values, times=length(ranges)))
  
  TukeyCorrectedPowerBlack <-
    data.frame(stratum = rep("Black", each=length(ranges)*length(n_values)),
               range = rep(ranges, each=length(n_values)),
               n = rep(n_values, times=length(ranges)))
  
  numCores <- detectCores()
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  results <- foreach(i = 1:nrow(TukeyCorrectedPowerBlack), .combine='rbind') %dopar% {
    range <- TukeyCorrectedPowerBlack$range[i]
    n <- TukeyCorrectedPowerBlack$n[i]
    
    splitted_range <- strsplit(range, split="-")[[1]]
    lower <- as.numeric(splitted_range[1])
    upper <- as.numeric(splitted_range[2])
    
    use_data_white <- subsetMetabolicBands(data_white, whitePredictions, 
                                           band_lower=lower, band_upper=upper)
    power_white <- powerCorrTest(use_data_white, n, nsim0=200, nsimP=100)
    
    use_data_black <- subsetMetabolicBands(data_black, blackPredictions, 
                                           band_lower=lower, band_upper=upper)
    power_black <- powerCorrTest(use_data_black, n, nsim0=200, nsimP=100)
    
    return(c(power_white, power_black))
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # write the results to data frames
  TukeyCorrectedPowerWhite$power <- results[,1]
  TukeyCorrectedPowerBlack$power <- results[,2]
  
  # save direly calculated results
  saveRDS(TukeyCorrectedPowerWhite, "TukeyCorrectedPowerWhite.rds")
  saveRDS(TukeyCorrectedPowerBlack, "TukeyCorrectedPowerBlack.rds")
  
}

TukeyCorrectedPowerWhite <- readRDS("TukeyCorrectedPowerWhite.rds")
TukeyCorrectedPowerBlack <- readRDS("TukeyCorrectedPowerBlack.rds")



## =============================================================================
# Output results of discussion

#source("discussion.R")







colors <- c("grey", "red", "green", "blue")[1+(whitePredictions>22&whitePredictions<25)+
                                              2*(whitePredictions>26&whitePredictions<29)+
                                              3*(whitePredictions>30&whitePredictions<35)]
plot(whitePredictions, exp(data_white$BMI),
     col=colors)
abline(v=c(22, 25), col="red", lty="dashed")
text(x=23.5, y=45, labels=length(which(colors=="red")), col="red")
abline(v=c(26, 29), col="green", lty="dashed")
text(x=27.5, y=45, labels=length(which(colors=="green")), col="green")
abline(v=c(30, 35), col="blue", lty="dashed")
text(x=32.5, y=45, labels=length(which(colors=="blue")), col="blue")

colors <- c("grey", "red", "green", "blue")[1+(blackPredictions>22&blackPredictions<25)+
                                              2*(blackPredictions>26&blackPredictions<29)+
                                              3*(blackPredictions>30&blackPredictions<35)]
plot(blackPredictions, exp(data_black$BMI),
     col=colors)
abline(v=c(22, 25), col="red", lty="dashed")
text(x=23.5, y=45, labels=length(which(colors=="red")), col="red")
abline(v=c(26, 29), col="green", lty="dashed")
text(x=27.5, y=45, labels=length(which(colors=="green")), col="green")
abline(v=c(30, 35), col="blue", lty="dashed")
text(x=32.5, y=45, labels=length(which(colors=="blue")), col="blue")
