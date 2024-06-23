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
library(pROC)
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
source("reclassify.R")

## =============================================================================
# Output results of the data exploration

#source("dataExploration.R")

# generate patient counts for every ethnicity and BMI class
biometrics <- tableBiometrics(df)
convertToTexTable(biometrics, "explore_biometrics.tex", rows.named=TRUE,
                  caption="Sample patient count for every ethnicity and obesity class.", 
                  reflabel="sample_biometrics")

# tabulate number of smoking patients 
smokingCounts <- tableSmokingStatus(df)
convertToTexTable(smokingCounts, "explore_smoking.tex",
                  caption="Number of smoking and non-smoking patients for every ethnicity.",
                  reflabel="explore_smoking")

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

data_white <- subset(data_regular, subset = Race=="White"&(!Smoking))
data_white_smoking <- subset(data_regular, subset = Race=="White"&Smoking)
set.seed(2)
valSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
data_white_val <- data_white[valSelectA,]
data_white_train <- data_white[-valSelectA,]


data_white_train_balanced <- oversample(data_white_train)


# black ethnicity data

data_black <- subset(data_regular, subset = Race=="Black"&(!Smoking))
data_black_smoking <- subset(data_regular, subset = Race=="Black"&Smoking)
set.seed(2)
valSelectB <- sample(x=1:nrow(data_black), size=round(0.2*nrow(data_black)))
data_black_val <- data_black[valSelectB,]
data_black_train <- data_black[-valSelectB,]



# tabulate cross-validated performance of all model recipes

formulations <- expand.grid(effects=c("main", "interaction"), 
                            types=c("OLS", "Ridge", "LASSO"), 
                            transformations=c("Log", "Inv"),
                            stringsAsFactors=FALSE)
balancing <- rep("Balanced", times=12)

if (!file.exists("trainWhite.tex")) {
  trainWhite <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations, ethnicity="White", 
                                   balancing, data_white_train)
  trainBlack <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations, ethnicity="Black", 
                                   new.data=data_black_train)
  convertToTexTable(trainWhite, "trainWhite.tex", 
                    caption="Cross-validated prediction evaluation of the models on white ethnicity training data.", 
                    reflabel="trainWhite")
  convertToTexTable(trainBlack, "trainBlack.tex", 
                    caption="Cross-validated prediction evaluation of the models on black ethnicity training data.", 
                    reflabel="trainBlack")
}


modelWhiteFull <- 
  trainRidgeLASSO(effect="interaction", type="Ridge", transformation="Inv", new.data=oversample(data_white))
modelBlackFull <- 
  trainRidgeLASSO(effect="main", type="Ridge", transformation="Inv", new.data=data_black)

data_white$predicted <- 1/predict(object=modelWhiteFull, 
                                  newx=makeMatrix(data_white, includeInteraction=TRUE)$mat,
                                  type="response")[,"s0"]
allPredsWhite <- data_white$predicted
allPredsClassWhite <- c("predicted: Normal weight", "predicted: Overweight", 
                        "predicted: Obese") [1 + (allPredsWhite>25) + (allPredsWhite>30)]
allPredsClassWhite <- factor(allPredsClassWhite, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
allObsClassWhite <- data_white$ObesityClass
allPredictTableWhite <- table(allPredsClassWhite, allObsClassWhite)
convertToTexTable(allPredictTableWhite, "allPredictTableWhite.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for all white ethnicity patients.",
                  reflabel="allPredictTableWhite")

data_black$predicted <- 1/predict(object=modelBlackFull, 
                           newx=makeMatrix(data_black, includeInteraction=FALSE)$mat,
                           type="response")[,"s0"]
allPredsBlack <- data_black$predicted
allPredsClassBlack <- c("predicted: Normal weight", "predicted: Overweight", 
                        "predicted: Obese") [1 + (allPredsBlack>25) + (allPredsBlack>30)]
allPredsClassBlack <- factor(allPredsClassBlack, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
allObsClassBlack <- data_black$ObesityClass
allPredictTableBlack <- table(allPredsClassBlack, allObsClassBlack)
convertToTexTable(allPredictTableBlack, "allPredictTableBlack.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for all black ethnicity patients.",
                  reflabel="allPredictTableBlack")


# make posters of model diagnostics of chosen models

modelWhiteTrain <- 
  trainRidgeLASSO(effect="interaction", type="Ridge", transformation="Inv", 
                  new.data=oversample(data_white_train))
modelBlackTrain <- 
  trainRidgeLASSO(effect="main", type="Ridge", transformation="Inv", 
                  new.data=data_black_train)

modelDiagnostics(effects="interaction", types="Ridge", transformations="Inv", balancing="Balanced", 
                 ethnicity="White", model=modelWhiteTrain, new.data=data_white_train)
modelDiagnostics(effects="main", types="Ridge", transformations="Inv", 
                 ethnicity="Black", model=modelBlackTrain, new.data=data_black_train)





## =============================================================================
# Output results of analysis on metabolic BMI

#source("metabolicAnalysis.R")

# Calculate scaled model coefficients + 95% CI from bootstrap replication

if (!file.exists("CIcoeffsWhite.pdf")) {
  
  bootstrapCoeffsWhite <- scaledEffects(data=data_white, effect="interaction", 
                                        type="Ridge", transformation="Inv", 
                                        balancing="Balanced", boot.n=200)
  bootstrapCoeffsWhite <- pivot_longer(bootstrapCoeffsWhite, 
                                       cols=colnames(bootstrapCoeffsWhite),
                                       names_to="met", values_to="coeff")
  CIcoeffsWhite <- bootstrapCoeffsWhite %>% 
    group_by(met) %>% summarise("est."=quantile(x=coeff, probs=0.5),
                                "95pc lCI"=quantile(x=coeff, probs=0.025),
                                "95pc uCI"=quantile(x=coeff, probs=0.975))
  CIcoeffsWhite <- CIcoeffsWhite[order(abs(CIcoeffsWhite$est.), decreasing=TRUE),]
  CIcoeffsWhite["direction"] <- ifelse(CIcoeffsWhite$est.>0,"Pos","Neg")
  CIcoeffsWhite <- CIcoeffsWhite[1:50,]
  p_CIcoeffsWhite <- ggplot(CIcoeffsWhite,aes(x=reorder(met,abs(est.)), y=est., 
                                              ymin=`95pc lCI`, ymax=`95pc uCI`, 
                                              col=direction)) +
    geom_pointrange(stat="identity") +
    geom_abline(slope=0, intercept=0, lty="dashed") +
    theme(axis.text.x = element_text(angle = 90),
          legend.title = element_blank(), legend.position="none") +
    xlab("") + ylab("coefficient") +
    coord_flip() +
    labs(title="White ethnicity: coefficients + 95% CI")
  
  bootstrapCoeffsBlack <- scaledEffects(data=data_black, effect="main", 
                                        type="Ridge", transformation="Inv", 
                                        balancing="", boot.n=200)
  bootstrapCoeffsBlack <- pivot_longer(bootstrapCoeffsBlack, 
                                       cols=colnames(bootstrapCoeffsBlack),
                                       names_to="met", values_to="coeff")
  CIcoeffsBlack <- bootstrapCoeffsBlack %>% 
    group_by(met) %>% summarise("est."=quantile(x=coeff, probs=0.5),
                                "95pc lCI"=quantile(x=coeff, probs=0.025),
                                "95pc uCI"=quantile(x=coeff, probs=0.975))
  CIcoeffsBlack <- CIcoeffsBlack[order(abs(CIcoeffsBlack$est.), decreasing=TRUE),]
  CIcoeffsBlack["direction"] <- ifelse(CIcoeffsBlack$est.>0,"Pos","Neg")
  CIcoeffsBlack <- CIcoeffsBlack[1:50,]
  p_CIcoeffsBlack <- ggplot(CIcoeffsBlack,aes(x=reorder(met,abs(est.)), y=est., 
                                              ymin=`95pc lCI`, ymax=`95pc uCI`, 
                                              col=direction))+
    geom_pointrange(stat="identity") +
    geom_abline(slope=0, intercept=0, lty="dashed") +
    theme(axis.text.x = element_text(angle = 90),
          legend.title = element_blank(), legend.position="none") +
    xlab("") + ylab("coefficient") +
    coord_flip() +
    labs(title="Black ethnicity: coefficients + 95% CI")
  
  ggsave("CIcoeffsWhite.pdf", p_CIcoeffsWhite)
  ggsave("CIcoeffsBlack.pdf", p_CIcoeffsBlack)
  
}


# Calculation of predictions of all patients.
# Predict 1/4 of patients with a model trained on the other 3/4

sampleID <- sample(1:4, size=nrow(data_white), replace=TRUE)
predict_white <- data.frame(ID=c(), predicted=c())
for (i in 1:4) {
  w <- which(sampleID==i)
  modelWhiteCross <- 
    trainRidgeLASSO(effect="interaction", type="Ridge", transformation="Inv", 
                    new.data=oversample(data_white[-w,]))
  whitePredictionsCross <- 
    1/predict(modelWhiteCross, newx=makeMatrix(data_white[w,], includeInteraction=TRUE)$mat)[,"s0"]
  predict_white <- bind_rows(predict_white, data.frame(ID=data_white$ID[w], 
                                                  predicted=whitePredictionsCross))
}
data_white <- merge(data_white, predict_white, by="ID")

sampleID <- sample(1:4, size=nrow(data_black), replace=TRUE)
predict_black <- data.frame(ID=c(), predicted=c())
for (i in 1:4) {
  w <- which(sampleID==i)
  modelBlackCross <- 
    trainRidgeLASSO(effect="main", type="Ridge", transformation="Inv", 
                    new.data=data_black[-w,])
  blackPredictionsCross <- 
    1/predict(modelBlackCross, newx=makeMatrix(data_black[w,], includeInteraction=FALSE)$mat)[,"s0"]
  predict_black <- bind_rows(predict_black, data.frame(ID=data_black$ID[w], 
                                                  predicted=blackPredictionsCross))
}
data_black <- merge(data_black, predict_black, by="ID")


# the different prediction groups are OO, NN, ON and NO
data_white <- data_white %>% 
  mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<25, "NN", 
                                  ifelse(ObesityClass=="Normal weight"&predicted>30, "NO", 
                                         ifelse(ObesityClass=="Obese"&predicted>30, "OO", 
                                                ifelse(ObesityClass=="Obese"&predicted<25, "ON", NA)))))
data_black <- data_black %>% 
  mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<25, "NN", 
                                  ifelse(ObesityClass=="Normal weight"&predicted>30, "NO", 
                                         ifelse(ObesityClass=="Obese"&predicted>30, "OO", 
                                                ifelse(ObesityClass=="Obese"&predicted<25, "ON", NA)))))

# The difference in circulating metabolite levels was estimated with a Tukey-corrected ANOVA
whitePlotsDifference <- plotANOVA(data_white, ethnicity="white")
blackPlotsDifference <- plotANOVA(data_black, ethnicity="black")
ggsave("levelsNOwhite.pdf", whitePlotsDifference$normal)
ggsave("levelsONwhite.pdf", whitePlotsDifference$obese)
ggsave("levelsNOblack.pdf", blackPlotsDifference$normal)
ggsave("levelsONblack.pdf", blackPlotsDifference$obese)


## =============================================================================
# Output results of reclassifying overweight

#source("reclassify.R")


# training and validation set data were predicted with models trained on the full training data set

data_white_train$predicted <- 
  1/predict(object=modelWhiteTrain, 
            newx=makeMatrix(data_white_train, includeInteraction=TRUE)$mat,
            type="response")[,"s0"]
data_black_train$predicted <- 
  1/predict(object=modelBlackTrain, 
            newx=makeMatrix(data_black_train, includeInteraction=FALSE)$mat,
            type="response")[,"s0"]

data_white_val$predicted <- 
  1/predict(object=modelWhiteTrain, 
            newx=makeMatrix(data_white_val, includeInteraction=TRUE)$mat,
            type="response")[,"s0"]
data_black_val$predicted <- 
  1/predict(object=modelBlackTrain, 
            newx=makeMatrix(data_black_val, includeInteraction=FALSE)$mat,
            type="response")[,"s0"]


# data of interest for the two analyses
trainPredictions <- bind_rows(data_white_train, data_black_train)
rocTrain <- roc(trainPredictions$ObesityClass, trainPredictions$predicted, 
                levels=c("Normal weight", "Obese"))
cutoff <- coords(rocTrain, x="best")[1,"threshold"]

valPredictions <- bind_rows(data_white_val, data_black_val)
valPredictions <- valPredictions %>% 
  mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
  mutate(Reclassified=c("Normal", "Obese")[1+(predicted>cutoff)])

plotSelectionROC <- ggplot(valPredictions, aes(x=predicted, y=exp(BMI))) + 
                      geom_point(aes(col=selectROC, pch=Race)) +
                      scale_color_manual(values=c("#FF0000", "#00CCCC", "lightgray")) +
                      labs(title="BMI predictions in validation set", 
                           subtitle="selection for ROC analysis",
                           x="predicted BMI", y="observed BMI")
ggsave("plotSelectionROC.pdf", plotSelectionROC)

plotSelectionReclassify <- ggplot(valPredictions, aes(x=predicted, y=exp(BMI))) + 
                             geom_point(aes(col=Reclassified, pch=Race)) +
                             scale_color_manual(values=c("#00CCCC", "#FF0000")) +
                             geom_hline(yintercept=25, lty="dashed") + 
                             geom_hline(yintercept=30, lty="dashed") +  
                             labs(title="BMI predictions in validation set", 
                                  subtitle="selection for reclassification analysis",
                                  x="predicted BMI", y="observed BMI")
ggsave("plotSelectionReclassify.pdf", plotSelectionReclassify)

# ROC curve analysis of observed normal weight versus obese
rocAllBMIclass <- roc(formula=ObesityClass~predicted, data=valPredictions, levels=c("Normal weight", "Obese"))
aucAll <- ci.auc(rocAllBMIclass)
statsAll <- ci.coords(rocAllBMIclass, x=cutoff, input="threshold", 
                      ret=c("sensitivity", "specificity"))

rocWhiteBMIclass <- roc(formula=ObesityClass~predicted, data=data_white_val, levels=c("Normal weight", "Obese"))
aucWhite <- ci.auc(rocWhiteBMIclass)
statsWhite <- ci.coords(rocWhiteBMIclass, x=cutoff, input="threshold", 
                        ret=c("sensitivity", "specificity"))

rocBlackBMIclass <- roc(formula=ObesityClass~predicted, data=data_black_val, levels=c("Normal weight", "Obese"))
aucBlack <- ci.auc(rocBlackBMIclass)
statsBlack <- ci.coords(rocBlackBMIclass, x=cutoff, input="threshold", 
                        ret=c("sensitivity", "specificity"))

ROCstats <- matrix(ncol=4, nrow=3)
colnames(ROCstats) <- c("AUC", "threshold", "sensitivity", "specificity")
rownames(ROCstats) <- c("White", "Black", "All")
ROCstats[c("White", "Black", "All"), "threshold"] <- rep(round(cutoff, digits=1), times=3)
ROCstats["White", "sensitivity"] <- sprintf("%.2f[%.2f-%.2f]", statsWhite$sensitivity[1,2], 
                                            statsWhite$sensitivity[1,1], statsWhite$sensitivity[1,3])
ROCstats["Black", "sensitivity"] <- sprintf("%.2f[%.2f-%.2f]", statsBlack$sensitivity[1,2], 
                                            statsBlack$sensitivity[1,1], statsBlack$sensitivity[1,3])
ROCstats["All", "sensitivity"] <- sprintf("%.2f[%.2f-%.2f]", statsAll$sensitivity[1,2], 
                                            statsAll$sensitivity[1,1], statsAll$sensitivity[1,3])
ROCstats["White", "specificity"] <- sprintf("%.2f[%.2f-%.2f]", statsWhite$specificity[1,2], 
                                            statsWhite$specificity[1,1], statsWhite$specificity[1,3])
ROCstats["Black", "specificity"] <- sprintf("%.2f[%.2f-%.2f]", statsBlack$specificity[1,2], 
                                            statsBlack$specificity[1,1], statsBlack$specificity[1,3])
ROCstats["All", "specificity"] <- sprintf("%.2f[%.2f-%.2f]", statsAll$specificity[1,2], 
                                            statsAll$specificity[1,1], statsAll$specificity[1,3])
ROCstats["White", "AUC"] <- sprintf("%.2f[%.2f-%.2f]", aucWhite[2], aucWhite[1], aucWhite[3])
ROCstats["Black", "AUC"] <- sprintf("%.2f[%.2f-%.2f]", aucBlack[2], aucBlack[1], aucBlack[3])
ROCstats["All", "AUC"] <- sprintf("%.2f[%.2f-%.2f]", aucAll[2], aucAll[1], aucAll[3])

convertToTexTable(ROCstats, "ROCstats.tex", rows.named=TRUE, minipage=TRUE,
                  caption="ROC performance of metabolic BMI predicting observed BMI classes normal weight versus obese.", 
                  reflabel="ROCstats")

trueRatesAll <- data.frame(stratum="All", TPR=rocAllBMIclass$sensitivities, 
                           TNR=1-rocAllBMIclass$specificities)
trueRatesWhite <- data.frame(stratum="White", TPR=rocWhiteBMIclass$sensitivities, 
                             TNR=1-rocWhiteBMIclass$specificities)
trueRatesBlack <- data.frame(stratum="Black", TPR=rocBlackBMIclass$sensitivities, 
                             TNR=1-rocBlackBMIclass$specificities)
trueRates <- bind_rows(trueRatesAll, trueRatesWhite, trueRatesBlack)

ROCcurves <- ggplot(trueRates, aes(x=TNR, y=TPR, by=stratum)) + 
               geom_path(aes(color=stratum)) +
               scale_color_manual(values=c("black", "blue", "red")) +
               geom_segment(aes(x=0, y=0, xend=1, yend=1), lty="dashed") +
               labs(title="Detection rates observed normal weight versus obese")
ggsave("ROCcurves.pdf", ROCcurves)

# table of reclassified normal weight and obese among the different observed BMI classes
classCounts1 <- valPredictions %>% group_by(ObesityClass, Reclassified) %>% summarise(n=n())
classCounts1 <- pivot_wider(classCounts1, names_from="Reclassified", values_from="n")
classCounts1 <- classCounts1 %>% 
  mutate(Race="All") %>%
  mutate(ObesityClass=as.character(ObesityClass)) %>%
  mutate(`pred.: Normal`=Normal) %>%
  mutate(`pred.: Obese`=Obese) %>%
  mutate("fraction pred. obese"=calculateFraction(Normal, Obese))
classCounts1$Normal <- NULL
classCounts1$Obese <- NULL

classCounts2 <- valPredictions %>% group_by(Race, ObesityClass, Reclassified) %>% summarise(n=n())
classCounts2 <- pivot_wider(classCounts2, names_from="Reclassified", values_from="n")
classCounts2 <- classCounts2 %>% 
  mutate(Race=as.character(Race)) %>%
  mutate(ObesityClass=as.character(ObesityClass)) %>%
  mutate(`pred.: Normal`=Normal) %>%
  mutate(`pred.: Obese`=Obese) %>%
  mutate("fraction pred. obese"=calculateFraction(Normal, Obese))
classCounts2$Normal <- NULL
classCounts2$Obese <- NULL

classCounts <- bind_rows(classCounts2, classCounts1)

# test for significant difference in reclassified patients between two ethnicities
relNorm <- rep(classCounts1$`pred.: Normal`/(classCounts1$`pred.: Normal`+classCounts1$`pred.: Obese`), times=2)
relObese <- rep(classCounts1$`pred.: Obese`/(classCounts1$`pred.: Normal`+classCounts1$`pred.: Obese`), times=2)

indepNorm <- relNorm*(classCounts2$`pred.: Normal`+classCounts2$`pred.: Obese`)
indepObese <- relObese*(classCounts2$`pred.: Normal`+classCounts2$`pred.: Obese`)

X2 <- sum((classCounts2$`pred.: Normal`-indepNorm)**2/indepNorm) + sum((classCounts2$`pred.: Obese`-indepObese)**2/indepObese)
p_indep <- 1-pchisq(X2, df=3)

convertToTexTable(classCounts, "classCounts.tex", minipage=TRUE,
                  caption=sprintf("Number of patients reclassified as normal weight and obese, the fraction reclassified as obese. The fractions classified as obese are independent of ethnicity (p=%.2f)", p_indep),
                  reflabel="classCounts")


## =============================================================================
# Output results of discussion

#source("discussion.R")

data_outlier_white <- subset(data_outlier, subset= Race=="White")
data_outlier_white$predicted <- 1/predict(object=modelWhiteFull, 
                                          newx=makeMatrix(data_outlier_white, includeInteraction=TRUE)$mat,
                                          type="response")[,"s0"]

data_outlier_black <- subset(data_outlier, subset= Race=="Black")
data_outlier_black$predicted <- 1/predict(object=modelBlackFull, 
                                          newx=makeMatrix(data_outlier_black, includeInteraction=FALSE)$mat,
                                          type="response")[,"s0"]

data_outlier <- bind_rows(data_outlier_white, data_outlier_black)

outlierPrediction <- 
  ggplot(data=data_outlier, aes(x=predicted, y=exp(BMI))) + 
    geom_point(aes(col=Race)) + 
    scale_color_manual(values=c("red", "blue")) +
    geom_abline(slope=1, intercept=0, lty="dashed") +
    labs(title="BMI prediction of removed outliers", x="predicted BMI", y="observed BMI")

ggsave(filename="outlierPrediction.pdf", outlierPrediction)


