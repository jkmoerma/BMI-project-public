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
library(inspectdf)
library(mice)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(corrplot)
library(qgraph)
library(DMwR)
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

df$`met_093/met_012` <- NULL
df$`met_047/met_012` <- NULL

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

# tabulate number of patients for every ethnicity and BMI class
biometrics <- tableBiometrics(df)
convertToTexTable(biometrics, "explore_biometrics.tex", rows.named=TRUE,
                  caption="Sample patient count for every ethnicity and obesity class.", 
                  reflabel="sample_biometrics")

# tabulate number of smoking patients 
smokingCounts <- tableSmokingStatus(df)
convertToTexTable(smokingCounts, "explore_smoking.tex", rows.named=TRUE,
                  caption="Number of smoking and non-smoking patients for every ethnicity.",
                  reflabel="explore_smoking")



# investigate erase outliers
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
w <- which(df_complete$ID==df$ID[which(is.na(df$met_010))])
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

# investigate metabolites associated with smoking with a bayesian network
if (FALSE) {
  parentAge <- data.frame(from=c("BMI","Race","Smoking",metabolites), 
                          to=rep("Age", times=length(metabolites)+3))
  parentRace <- data.frame(from=c("BMI","Age","Smoking",metabolites), 
                           to=rep("Race", times=length(metabolites)+3))
  parentSmoking <- data.frame(from=c("BMI","Age","Race", metabolites), 
                              to=rep("Smoking", times=length(metabolites)+3))
  blacklist <- bind_rows(parentAge, parentRace, parentSmoking)
  networkStructure <- hc(data_model[,c("BMI","Age","Race", "Smoking", metabolites)], score="bic-cg", blacklist=blacklist)
  
  # explore smoking
  smokeMets <- networkStructure$nodes$Smoking$children
  for (met in smokeMets) {
    boxplot(eval(parse(text=paste(met, "~ Smoking"))), data=data_model)
  }
  
  # explore ethnicity
  raceMets <- networkStructure$nodes$Race$children
  for (met in raceMets) {
    boxplot(eval(parse(text=paste(met, "~ Race"))), data=data_model, 
            subset = ObesityClass=="Normal weight")
  }
  
  # explore age
  ageMets <- networkStructure$nodes$Age$children
  for (met in ageMets) {
    trend <- lm(eval(parse(text=paste(met, "~ Age*Race"))), data=data_model, 
                subset = ObesityClass=="Normal weight")
    plot(eval(parse(text=paste(met, "~ Age"))), data=data_model, 
         subset = ObesityClass=="Normal weight")
    abline(a=trend$coefficients["(Intercept)"], b=trend$coefficients["Age"], col="red")
  }
  
}

boxplotSmoking062 <- ggplot(data=data_model, 
                         aes(x=as.factor(c("non-smoking", "smoking")[1+Smoking]), 
                             y=met_062)) + 
  geom_boxplot(aes(col=as.factor(Smoking))) +
  theme(axis.text.x = element_text(angle = 0),
        legend.title = element_blank(), legend.position="none") +
  xlab("") + ylab("log(met_062)") + labs(title="Circulating levels of met_062")
ggsave("explore_boxplotSmoking062.pdf", boxplotSmoking062)

boxplotSmokingBMIwhite <- ggplot(data=subset(data_model, subset = Race=="White"), 
                            aes(x=as.factor(c("non-smoking", "smoking")[1+Smoking]), 
                                y=exp(BMI))) + 
  geom_boxplot(aes(col=as.factor(Smoking))) +
  theme(axis.text.x = element_text(angle = 0),
        legend.title = element_blank(), legend.position="none") +
  xlab("") + ylab("BMI") + labs(title="BMI smoking/non-smoking (White)")
ggsave("explore_boxplotSmokingBMIwhite.pdf", boxplotSmokingBMIwhite)

boxplotSmokingBMIblack <- ggplot(data=subset(data_model, subset = Race=="Black"), 
                                 aes(x=as.factor(c("non-smoking", "smoking")[1+Smoking]), 
                                     y=exp(BMI))) + 
  geom_boxplot(aes(col=as.factor(Smoking))) +
  theme(axis.text.x = element_text(angle = 0),
        legend.title = element_blank(), legend.position="none") +
  xlab("") + ylab("BMI") + labs(title="BMI smoking/non-smoking (Black)")
ggsave("explore_boxplotSmokingBMIblack.pdf", boxplotSmokingBMIblack)

raceMets <- c("met_002", "met_034", "met_068", "met_149")
summarizeRaceDiffs <- c()
for (met in raceMets) {
  formulastring <- paste0("scale(", met, ")~Race")
  differences <- TukeyHSD(aov(eval(parse(text=formulastring)), data=data_model, 
                              subset = ObesityClass=="Normal weight"))$Race
  differences <- differences[c("Black-White", "South Asian-White", 
                               "East Asian-White", "Mixed-White"),]
  metEntry <- sprintf("%.2f [%.2f %.2f]", differences[,"diff"],
                      differences[,"lwr"], differences[,"upr"])
  names(metEntry) <- rownames(differences)
  summarizeRaceDiffs <- bind_rows(summarizeRaceDiffs, c(met=met, metEntry))
}
summarizeRaceDiffsMatrix <- as.matrix(summarizeRaceDiffs[,-1])
rownames(summarizeRaceDiffsMatrix) <- summarizeRaceDiffs$met
convertToTexTable(summarizeRaceDiffsMatrix, "explore_ethnicity.tex", rows.named=TRUE,
                  caption="Difference in standardized log concentration of metabolites between ethnicities among normal weight patients.",
                  reflabel="explore_ethnicity")

summarizeAge <- df %>% group_by(Race, ObesityClass) %>% 
                  summarise("[IQR]"=sprintf("%.1f [%.1f %.1f]", median(Age), 
                                                       quantile(Age, probs=0.25), 
                                                       quantile(Age, probs=0.75)))
summarizeAge <- pivot_wider(summarizeAge, names_from="Race", values_from="[IQR]",
                            names_glue = "{Race} {.value}")
summarizeAgeMatrix <- as.matrix(summarizeAge[,-1])
rownames(summarizeAgeMatrix) <- summarizeAge$ObesityClass
convertToTexTable(summarizeAgeMatrix, "explore_age.tex", rows.named=TRUE,
                  caption="Median + IQR of patient ages. Grouping by ethnicity and BMI class.",
                  reflabel="explore_age")

# explore correlations
cor_matrix_white <- cor(subset(data_model, subset = Race=="White", 
                         select = metabolites))
qgraph(cor_matrix_white, layout = "spring", threshold = 0.0, 
       labels = gsub(pattern="met_", replacement="", fixed=TRUE, 
                     x=colnames(cor_matrix_white)),
       label.scale.equal=TRUE, label.cex=2.5,
       filetype="pdf", filename="explore_correlationsWhite.pdf",
       mar = c(0.5, 1, 0.5, 0.5))

cor_matrix_black <- cor(subset(data_model, subset = Race=="Black", 
                               select = metabolites))
qgraph(cor_matrix_black, layout = "spring", threshold = 0.0, 
       labels = gsub(pattern="met_", replacement="", fixed=TRUE, 
                     x=colnames(cor_matrix_black)),
       label.scale.equal=TRUE, label.cex=2.5,
       filetype="pdf", filename="explore_correlationsBlack.pdf",
       mar = c(0.5, 0.5, 0.5, 1))

stratified_correlations <- data_model %>% group_by(Race) %>% inspect_cor() %>% 
                             filter(abs(corr)>0.85)
strongest <- unique(c(stratified_correlations$col_1, stratified_correlations$col_2))
plot_correlations <- c()
for (met in strongest) {
  met_correlations1 <- stratified_correlations %>% filter(col_1==met)
  met_correlations2 <- stratified_correlations %>% filter(col_2==met)
  
  names(met_correlations1)[c(2,3)] <- names(met_correlations1)[c(3,2)]
  
  plot_correlations <- bind_rows(plot_correlations, met_correlations2)
  plot_correlations <- bind_rows(plot_correlations, met_correlations1)
}
plot_correlations$posX <- 1:nrow(plot_correlations)
plot_correlations$labX <- paste(plot_correlations$col_1, 
                                plot_correlations$col_2, sep=" & ")
p_correlations <- ggplot(plot_correlations,aes(x=reorder(labX,posX), y=corr, 
                                            by=Race, col=col_2)) +
  geom_point(aes(pch=Race)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=-0.85, col="red") +
  geom_hline(yintercept=0.85, col="red") +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank(), 
        legend.position="none") +
  xlab("") + ylab("Pearson correlation") +
  coord_flip() +
  labs(title="Selection of correlations > 0.85")
ggsave("explore_correlations.pdf", p_correlations, width=7, height=9)

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

if (!file.exists("regression_trainWhite.tex")) {
  trainWhite <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations, ethnicity="White", 
                                   balancing, data_white_train)
  trainBlack <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations, ethnicity="Black", 
                                   new.data=data_black_train)
  convertToTexTable(trainWhite, "regression_trainWhite.tex", 
                    caption="Cross-validated prediction evaluation of the models on white ethnicity training data.", 
                    reflabel="regression_trainWhite")
  convertToTexTable(trainBlack, "regression_trainBlack.tex", 
                    caption="Cross-validated prediction evaluation of the models on black ethnicity training data.", 
                    reflabel="regression_trainBlack")
}


modelWhiteFull <- 
  trainRidgeLASSO(effect="interaction", type="Ridge", transformation="Inv", new.data=oversample(data_white))
modelBlackFull <- 
  trainRidgeLASSO(effect="main", type="Ridge", transformation="Inv", new.data=data_black)



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

if (!file.exists("metabolic_CIcoeffsMain.pdf") | 
    !file.exists("metabolic_CIcoeffsIntMain.pdf") | 
    !file.exists("metabolic_CIcoeffsIntInt.pdf")) {
  
  # generate coefficient estimates for main effects modeling
  bootstrapCoeffsWhiteMain <- scaledEffects(data=data_white, effect="main", 
                                        type="Ridge", transformation="Inv", 
                                        balancing="Balanced", boot.n=200)
  bootstrapCoeffsWhiteMain <- pivot_longer(bootstrapCoeffsWhiteMain, 
                                       cols=colnames(bootstrapCoeffsWhiteMain),
                                       names_to="met", values_to="coeff")
  bootstrapCoeffsBlackMain <- scaledEffects(data=data_black, effect="main", 
                                        type="Ridge", transformation="Inv", 
                                        balancing="", boot.n=200)
  bootstrapCoeffsBlackMain <- pivot_longer(bootstrapCoeffsBlackMain, 
                                       cols=colnames(bootstrapCoeffsBlackMain),
                                       names_to="met", values_to="coeff")
  
  # store plot of beta-coefficients main effect modeling
  p_CIcoeffsMain <- plotScaledEffects(bootstrapCoeffsWhiteMain, 
                                      bootstrapCoeffsBlackMain, effect="main")
  ggsave("metabolic_CIcoeffsMain.pdf", p_CIcoeffsMain)
  
  
  
  # generate coefficient estimates for interaction effects modeling
  bootstrapCoeffsWhiteInt <- scaledEffects(data=data_white, effect="interaction", 
                                            type="Ridge", transformation="Inv", 
                                            balancing="Balanced", boot.n=200)
  bootstrapCoeffsWhiteInt <- pivot_longer(bootstrapCoeffsWhiteInt, 
                                           cols=colnames(bootstrapCoeffsWhiteInt),
                                           names_to="met", values_to="coeff")
  bootstrapCoeffsBlackInt <- scaledEffects(data=data_black, effect="interaction", 
                                            type="Ridge", transformation="Inv", 
                                            balancing="", boot.n=200)
  bootstrapCoeffsBlackInt <- pivot_longer(bootstrapCoeffsBlackInt, 
                                           cols=colnames(bootstrapCoeffsBlackInt),
                                           names_to="met", values_to="coeff")
  
  # store plots of beta-coefficients interaction effects modeling
  p_CIcoeffsInt <- plotScaledEffects(bootstrapCoeffsWhiteInt, 
                                     bootstrapCoeffsBlackInt, effect="interaction")
  ggsave("metabolic_CIcoeffsIntMain.pdf", p_CIcoeffsInt$main)
  ggsave("metabolic_CIcoeffsIntInt.pdf", p_CIcoeffsInt$interaction)
  
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

# make plot for illustration of these prediction groups
countPredictionGroup <- bind_rows(table(data_white$predictionGroup), table(data_black$predictionGroup))
countPredictionGroup <- as.matrix(countPredictionGroup)
rownames(countPredictionGroup) <- c("White", "Black")
convertToTexTable(countPredictionGroup, "metabolic_predictionGroup.tex", rows.named=TRUE,
                  caption="Number of patients in every prediction group for the ANOVA on metabolic outliers.", 
                  reflabel="metabolic_predictionGroup")
illustrateGroups <- ggplot(data=bind_rows(data_white, data_black),
                           aes(x=predicted, y=exp(BMI))) +
                      geom_point(aes(col=predictionGroup, pch=Race)) +
                      labs(title="Metabolic normals (NN, OO) and outliers (NO, ON)", 
                           x="predicted BMI", y="observed BMI")
ggsave("metabolic_predictionGroupDemonstration.pdf", illustrateGroups)


# The difference in circulating metabolite levels was estimated with a Tukey-corrected ANOVA
whitePlotsDifference <- plotANOVA(data_white, ethnicity="white")
blackPlotsDifference <- plotANOVA(data_black, ethnicity="black")
ggsave("metabolic_levelsNOwhite.pdf", whitePlotsDifference$normal)
ggsave("metabolic_levelsONwhite.pdf", whitePlotsDifference$obese)
ggsave("metabolic_levelsNOblack.pdf", blackPlotsDifference$normal)
ggsave("metabolic_levelsONblack.pdf", blackPlotsDifference$obese)


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

data_white_smoking$predicted <- 
  1/predict(object=modelWhiteTrain, 
            newx=makeMatrix(data_white_smoking, includeInteraction=TRUE)$mat,
            type="response")[,"s0"]
data_black_smoking$predicted <- 
  1/predict(object=modelBlackTrain, 
            newx=makeMatrix(data_black_smoking, includeInteraction=FALSE)$mat,
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
ggsave("reclassify_plotSelectionROC.pdf", plotSelectionROC)

plotSelectionReclassify <- ggplot(valPredictions, aes(x=predicted, y=exp(BMI))) + 
                             geom_point(aes(col=Reclassified, pch=Race)) +
                             scale_color_manual(values=c("#00CCCC", "#FF0000")) +
                             geom_hline(yintercept=25, lty="dashed") + 
                             geom_hline(yintercept=30, lty="dashed") +  
                             labs(title="BMI predictions in validation set", 
                                  subtitle="selection for reclassification analysis",
                                  x="predicted BMI", y="observed BMI")
ggsave("reclassify_plotSelectionReclassify.pdf", plotSelectionReclassify)

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

convertToTexTable(ROCstats, "reclassify_ROCstats.tex", rows.named=TRUE, minipage=TRUE,
                  caption="ROC performance of metabolic BMI predicting observed BMI classes normal weight versus obese.", 
                  reflabel="reclassify_ROCstats")

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
ggsave("reclassify_ROCcurves.pdf", ROCcurves)

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

convertToTexTable(classCounts, "reclassify_classCounts.tex", minipage=TRUE,
                  caption=sprintf("Number of patients reclassified as normal weight and obese, the fraction reclassified as obese. The fractions classified as obese are independent of ethnicity (p=%.2e)", p_indep),
                  reflabel="reclassify_classCounts")

# comparing non-smokers of the validation set with smokers
all_data <- bind_rows(data_white_val, data_white_smoking, 
                      data_black_val, data_black_smoking)

p_smokeWhite <- ggplot(data=bind_rows(data_white_val, data_white_smoking), 
                       aes(x=predicted, y=exp(BMI))) + 
                  geom_point(aes(col=as.factor(Smoking))) +
                  scale_color_manual(values=c("green", "brown")) +
                  geom_hline(yintercept=25, lty="dashed") + 
                  geom_hline(yintercept=30, lty="dashed") +
                  labs(title="Non-smokers vs. smokers (White)", 
                       xlab="predicted BMI", ylab="observed BMI")
p_smokeBlack <- ggplot(data=bind_rows(data_black_val, data_black_smoking), 
                       aes(x=predicted, y=exp(BMI))) + 
                  geom_point(aes(col=as.factor(Smoking))) +
                  scale_color_manual(values=c("green", "brown")) +
                  geom_hline(yintercept=25, lty="dashed") + 
                  geom_hline(yintercept=30, lty="dashed") +
                  labs(title="Non-smokers vs. smokers (Black)", 
                       xlab="predicted BMI", ylab="observed BMI")
ggsave("reclassify_compareSmokingWhite.pdf", p_smokeWhite)
ggsave("reclassify_compareSmokingBlack.pdf", p_smokeBlack)


countSample <- all_data %>% group_by(Race, ObesityClass, Smoking) %>% summarise(n=n())
countSample <- pivot_wider(countSample, values_from="n", names_from=Smoking)

confIntMedian <- function(wilcoxTest) {
  sprintf("%.2f[%.2f %.2f]", wilcoxTest$estimate, wilcoxTest$conf.int[1], 
          wilcoxTest$conf.int[2])
}

table_smoke <- all_data %>% group_by(Race, ObesityClass) %>% 
                 summarise("pred. BMI (non smoking-smoking)"= confIntMedian(wilcox.test(predicted~Smoking, alternative="two.sided", conf.int=TRUE)), 
                           "p-val."=wilcox.test(predicted~Smoking, alternative="two.sided", conf.int=TRUE)$p.value)
table_smoke$Race <- as.character(table_smoke$Race)
table_smoke$ObesityClass <- as.character(table_smoke$ObesityClass)
table_smoke$`p-val.` <- signif(table_smoke$`p-val.`, digits=2)
convertToTexTable(table_smoke, "reclassify_compareSmokingTable.tex",
                  caption="Difference in median predicted BMI + 95\\% CI (non smoking - smoking). The prediction model trained on the training data was used to predict validation set patients and smoking patients.",
                  reflabel="reclassify_compareSmokingTable")




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


