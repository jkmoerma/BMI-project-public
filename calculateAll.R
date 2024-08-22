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

# set default width adn height for saving ggplot objects
ggplotWidth <- 5
ggplotHeight <- 5

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
source("reclassify.R")
source("metabolicAnalysis.R")
source("outliers.R")
source("smoking.R")

## =============================================================================
# Output results of the data exploration

#source("dataExploration.R")

# tabulate number of patients for every ethnicity and BMI class
biometrics <- tableBiometrics(df)
convertToTexTable(biometrics, "explore_biometrics.tex", rows.named=TRUE,
                  caption="Sample patient count for every ethnicity and obesity class.", 
                  reflabel="explore_biometrics")

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
    boxplot(eval(parse(text=paste(met, "~ Smoking*Race"))), data=data_model)
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
  theme(legend.title = element_blank(), legend.position="none",
        title =element_text(size=16, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=17),
        axis.title.y = element_text(size=12,color='black',face='bold')) +
  labs(title="Circulating levels of met_062", x="", y="log(met_062)")
ggsave("explore_boxplotSmoking062.pdf", boxplotSmoking062, 
       width=ggplotWidth, height=ggplotHeight)

boxplotSmokingBMIwhite <- ggplot(data=subset(data_model, subset = Race=="White"), 
                            aes(x=as.factor(c("non-smoking", "smoking")[1+Smoking]), 
                                y=exp(BMI))) + 
  geom_boxplot(aes(col=as.factor(Smoking))) +
  theme(legend.title = element_blank(), legend.position="none",
        title =element_text(size=15, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=17),
        axis.title.y = element_text(size=12,color='black',face='bold')) +
  labs(title="BMI smoking/non-smoking (White)", x="", y="BMI")
ggsave("explore_boxplotSmokingBMIwhite.pdf", boxplotSmokingBMIwhite, 
       width=ggplotWidth, height=ggplotHeight)

boxplotSmokingBMIblack <- ggplot(data=subset(data_model, subset = Race=="Black"), 
                                 aes(x=as.factor(c("non-smoking", "smoking")[1+Smoking]), 
                                     y=exp(BMI))) + 
  geom_boxplot(aes(col=as.factor(Smoking))) +
  theme(legend.title = element_blank(), legend.position="none",
        title =element_text(size=15, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=17),
        axis.title.y = element_text(size=12,color='black',face='bold')) +
  labs(title="BMI smoking/non-smoking (Black)", x="", y="BMI")
ggsave("explore_boxplotSmokingBMIblack.pdf", boxplotSmokingBMIblack, 
       width=ggplotWidth, height=ggplotHeight)

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
                  caption="Difference in standardized log concentration of metabolites between ethnicities among normal weight patients, estimate + 95\\% CI.",
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
                         select = c("Age", metabolites)))
qgraph(cor_matrix_white, layout = "spring", threshold = 0.0, 
       labels = gsub(pattern="met_", replacement="", fixed=TRUE, 
                     x=colnames(cor_matrix_white)),
       label.scale.equal=TRUE, label.cex=2.5,
       filetype="pdf", filename="explore_correlationsWhite.pdf",
       mar = c(0.5, 1, 0.5, 0.5))

cor_matrix_black <- cor(subset(data_model, subset = Race=="Black", 
                               select = c("Age", metabolites)))
qgraph(cor_matrix_black, layout = "spring", threshold = 0.0, 
       labels = gsub(pattern="met_", replacement="", fixed=TRUE, 
                     x=colnames(cor_matrix_black)),
       label.scale.equal=TRUE, label.cex=2.5,
       filetype="pdf", filename="explore_correlationsBlack.pdf",
       mar = c(0.5, 0.5, 0.5, 1))


# filter out outliers of data set
filtered_data <- filterOutliers(data_model, outliers=distributionList$outliers, 
                                mets=c("met_011", "met_015", "met_017", "met_027", 
                                       "met_029", "met_031", "met_034", "met_038",
                                       "met_049", "met_059", "met_060", "met_064", 
                                       "met_066", "met_073", "met_084", "met_085",
                                       "met_088", "met_090", "met_132", "met_134"))
data_regular <- filtered_data$regulars
data_influential <- filtered_data$outliers
  

## =============================================================================
# Output results of the regression analysis

#source("linearRegression.R")

# white ethnicity data

data_white <- subset(data_regular, subset = Race=="White"&(!Smoking))
data_white_smoking <- subset(data_regular, subset = Race=="White"&Smoking)
set.seed(23/8/2024)
valSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
data_white_test <- data_white[valSelectA,]
data_white_train <- data_white[-valSelectA,]

data_white_train_balanced <- oversample(data_white_train)


# black ethnicity data

data_black <- subset(data_regular, subset = Race=="Black"&(!Smoking))
data_black_smoking <- subset(data_regular, subset = Race=="Black"&Smoking)
set.seed(23/8/2024)
valSelectB <- sample(x=1:nrow(data_black), size=round(0.2*nrow(data_black)))
data_black_test <- data_black[valSelectB,]
data_black_train <- data_black[-valSelectB,]


# all ethnicity data

data_all <- bind_rows(data_white, data_black)
data_all_test <- bind_rows(data_white_test, data_black_test)
data_all_train <- bind_rows(data_white_train, data_black_train)


# tabulate cross validated performance of all model recipes

formulations <- expand.grid(effects=c("main", "interaction"), 
                            types=c("OLS", "Ridge", "LASSO"), 
                            transformations=c("Log", "Inv"),
                            stringsAsFactors=FALSE)
balancing <- rep("Balanced", times=nrow(formulations))

if (!file.exists("regression_trainWhite.tex")) {
  trainWhite <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations,
                                   balancing, data_white_train)
  trainBlack <- tabulateValidation(formulations$effects, formulations$types, 
                                   formulations$transformations, 
                                   new.data=data_black_train)
  
  # choose one formulation for both ethnicities
  w <- which.max(trainWhite$AUC.cv+trainBlack$AUC.cv)
  formulationWhite <- list(effect=formulations$effects[w], 
                           type=formulations$types[w],
                           transformation=formulations$transformations[w])
  formulationBlack <- list(effect=formulations$effects[w], 
                           type=formulations$types[w],
                           transformation=formulations$transformations[w])
  
  # highlight chosen formulation in table
  trainWhite[w, 1] <- paste("\\rowcolor{green}", trainWhite[w, 1])
  trainBlack[w, 1] <- paste("\\rowcolor{green}", trainBlack[w, 1])
  
  convertToTexTable(trainWhite, "regression_trainWhite.tex", 
                    caption="Cross validated prediction evaluation of the models on white ethnicity training data.", 
                    reflabel="regression_trainWhite")
  convertToTexTable(trainBlack, "regression_trainBlack.tex", 
                    caption="Cross validated prediction evaluation of the models on black ethnicity training data.", 
                    reflabel="regression_trainBlack")
  
  
  saveRDS(formulationWhite, "regression_formulationWhite.rds")
  saveRDS(formulationBlack, "regression_formulationBlack.rds")
  
}

formulationWhite <- readRDS("regression_formulationWhite.rds")
formulationBlack <- readRDS("regression_formulationBlack.rds")

modelWhiteTrain <- 
  trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                  transformation=formulationWhite$transformation, 
                  new.data=oversample(data_white_train))
modelBlackTrain <- 
  trainRidgeLASSO(effect=formulationBlack$effect, type=formulationBlack$type, 
                  transformation=formulationBlack$transformation, 
                  new.data=data_black_train)
modelAllTrain <- 
  trainRidgeLASSO(effect=formulationBlack$effect, type=formulationBlack$type, 
                  transformation=formulationBlack$transformation, 
                  new.data=oversample(data_all_train))


# illlustrate the effect of oversampling to correct for the BMI class imbalance in white ethnicity patients

modelWhiteUnbalanced <- 
  trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                  transformation=formulationWhite$transformation, 
                  new.data=data_white_train)

whitePredictionsUnbalanced <- 
  exp(predict(modelWhiteUnbalanced, 
              newx=makeMatrix(data_white_train, includeInteraction=formulationWhite$effect=="interaction")$mat)[,"s0"])
whitePredictionsBalanced <- 
  exp(predict(modelWhiteTrain, 
              newx=makeMatrix(data_white_train, includeInteraction=formulationWhite$effect=="interaction")$mat)[,"s0"])

unbalancedPredictions <- data.frame("observed BMI"=exp(data_white_train$BMI),
                                    "predicted BMI"=whitePredictionsUnbalanced,
                                    "balanced"=rep("no", times=nrow(data_white_train)),
                                    "end segment"=whitePredictionsBalanced,
                                    check.names=FALSE)
balancedPredictions <- data.frame("observed BMI"=exp(data_white_train$BMI),
                                  "predicted BMI"=whitePredictionsBalanced,
                                  "balanced"=rep("yes", times=nrow(data_white_train)),
                                  "end segment"=whitePredictionsUnbalanced,
                                  check.names=FALSE)
p_balancing <- ggplot(data=bind_rows(unbalancedPredictions, balancedPredictions),
                      aes(x=`predicted BMI`, y=`observed BMI`)) +
  geom_point(aes(col=balanced)) + 
  geom_abline(intercept=0, slope=1, lty="dashed") + 
  geom_segment(data=unbalancedPredictions,
               arrow=arrow(angle=30, length=unit(0.1, "inches"),
                           ends="last", type="open"),
               aes(x=`predicted BMI`, y=`observed BMI`, 
                   xend=`end segment`, yend=`observed BMI`)) + 
  labs(title="BMI predictions: unbalanced and balanced with SMOTE", 
       subtitle=paste0("effects: ", formulationWhite$effect, ", type: ", 
                       formulationWhite$type, ", BMI transformation: ", 
                       formulationWhite$transformation))
ggsave("regression_balancing.pdf", p_balancing, width=ggplotWidth, height=ggplotHeight)

 
# visualize linearity assumption for the different BMI transformations
modelWhiteLog <- 
  trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                  transformation="Log", 
                  new.data=oversample(data_white_train))
modelWhiteInv <- 
  trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                  transformation="Inv", 
                  new.data=oversample(data_white_train))
modelBlackLog <- 
  trainRidgeLASSO(effect=formulationBlack$effect, type=formulationBlack$type, 
                  transformation="Log", 
                  new.data=data_black_train)
modelBlackInv <- 
  trainRidgeLASSO(effect=formulationBlack$effect, type=formulationBlack$type, 
                  transformation="Inv", 
                  new.data=data_black_train)

whitePredictionsLog <- 
  predict(modelWhiteTrain, 
          newx=makeMatrix(data_white_train, 
                          includeInteraction = formulationWhite$effect=="interaction")$mat)[,"s0"]
whitePredictionsInv <- 
  predict(modelWhiteInv, 
          newx=makeMatrix(data_white_train, 
                          includeInteraction = formulationBlack$effect=="interaction")$mat)[,"s0"]
p_linearityWhiteLog <- ggplot(data.frame("predicted log(BMI)"=whitePredictionsLog,
                                         "residual"=data_white_train$BMI-whitePredictionsLog,
                                         check.names=FALSE), 
                              aes(x=`predicted log(BMI)`, y=residual)) + 
  geom_point() + geom_smooth() +
  geom_hline(yintercept=0, lty="dashed") +
  labs(title="Linearity check: transformation log(BMI)", 
       subtitle=paste0("Ethnicity: white, effects: ", formulationWhite$effect, 
                       ", type: ", formulationWhite$type))
p_linearityWhiteInv <- ggplot(data.frame("predicted 1/BMI"=whitePredictionsInv,
                                         "residual"=exp(-data_white_train$BMI)-whitePredictionsInv,
                                         check.names=FALSE), 
                              aes(x=`predicted 1/BMI`, y=residual)) + 
  geom_point() + geom_smooth() + 
  geom_hline(yintercept=0, lty="dashed") +
  labs(title="Linearity check: transformation 1/BMI", 
       subtitle=paste0("Ethnicity: white, effects: ", formulationWhite$effect, 
                       ", type: ", formulationWhite$type))
ggsave("regression_linearityWhiteLog.pdf", p_linearityWhiteLog, width=ggplotWidth, height=ggplotHeight)
ggsave("regression_linearityWhiteInv.pdf", p_linearityWhiteInv, width=ggplotWidth, height=ggplotHeight)

blackPredictionsLog <- 
  predict(modelBlackTrain, 
          newx=makeMatrix(data_black_train, includeInteraction=formulationBlack$effect=="interaction")$mat)[,"s0"]
blackPredictionsInv <- 
  predict(modelBlackInv, 
          newx=makeMatrix(data_black_train, includeInteraction=formulationBlack$effect=="interaction")$mat)[,"s0"]
p_linearityBlackLog <- ggplot(data.frame("predicted log(BMI)"=blackPredictionsLog,
                                         "residual"=data_black_train$BMI-blackPredictionsLog,
                                         check.names=FALSE), 
                              aes(x=`predicted log(BMI)`, y=residual)) + 
  geom_point() + geom_smooth() + geom_hline(yintercept=0, lty="dashed") +
  labs(title="Linearity check: transformation log(BMI)", 
       subtitle=paste0("Ethnicity: black, effects: ", formulationBlack$effect, 
                       ", type: ", formulationBlack$type))
p_linearityBlackInv <- ggplot(data.frame("predicted 1/BMI"=blackPredictionsInv,
                                         "residual"=exp(-data_black_train$BMI)-blackPredictionsInv,
                                         check.names=FALSE), 
                              aes(x=`predicted 1/BMI`, y=residual)) + 
  geom_point() + geom_smooth() + geom_hline(yintercept=0, lty="dashed") +
  labs(title="Linearity check: transformation 1/BMI", 
       subtitle=paste0("Ethnicity: black, effects: ", formulationBlack$effect, 
                       ", type: ", formulationBlack$type))
ggsave("regression_linearityBlackLog.pdf", p_linearityBlackLog, width=ggplotWidth, height=ggplotHeight)
ggsave("regression_linearityBlackInv.pdf", p_linearityBlackInv, width=ggplotWidth, height=ggplotHeight)


# demonstrate lambda parameter tuning for ridge regression
lambdaRidgeWhite <- 
  tuneLambda(data_matrix=makeMatrix(data_white_train_balanced, 
                                    includeInteraction=formulationWhite$effect=="interaction")$mat, 
             alpha=0, new.y=data_white_train_balanced$BMI, 
             transformation=formulationWhite$transformation, 
             effect=formulationWhite$effect, returnAll=TRUE)
bestLambdaRidgeWhite <- lambdaRidgeWhite$lambda[which.max(lambdaRidgeWhite$AUC)]
p_tuneLambdaRidgeWhite <- ggplot(lambdaRidgeWhite, aes(x=lambda, y=AUC)) +
  geom_point() + scale_x_log10() +
  geom_vline(xintercept=bestLambdaRidgeWhite, col="blue") +
  labs(title="Demonstration of parameter tuning: Ridge regression",
       subtitle=paste0("Ethnicity: white, effects: ", formulationWhite$effect, 
                       ", BMI transformation: ", formulationWhite$transformation))
ggsave("regression_tuneLambdaRidgeWhite.pdf", p_tuneLambdaRidgeWhite, width=ggplotWidth, height=ggplotHeight)

lambdaRidgeBlack <- 
  tuneLambda(data_matrix=makeMatrix(data_black_train, 
                                    includeInteraction=formulationBlack$effect=="interaction")$mat, 
             alpha=0, new.y=data_black_train$BMI, 
             transformation=formulationBlack$transformation, 
             effect=formulationBlack$effect, returnAll=TRUE)
bestLambdaRidgeBlack <- lambdaRidgeBlack$lambda[which.max(lambdaRidgeBlack$AUC)]
p_tuneLambdaRidgeBlack <- ggplot(lambdaRidgeBlack, aes(x=lambda, y=AUC)) +
  geom_point() + scale_x_log10() +
  geom_vline(xintercept=bestLambdaRidgeBlack, col="blue") +
  labs(title="Demonstration of parameter tuning: Ridge regression",
       subtitle=paste0("Ethnicity: black, effects: ", formulationBlack$effect, 
                       ", BMI transformation: ", formulationBlack$transformation))
ggsave("regression_tuneLambdaRidgeBlack.pdf", p_tuneLambdaRidgeBlack, width=ggplotWidth, height=ggplotHeight)


# demonstrate lambda parameter tuning for LASSO regression
lambdaLASSOWhite <- 
  tuneLambda(data_matrix=makeMatrix(data_white_train_balanced, 
                                    includeInteraction=formulationWhite$effect=="interaction")$mat, 
             alpha=1, new.y=data_white_train_balanced$BMI, 
             transformation=formulationWhite$transformation, 
             effect=formulationWhite$effect, returnAll=TRUE)
bestLambdaLASSOWhite <- lambdaLASSOWhite$lambda[which.max(lambdaLASSOWhite$AUC)]
p_tuneLambdaLASSOWhite <- ggplot(lambdaLASSOWhite, aes(x=lambda, y=AUC)) +
  geom_point() + scale_x_log10() +
  geom_vline(xintercept=bestLambdaLASSOWhite, col="blue") +
  labs(title="Demonstration of parameter tuning: LASSO regression",
       subtitle=paste0("Ethnicity: white, effects: ", formulationWhite$effect, 
                       ", BMI transformation: ", formulationWhite$transformation))
ggsave("regression_tuneLambdaLASSOWhite.pdf", p_tuneLambdaLASSOWhite, width=ggplotWidth, height=ggplotHeight)

lambdaLASSOBlack <- 
  tuneLambda(data_matrix=makeMatrix(data_black_train, 
                                    includeInteraction=formulationBlack$effect=="interaction")$mat, 
             alpha=1, new.y=data_black_train$BMI, 
             transformation=formulationBlack$transformation, 
             effect=formulationBlack$effect, returnAll=TRUE)
bestLambdaLASSOBlack <- lambdaLASSOBlack$lambda[which.max(lambdaLASSOBlack$AUC)]
p_tuneLambdaLASSOBlack <- ggplot(lambdaLASSOBlack, aes(x=lambda, y=AUC)) +
  geom_point() + scale_x_log10() +
  geom_vline(xintercept=bestLambdaLASSOBlack, col="blue") +
  labs(title="Demonstration of parameter tuning: LASSO regression",
       subtitle=paste0("Ethnicity: black, effects: ", formulationBlack$effect, 
                       ", BMI transformation: ", formulationBlack$transformation))
ggsave("regression_tuneLambdaLASSOBlack.pdf", p_tuneLambdaLASSOBlack, width=ggplotWidth, height=ggplotHeight)



# make posters of model diagnostics of chosen models

modelDiagnostics(effects=formulationWhite$effect, types=formulationWhite$type, 
                 transformations=formulationWhite$transformation, 
                 balancing="Balanced", ethnicity="White", model=modelWhiteTrain, 
                 new.data=data_white_train)
modelDiagnostics(effects=formulationBlack$effect, types=formulationBlack$type, 
                 transformations=formulationBlack$transformation, 
                 ethnicity="Black", model=modelBlackTrain, 
                 new.data=data_black_train)



## =============================================================================
# Output results of reclassifying overweight

#source("reclassify.R")


# predict training set patients with same ethnicity predictive model
data_white_train$predicted <- 
  exp(predict(object=modelWhiteTrain, 
              newx=makeMatrix(data_white_train, 
                              includeInteraction=formulationWhite$effect=="interaction")$mat,
              type="response")[,"s0"])
data_black_train$predicted <- 
  exp(predict(object=modelBlackTrain, 
              newx=makeMatrix(data_black_train, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])
data_all_train$predicted <- 
  exp(predict(object=modelAllTrain, 
              newx=makeMatrix(data_all_train, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])

# predict test set patients with same ethnicity predictive model
data_white_test$predicted <- 
  exp(predict(object=modelWhiteTrain, 
              newx=makeMatrix(data_white_test, includeInteraction=formulationWhite$effect=="interaction")$mat,
              type="response")[,"s0"])
data_black_test$predicted <- 
  exp(predict(object=modelBlackTrain, 
              newx=makeMatrix(data_black_test, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])
data_all_test$predicted <- 
  exp(predict(object=modelAllTrain, 
              newx=makeMatrix(data_all_test, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])



# data of interest for the two analyses
rocTrainWhite <- roc(data_white_train$ObesityClass, data_white_train$predicted, 
                     levels=c("Normal weight", "Obese"))
cutoffWhite <- coords(rocTrainWhite, x="best")[1,"threshold"]
data_white_test <- data_white_test %>% 
  mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
  mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffWhite)])
aucWhiteTrain <- ci.auc(rocTrainWhite)

rocTrainBlack <- roc(data_black_train$ObesityClass, data_black_train$predicted, 
                     levels=c("Normal weight", "Obese"))
cutoffBlack <- coords(rocTrainBlack, x="best")[1,"threshold"]
data_black_test <- data_black_test %>% 
  mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
  mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffBlack)])
aucBlackTrain <- ci.auc(rocTrainBlack)

rocTrainAll <- roc(data_all_train$ObesityClass, data_all_train$predicted, 
                   levels=c("Normal weight", "Obese"))
cutoffAll <- coords(rocTrainAll, x="best")[1,"threshold"]
data_all_test <- data_all_test %>% 
  mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
  mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffAll)])
aucAllTrain <- ci.auc(rocTrainAll)

plotSelectionROC <- ggplot(bind_rows(data_white_test, data_black_test), 
                           aes(x=predicted, y=exp(BMI))) + 
  geom_point(aes(col=selectROC, pch=Race)) +
  scale_color_manual(values=c("#FF0000", "#00CCCC", "lightgray")) +
  labs(title="Selection for ROC analysis", 
       subtitle=paste0("stratified predictions, effects: ", formulationBlack$effect, 
                       ", type: ", formulationBlack$type,
                       ", BMI transformation: ", formulationBlack$transformation),
       x="predicted BMI", y="observed BMI")
ggsave("reclassify_plotSelectionROC.pdf", plotSelectionROC, width=ggplotWidth, height=ggplotHeight)

plotSelectionReclassify <- ggplot(bind_rows(data_white_test, data_black_test), 
                                  aes(x=predicted, y=exp(BMI))) + 
  geom_point(aes(col=Reclassified, pch=Race)) +
  scale_color_manual(values=c("#00CCCC", "#FF0000")) +
  geom_hline(yintercept=25, lty="dashed") + 
  geom_hline(yintercept=30, lty="dashed") +  
  labs(title="Selection for reclassification analysis",
       subtitle=paste0("stratified predictions, effects: ", formulationBlack$effect, 
                       ", type: ", formulationBlack$type,
                       ", BMI transformation: ", formulationBlack$transformation),
       x="predicted BMI", y="observed BMI")
ggsave("reclassify_plotSelectionReclassify.pdf", plotSelectionReclassify, width=ggplotWidth, height=ggplotHeight)

# ROC curve analysis of observed normal weight versus obese
rocAllBMIclass <- roc(formula=ObesityClass~predicted, data=data_all_test, 
                      levels=c("Normal weight", "Obese"))
aucAll <- ci.auc(rocAllBMIclass)
statsAll <- ci.coords(rocAllBMIclass, x=cutoffAll, input="threshold", 
                      ret=c("sensitivity", "specificity"))

rocWhiteBMIclass <- roc(formula=ObesityClass~predicted, data=data_white_test, 
                        levels=c("Normal weight", "Obese"))
aucWhite <- ci.auc(rocWhiteBMIclass)
statsWhite <- ci.coords(rocWhiteBMIclass, x=cutoffWhite, input="threshold", 
                        ret=c("sensitivity", "specificity"))

rocBlackBMIclass <- roc(formula=ObesityClass~predicted, data=data_black_test, 
                        levels=c("Normal weight", "Obese"))
aucBlack <- ci.auc(rocBlackBMIclass)
statsBlack <- ci.coords(rocBlackBMIclass, x=cutoffBlack, input="threshold", 
                        ret=c("sensitivity", "specificity"))

ROCstats <- matrix(ncol=3, nrow=3)
colnames(ROCstats) <- c("AUC train", "AUC test", "threshold")
rownames(ROCstats) <- c("White", "Black", "All")
ROCstats["White", "threshold"] <- sprintf("%.1f", cutoffWhite)
ROCstats["Black", "threshold"] <- sprintf("%.1f", cutoffBlack)
ROCstats["All", "threshold"] <- sprintf("%.1f", cutoffAll)

ROCstats["White", "AUC train"] <- sprintf("%.2f[%.2f-%.2f]", aucWhiteTrain[2], aucWhiteTrain[1], aucWhiteTrain[3])
ROCstats["Black", "AUC train"] <- sprintf("%.2f[%.2f-%.2f]", aucBlackTrain[2], aucBlackTrain[1], aucBlackTrain[3])
ROCstats["All", "AUC train"] <- sprintf("%.2f[%.2f-%.2f]", aucAllTrain[2], aucAllTrain[1], aucAllTrain[3])

ROCstats["White", "AUC test"] <- sprintf("%.2f[%.2f-%.2f]", aucWhite[2], aucWhite[1], aucWhite[3])
ROCstats["Black", "AUC test"] <- sprintf("%.2f[%.2f-%.2f]", aucBlack[2], aucBlack[1], aucBlack[3])
ROCstats["All", "AUC test"] <- sprintf("%.2f[%.2f-%.2f]", aucAll[2], aucAll[1], aucAll[3])

convertToTexTable(ROCstats, "reclassify_ROCstats.tex", rows.named=TRUE, minipage=TRUE,
                  caption="ROC performance of metabolic BMI predicting observed BMI classes normal weight versus obese, AUC + 95\\% CI.", 
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
ggsave("reclassify_ROCcurves.pdf", ROCcurves, width=ggplotWidth, height=ggplotHeight)

# table of reclassified normal weight and obese among the different observed BMI classes
classCountsWhite <- countClasses(data_white_test, race="White", strat="Stratified")
classCountsBlack <- countClasses(data_black_test, race="Black", strat="Stratified")
#classCountsAll <- countClasses(data_all_test, race="All", strat="Unstratified")

#classCounts <- bind_rows(classCountsWhite, classCountsBlack, classCountsAll)

classCountsStratified <- bind_rows(classCountsWhite, classCountsBlack)

# test for significant difference in reclassified patients between two ethnicities
p_indep <- testReclassificationIndependence(classCountsWhite, classCountsBlack)
convertToTexTable(classCountsStratified, "reclassify_classCounts.tex", minipage=TRUE,
                  caption="Number of patients reclassified as healthy and syndromatic, the fraction reclassified as syndromatic.",
                  reflabel="reclassify_classCounts")

# reclassify patients with unstratified model
rocTrainWhiteUnstratified <- roc(formula=ObesityClass~predicted, 
                                 data=subset(data_all_train, subset = Race=="White"),
                                 levels=c("Normal weight", "Obese"))
cutoffWhiteUnstratified <- coords(rocTrainWhiteUnstratified, x="best")[1,"threshold"]

rocTrainBlackUnstratified <- roc(formula=ObesityClass~predicted, 
                                 data=subset(data_all_train, subset = Race=="Black"), 
                                 levels=c("Normal weight", "Obese"))
cutoffBlackUnstratified <- coords(rocTrainBlackUnstratified, x="best")[1,"threshold"]

data_all_test <- data_all_test %>% 
  mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
  #mutate(cutoff=ifelse(Race=="White", cutoffWhiteUnstratified, cutoffBlackUnstratified)) %>%
  #mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoff)])
  mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffAll)])

# calculate reclassification fractions
classCountsWhite <- 
  countClassesNumeric(data_white_test, race="White", strat="Stratified")
classCountsBlack <- 
  countClassesNumeric(data_black_test, race="Black", strat="Stratified")

classCountsWhiteUnstratified <- 
  countClassesNumeric(subset(data_all_test, subset = Race=="White"), race="White", strat="Unstratified")
classCountsBlackUnstratified <- 
  countClassesNumeric(subset(data_all_test, subset = Race=="Black"), race="Black", strat="Unstratified")


# test independence of reclassification

compareStrat <- matrix(nrow=2, ncol=1)
colnames(compareStrat) <- "strat. indep."
rownames(compareStrat) <- c("White", "Black")
compareStrat["White", "strat. indep."] <- 
  signif(testReclassificationIndependence(classCountsWhite, classCountsWhiteUnstratified), digits=2)
compareStrat["Black", "strat. indep."] <- 
  signif(testReclassificationIndependence(classCountsBlack, classCountsBlackUnstratified), digits=2)

compareRace <- matrix(nrow=2, ncol=1)
colnames(compareRace) <- "ethn. indep."
rownames(compareRace) <- c("unstratified", "stratified")
compareRace["unstratified", "ethn. indep."] <- 
  signif(testReclassificationIndependence(classCountsWhiteUnstratified, 
                                          classCountsBlackUnstratified), digits=2)
compareRace["stratified", "ethn. indep."] <- 
  signif(testReclassificationIndependence(classCountsWhite, classCountsBlack), digits=2)

convertToTexTable(compareStrat, "reclassify_compareStrat.tex", minipage=TRUE, rows.named=TRUE,
                  caption="Stratification independence test p-values.",
                  reflabel="reclassify_compareStrat")
convertToTexTable(compareRace, "reclassify_compareRace.tex", minipage=TRUE, rows.named=TRUE,
                  caption="Ethnicity independence test p-values.",
                  reflabel="reclassify_compareRace")

# visualize computed reclassification fractions
reclassifyFractions <- bind_rows(classCountsWhite, classCountsBlack, 
                                 classCountsWhiteUnstratified, 
                                 classCountsBlackUnstratified) %>%
                         mutate(grouping=paste0(Race, model, sep=""))
p_fractions <- ggplot(data=reclassifyFractions, 
                      aes(x=factor(ObesityClass, levels=c("Normal weight", "Overweight", "Obese")), 
                          y=`frac. syndromatic (point)`,
                          ymin=`frac. syndromatic (lCI)`, 
                          ymax=`frac. syndromatic (uCI)`, 
                          col=Race, pch=model)) +
                 geom_pointrange(stat="identity", position=position_dodge2(width=0.25, padding=0.5)) +
                 geom_line(aes(group=grouping), position=position_dodge2(width=0.25, padding=0.5)) +
                 geom_hline(yintercept=0, lty="dashed") +
                 geom_hline(yintercept=1, lty="dashed") +
                 scale_color_manual(values=c("blue", "red")) +
                 labs(x="", y="fraction reclassified syndromatic", 
                      title="Test set reclassifications")
ggsave("reclassify_fractions.pdf", p_fractions, width=ggplotWidth, height=ggplotHeight)



## =============================================================================
# Output results of analysis on metabolic BMI

#source("metabolicAnalysis.R")

# Calculate scaled model coefficients + 95% CI from bootstrap replication

if (!file.exists("metabolic_CIcoeffsMain.pdf")) {
  
  # generate coefficient estimates for interaction effects modeling
  bootstrapCoeffsWhite <- 
    scaledEffects(data=data_white, effect=formulationWhite$effect, 
                  type=formulationWhite$type, 
                  transformation=formulationWhite$transformation, 
                  balancing="Balanced", boot.n=1000)
  bootstrapCoeffsWhite <- pivot_longer(bootstrapCoeffsWhite, 
                                       cols=colnames(bootstrapCoeffsWhite),
                                       names_to="met", values_to="coeff")
  bootstrapCoeffsBlack <- 
    scaledEffects(data=data_black, effect=formulationBlack$effect, 
                  type=formulationBlack$type, 
                  transformation=formulationBlack$transformation, 
                  balancing="", boot.n=1000)
  bootstrapCoeffsBlack <- pivot_longer(bootstrapCoeffsBlack, 
                                       cols=colnames(bootstrapCoeffsBlack),
                                       names_to="met", values_to="coeff")
  #bootstrapCoeffsAll <- 
  #  scaledEffects(data=data_all, effect=formulationBlack$effect, 
  #                type=formulationBlack$type, 
  #                transformation=formulationBlack$transformation, 
  #                balancing="Balanced", boot.n=1000)
  #bootstrapCoeffsAll <- pivot_longer(bootstrapCoeffsAll, 
  #                                   cols=colnames(bootstrapCoeffsAll),
  #                                   names_to="met", values_to="coeff")
  
  # store plots of beta-coefficients interaction effects modeling
  p_CIcoeffs <- plotScaledEffects(bootstrapCoeffsWhite, bootstrapCoeffsBlack, 
                                  #bootstrapCoeffsAll, 
                                  effect=formulationBlack$effect, 
                                  type=formulationBlack$type, 
                                  transformation=formulationBlack$transformation)
  ggsave("metabolic_CIcoeffsMain.pdf", p_CIcoeffs, width = 6, height = 9)
  
}



## =============================================================================
# Output results of analysis on metabolites related to obesity

#source("outliers.R")

# illustrate the interpretation of metabolite level differences among outliers

illustrateNO <- illustrNO()
illustrateON <- illustrON()

ggsave("outliers_illustrateNO.pdf", illustrateNO, width=ggplotWidth, height=ggplotHeight)
ggsave("outliers_illustrateON.pdf", illustrateON, width=ggplotWidth, height=ggplotHeight)


# calculate cross validated stratified predictions for the full data sets

sampleID <- sample(1:4, size=nrow(data_white), replace=TRUE)
predict_white <- data.frame(ID=c(), predicted=c())
for (i in 1:4) {
  w <- which(sampleID==i)
  modelWhiteCross <- 
    trainRidgeLASSO(effect=formulationWhite$effect, 
                    type=formulationWhite$type, 
                    transformation=formulationWhite$transformation, 
                    new.data=oversample(data_white[-w,]), 
                    lambda=modelWhiteTrain$lambda)
  whitePredictionsCross <- 
    exp(predict(modelWhiteCross, 
                newx=makeMatrix(data_white[w,], 
                                includeInteraction=formulationWhite$effect=="interaction")$mat)[,"s0"])
  predict_white <- bind_rows(predict_white, data.frame(ID=data_white$ID[w], 
                                                  predicted=whitePredictionsCross))
}
data_white <- merge(data_white, predict_white, by="ID")

sampleID <- sample(1:4, size=nrow(data_black), replace=TRUE)
predict_black <- data.frame(ID=c(), predicted=c())
for (i in 1:4) {
  w <- which(sampleID==i)
  modelBlackCross <- 
    trainRidgeLASSO(effect=formulationBlack$effect, 
                    type=formulationBlack$type, 
                    transformation=formulationBlack$transformation, 
                    new.data=data_black[-w,], lambda=modelBlackTrain$lambda)
  blackPredictionsCross <- 
    exp(predict(modelBlackCross, 
                newx=makeMatrix(data_black[w,], 
                                includeInteraction=formulationBlack$effect=="interaction")$mat)[,"s0"])
  predict_black <- bind_rows(predict_black, data.frame(ID=data_black$ID[w], 
                                                  predicted=blackPredictionsCross))
}
data_black <- merge(data_black, predict_black, by="ID")


# the different prediction groups are OO, NN, ON and NO
data_white <- data_white %>% 
  mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<cutoffWhite, "NN", 
                                  ifelse(ObesityClass=="Normal weight"&predicted>cutoffWhite, "NO", 
                                         ifelse(ObesityClass=="Obese"&predicted>cutoffWhite, "OO", 
                                                ifelse(ObesityClass=="Obese"&predicted<cutoffWhite, "ON", NA)))))
data_black <- data_black %>% 
  mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<cutoffBlack, "NN", 
                                  ifelse(ObesityClass=="Normal weight"&predicted>cutoffBlack, "NO", 
                                         ifelse(ObesityClass=="Obese"&predicted>cutoffBlack, "OO", 
                                                ifelse(ObesityClass=="Obese"&predicted<cutoffBlack, "ON", NA)))))

# make plot for illustration of these prediction groups
countPredictionGroup <- bind_rows(table(data_white$predictionGroup), table(data_black$predictionGroup))
countPredictionGroup <- as.matrix(countPredictionGroup)
rownames(countPredictionGroup) <- c("White", "Black")
convertToTexTable(countPredictionGroup, "outliers_predictionGroup.tex", rows.named=TRUE,
                  caption="Number of patients in every prediction group for the ANOVA on metabolic outliers.", 
                  reflabel="outliers_predictionGroup")
illustrateGroups <- ggplot(data=bind_rows(data_white, data_black),
                           aes(x=predicted, y=exp(BMI))) +
                      geom_point(aes(col=predictionGroup, pch=Race)) +
                      labs(title="Metabolic normals (NN, OO) and outliers (NO, ON)", 
                           x="predicted BMI", y="observed BMI")
ggsave("outliers_predictionGroupDemonstration.pdf", illustrateGroups, width=ggplotWidth, height=ggplotHeight)


# The difference in circulating metabolite levels was estimated with a Tukey-corrected ANOVA
whitePlotsDifference <- plotANOVA(data_white, ethnicity="white")
blackPlotsDifference <- plotANOVA(data_black, ethnicity="black")
ggsave("outliers_levelsNOwhite.pdf", whitePlotsDifference$normal, width=ggplotWidth, height=ggplotHeight)
ggsave("outliers_levelsONwhite.pdf", whitePlotsDifference$obese, width=ggplotWidth, height=ggplotHeight)
ggsave("outliers_levelsNOblack.pdf", blackPlotsDifference$normal, width=ggplotWidth, height=ggplotHeight)
ggsave("outliers_levelsONblack.pdf", blackPlotsDifference$obese, width=ggplotWidth, height=ggplotHeight)





## =============================================================================
# Output results of analysis on smoking patients

#source("smoking.R")

# predict smoking patients with same ethnicity predictive model
data_white_smoking$predicted <- 
  exp(predict(object=modelWhiteTrain, 
              newx=makeMatrix(data_white_smoking, 
                              includeInteraction=formulationWhite$effect=="interaction")$mat,
              type="response")[,"s0"])
data_black_smoking$predicted <- 
  exp(predict(object=modelBlackTrain, 
              newx=makeMatrix(data_black_smoking, 
                              includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])

# predict white test set + smoking patients with black predictive model
data_white_test$predictedOther <- 
  exp(predict(object=modelBlackTrain, 
              newx=makeMatrix(data_white_test, includeInteraction=formulationWhite$effect=="interaction")$mat,
              type="response")[,"s0"])
data_white_smoking$predictedOther <- 
  exp(predict(object=modelBlackTrain, 
              newx=makeMatrix(data_white_smoking, includeInteraction=formulationWhite$effect=="interaction")$mat,
              type="response")[,"s0"])

# predict black test set + smoking patients with white predictive model
data_black_test$predictedOther <- 
  exp(predict(object=modelWhiteTrain, 
              newx=makeMatrix(data_black_test, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])
data_black_smoking$predictedOther <- 
  exp(predict(object=modelWhiteTrain, 
              newx=makeMatrix(data_black_smoking, includeInteraction=formulationBlack$effect=="interaction")$mat,
              type="response")[,"s0"])

# comparing non-smokers of the test set with smokers
p_smokeWhite <- ggplot(data=bind_rows(data_white_test, data_white_smoking) %>% 
                         rowwise() %>%
                         mutate(Status=c("non-smoking", "smoking")[1+Smoking]), 
                       aes(x=predicted, y=exp(BMI))) + 
  geom_point(aes(col=Status)) +
  scale_color_manual(values=c("green", "brown")) +
  geom_hline(yintercept=25, lty="dashed") + 
  geom_hline(yintercept=30, lty="dashed") +
  theme(title =element_text(size=16, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=14),
        axis.text.y = element_text(angle = 0, color='black', size=14),
        axis.title.x = element_text(size=12,color='black'),
        axis.title.y = element_text(size=12,color='black')) +
  labs(title="Non-smokers vs. smokers (White)", 
       x="predicted BMI", y="observed BMI")
p_smokeBlack <- ggplot(data=bind_rows(data_black_test, data_black_smoking) %>% 
                         rowwise() %>%
                         mutate(Status=c("non-smoking", "smoking")[1+Smoking]), 
                       aes(x=predicted, y=exp(BMI))) + 
  geom_point(aes(col=Status)) +
  scale_color_manual(values=c("green", "brown")) +
  geom_hline(yintercept=25, lty="dashed") + 
  geom_hline(yintercept=30, lty="dashed") +
  theme(title =element_text(size=16, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=14),
        axis.text.y = element_text(angle = 0, color='black', size=14),
        axis.title.x = element_text(size=12,color='black'),
        axis.title.y = element_text(size=12,color='black')) +
  labs(title="Non-smokers vs. smokers (Black)", 
       x="predicted BMI", y="observed BMI")
ggsave("smoking_compareSmokingWhite.pdf", p_smokeWhite, width=ggplotWidth, height=ggplotHeight)
ggsave("smoking_compareSmokingBlack.pdf", p_smokeBlack, width=ggplotWidth, height=ggplotHeight)


all_data <- bind_rows(data_white_test, data_white_smoking, 
                      data_black_test, data_black_smoking)

countSample <- all_data %>% group_by(Race, ObesityClass, Smoking) %>% summarise(n=n())
countSample <- pivot_wider(countSample, values_from="n", names_from=Smoking)
countSample$`non-smoking` <- countSample$`0`
countSample$`smoking` <- countSample$`1`
countSample$`0` <- NULL
countSample$`1` <- NULL

table_smoke <- all_data %>% group_by(Race, ObesityClass) %>% 
  summarise("pred. BMI (smoking-non smoking)"= confIntMedian(wilcox.test(predicted~Smoking, alternative="two.sided", conf.int=TRUE),
                                                             negative=TRUE), 
            "p-val."=wilcox.test(predicted~Smoking, alternative="two.sided", conf.int=TRUE)$p.value)
table_smoke$Race <- as.character(table_smoke$Race)
table_smoke$ObesityClass <- as.character(table_smoke$ObesityClass)
table_smoke$`p-val.` <- signif(table_smoke$`p-val.`, digits=2)

table_smoke <- merge(table_smoke, countSample, by=c("Race", "ObesityClass"), sort=FALSE)
table_smoke <- table_smoke[,c("Race", "ObesityClass", "non-smoking", "smoking", 
                              "pred. BMI (smoking-non smoking)", "p-val.")]

convertToTexTable(table_smoke, "smoking_compareSmokingTable.tex",
                  caption="Difference in median predicted BMI + 95\\% CI (smoking - non smoking), Wilcoxon test p-value. The prediction model trained on the training data was used to predict test set patients and smoking patients.",
                  reflabel="smoking_compareSmokingTable")

# evaluation of smokers vs. non-smokers using the other ethnicity prediction model
p_smokeWhiteOther <- ggplot(data=bind_rows(data_white_test, data_white_smoking) %>% 
                              rowwise() %>%
                              mutate(Status=c("non-smoking", "smoking")[1+Smoking]), 
                            aes(x=predictedOther, y=exp(BMI))) + 
  geom_point(aes(col=Status)) +
  scale_color_manual(values=c("green", "brown")) +
  geom_hline(yintercept=25, lty="dashed") + 
  geom_hline(yintercept=30, lty="dashed") +
  theme(title =element_text(size=16, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=14),
        axis.text.y = element_text(angle = 0, color='black', size=14),
        axis.title.x = element_text(size=12,color='black'),
        axis.title.y = element_text(size=12,color='black')) +
  labs(title="Non-smokers vs. smokers (White)", 
       subtitle="black predictive model",
       x="predicted BMI", y="observed BMI")
p_smokeBlackOther <- ggplot(data=bind_rows(data_black_test, data_black_smoking) %>% 
                              rowwise() %>%
                              mutate(Status=c("non-smoking", "smoking")[1+Smoking]), 
                            aes(x=predictedOther, y=exp(BMI))) + 
  geom_point(aes(col=Status)) +
  scale_color_manual(values=c("green", "brown")) +
  geom_hline(yintercept=25, lty="dashed") + 
  geom_hline(yintercept=30, lty="dashed") +
  theme(title =element_text(size=16, face='bold'),
        axis.text.x = element_text(angle = 0, color='black', size=14),
        axis.text.y = element_text(angle = 0, color='black', size=14),
        axis.title.x = element_text(size=12,color='black'),
        axis.title.y = element_text(size=12,color='black')) +
  labs(title="Non-smokers vs. smokers (Black)", 
       subtitle="white predictive model",
       x="predicted BMI", y="observed BMI")
ggsave("smoking_compareSmokingWhiteOther.pdf", p_smokeWhiteOther, width=ggplotWidth, height=ggplotHeight)
ggsave("smoking_compareSmokingBlackOther.pdf", p_smokeBlackOther, width=ggplotWidth, height=ggplotHeight)

table_smoke_other <- all_data %>% group_by(Race, ObesityClass) %>% 
  summarise("pred. BMI (smoking-non smoking)"= confIntMedian(wilcox.test(predictedOther~Smoking, alternative="two.sided", conf.int=TRUE),
                                                             negative=TRUE), 
            "p-val."=wilcox.test(predictedOther~Smoking, alternative="two.sided", conf.int=TRUE)$p.value)
table_smoke_other$Race <- as.character(table_smoke_other$Race)
table_smoke_other$ObesityClass <- as.character(table_smoke_other$ObesityClass)
table_smoke_other$`p-val.` <- signif(table_smoke_other$`p-val.`, digits=2)

table_smoke_other <- merge(table_smoke_other, countSample, by=c("Race", "ObesityClass"), sort=FALSE)
table_smoke_other <- table_smoke_other[,c("Race", "ObesityClass", "non-smoking", "smoking", 
                              "pred. BMI (smoking-non smoking)", "p-val.")]

convertToTexTable(table_smoke_other, "smoking_compareSmokingTableOther.tex",
                  caption="Difference in median predicted BMI + 95\\% CI (smoking - non smoking), Wilcoxon test p-value. The prediction model trained on the training data of the other ethnicity was used to predict test set patients and smoking patients.",
                  reflabel="smoking_compareSmokingTableOther")

# metabolic differences between smokers and non-smokers
whiteSmokingContr <- plotSmokingContributions(data_white_smoking, data_white_test, modelWhiteTrain) +
  labs(title="Contribution to prediction shift smokers", subtitle="White ethnicity")
blackSmokingContr <- plotSmokingContributions(data_black_smoking, data_black_test, modelBlackTrain) +
  labs(title="Contribution to prediction shift smokers", subtitle="Black ethnicity")

ggsave("smoking_fractionWhite.pdf", whiteSmokingContr, width=ggplotWidth, height=ggplotHeight)
ggsave("smoking_fractionBlack.pdf", blackSmokingContr, width=ggplotWidth, height=ggplotHeight)


## =============================================================================
# Output results of discussion

#source("discussion.R")
  
# is AUC of test set predictions inferior to AUC of training set?
# post-hoc analysis on data splitting

if (!file.exists("posthocDataWhite.rds")) {
  aucTrain <- vector(mode="numeric", length=100)
  aucTest <- vector(mode="numeric", length=100)
  numbersObeseTrain <- c()
  numbersObeseTest <- c()
  meanAgeTrain <- vector(mode="numeric", length=100)
  meanAgeTrainNormal <- vector(mode="numeric", length=100)
  meanAgeTrainObese <- vector(mode="numeric", length=100)
  meanAgeTest <- vector(mode="numeric", length=100)
  meanAgeTestNormal <- vector(mode="numeric", length=100)
  meanAgeTestObese <- vector(mode="numeric", length=100)
  
  for (i in c(1:101)) {
    if (i==101) {
      set.seed(23/8/2024)
    } else {
      set.seed(i)
    }
    valSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
    data_white_test <- data_white[valSelectA,]
    data_white_train <- data_white[-valSelectA,]
    modelWhiteTrain <- 
      trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                      transformation=formulationWhite$transformation, 
                      new.data=oversample(data_white_train))
    data_white_train$predicted <- 
      exp(predict(object=modelWhiteTrain, 
                  newx=makeMatrix(data_white_train, 
                                  includeInteraction=formulationWhite$effect=="interaction")$mat,
                  type="response")[,"s0"])
    data_white_test$predicted <- 
      exp(predict(object=modelWhiteTrain, 
                  newx=makeMatrix(data_white_test, includeInteraction=formulationWhite$effect=="interaction")$mat,
                  type="response")[,"s0"])
    rocTrainWhite <- roc(data_white_train$ObesityClass, data_white_train$predicted, 
                         levels=c("Normal weight", "Obese"))
    cutoffWhite <- coords(rocTrainWhite, x="best")[1,"threshold"]
    data_white_test <- data_white_test %>% 
      mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
      mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffWhite)])
    rocWhiteBMIclass <- roc(formula=ObesityClass~predicted, data=data_white_test, 
                            levels=c("Normal weight", "Obese"))
    aucTest[i] <- auc(rocWhiteBMIclass)
    statsWhite <- ci.coords(rocWhiteBMIclass, x=cutoffWhite, input="threshold", 
                            ret=c("sensitivity", "specificity"))
    aucTrain[i] <- auc(rocTrainWhite)
    
    classCountsTrain <- table(data_white_train$ObesityClass)
    numbersObeseTrain <- bind_rows(numbersObeseTrain, classCountsTrain)
    
    classCountsTest <- table(data_white_test$ObesityClass)
    numbersObeseTest <- bind_rows(numbersObeseTest, classCountsTest)
    
    meanAgeTrain[i] <- mean(data_white_train$Age)
    meanAgeTest[i] <- mean(data_white_test$Age)
    
    meanAgeTrainNormal[i] <- mean(subset(data_white_train, subset=ObesityClass=="Normal weight")$Age)
    meanAgeTestNormal[i] <- mean(subset(data_white_test, subset=ObesityClass=="Normal weight")$Age)
    
    meanAgeTrainObese[i] <- mean(subset(data_white_train, subset=ObesityClass=="Obese")$Age)
    meanAgeTestObese[i] <- mean(subset(data_white_test, subset=ObesityClass=="Obese")$Age)
    
  }
  posthocDataWhite <- bind_cols("AUC train"=aucTrain, 
                                "AUC test"=aucTest,
                                numbersObeseTrain, numbersObeseTest,
                                "age train"=meanAgeTrain,
                                "age test"=meanAgeTest,
                                "age train normal"=meanAgeTrainNormal,
                                "age test normal"=meanAgeTestNormal,
                                "age train obese"=meanAgeTrainObese,
                                "age test obese"=meanAgeTestObese,
                                seeds=c(rep("1 to 100", times=100), "23/8/2024"),
                                check.names=FALSE)
  saveRDS(posthocDataWhite, "posthocDataWhite.rds")
}

if (!file.exists("posthocDataBlack.rds")) {
  aucTrain <- vector(mode="numeric", length=100)
  aucTest <- vector(mode="numeric", length=100)
  numbersObeseTrain <- c()
  numbersObeseTest <- c()
  meanAgeTrain <- vector(mode="numeric", length=100)
  meanAgeTrainNormal <- vector(mode="numeric", length=100)
  meanAgeTrainObese <- vector(mode="numeric", length=100)
  meanAgeTest <- vector(mode="numeric", length=100)
  meanAgeTestNormal <- vector(mode="numeric", length=100)
  meanAgeTestObese <- vector(mode="numeric", length=100)
  
  for (i in c(1:101)) {
    if (i==101) {
      set.seed(23/8/2024)
    } else {
      set.seed(i)
    }
    valSelectA <- sample(x=1:nrow(data_black), size=round(0.2*nrow(data_black)))
    data_black_test <- data_black[valSelectA,]
    data_black_train <- data_black[-valSelectA,]
    modelblackTrain <- 
      trainRidgeLASSO(effect=formulationBlack$effect, type=formulationBlack$type, 
                      transformation=formulationBlack$transformation, 
                      new.data=data_black_train)
    data_black_train$predicted <- 
      exp(predict(object=modelBlackTrain, 
                  newx=makeMatrix(data_black_train, 
                                  includeInteraction=formulationBlack$effect=="interaction")$mat,
                  type="response")[,"s0"])
    data_black_test$predicted <- 
      exp(predict(object=modelBlackTrain, 
                  newx=makeMatrix(data_black_test, includeInteraction=formulationBlack$effect=="interaction")$mat,
                  type="response")[,"s0"])
    rocTrainBlack <- roc(data_black_train$ObesityClass, data_black_train$predicted, 
                         levels=c("Normal weight", "Obese"))
    cutoffBlack <- coords(rocTrainBlack, x="best")[1,"threshold"]
    data_black_test <- data_black_test %>% 
      mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
      mutate(Reclassified=c("Healthy", "Syndromatic")[1+(predicted>cutoffBlack)])
    rocBlackBMIclass <- roc(formula=ObesityClass~predicted, data=data_black_test, 
                            levels=c("Normal weight", "Obese"))
    aucTest[i] <- auc(rocBlackBMIclass)
    statsBlack <- ci.coords(rocBlackBMIclass, x=cutoffBlack, input="threshold", 
                            ret=c("sensitivity", "specificity"))
    aucTrain[i] <- auc(rocTrainBlack)
    
    classCountsTrain <- table(data_black_train$ObesityClass)
    numbersObeseTrain <- bind_rows(numbersObeseTrain, classCountsTrain)
    
    classCountsTest <- table(data_black_test$ObesityClass)
    numbersObeseTest <- bind_rows(numbersObeseTest, classCountsTest)
    
    meanAgeTrain[i] <- mean(data_black_train$Age)
    meanAgeTest[i] <- mean(data_black_test$Age)
    
    meanAgeTrainNormal[i] <- mean(subset(data_black_train, subset=ObesityClass=="Normal weight")$Age)
    meanAgeTestNormal[i] <- mean(subset(data_black_test, subset=ObesityClass=="Normal weight")$Age)
    
    meanAgeTrainObese[i] <- mean(subset(data_black_train, subset=ObesityClass=="Obese")$Age)
    meanAgeTestObese[i] <- mean(subset(data_black_test, subset=ObesityClass=="Obese")$Age)
    
  }
  posthocDataBlack <- bind_cols("AUC train"=aucTrain, 
                                "AUC test"=aucTest,
                                numbersObeseTrain, numbersObeseTest,
                                "age train"=meanAgeTrain,
                                "age test"=meanAgeTest,
                                "age train normal"=meanAgeTrainNormal,
                                "age test normal"=meanAgeTestNormal,
                                "age train obese"=meanAgeTrainObese,
                                "age test obese"=meanAgeTestObese,
                                seeds=c(rep("1 to 100", times=100), "23/8/2024"),
                                check.names=FALSE)
  saveRDS(posthocDataBlack, "posthocDataBlack.rds")
}

posthocDataWhite <- readRDS("posthocDataWhite.rds")

p_posthocWhite <- ggplot(data=posthocDataWhite,
                         mapping=aes(x=`AUC train`, y=`AUC test`, col=seeds)) +
                    scale_fill_continuous(type="gradient", palette="PuBuGn") +
                    geom_point() +
                    geom_abline(slope=1, intercept=0, lty="dashed")

posthocDataBlack <- readRDS("posthocDataBlack.rds")

p_posthocBlack <- ggplot(data=posthocDataBlack,
                         mapping=aes(x=`AUC train`, y=`AUC test`, col=seeds)) +
                    scale_fill_continuous(type="gradient", palette="PuBuGn") +
                    geom_point() +
                    geom_abline(slope=1, intercept=0, lty="dashed")

ggsave("posthoc_white.pdf", p_posthocWhite, width=ggplotWidth, height=ggplotHeight)
ggsave("posthoc_black.pdf", p_posthocBlack, width=ggplotWidth, height=ggplotHeight)

# store the median AUC values + 95pc CI of the post-hoc analysis
TrainTestAUC <- matrix(nrow=2, ncol=2)
colnames(TrainTestAUC) <- c("White", "Black")
rownames(TrainTestAUC) <- c("AUC train", "AUC test")

TrainTestAUC["AUC train", "White"] <- sprintf("%.3f[%.3f %.3f]",
                                              quantile(posthocDataWhite$`AUC train`, probs=0.5),
                                              quantile(posthocDataWhite$`AUC train`, probs=0.025),
                                              quantile(posthocDataWhite$`AUC train`, probs=0.975))
TrainTestAUC["AUC test", "White"] <- sprintf("%.3f[%.3f %.3f]",
                                              quantile(posthocDataWhite$`AUC test`, probs=0.5),
                                              quantile(posthocDataWhite$`AUC test`, probs=0.025),
                                              quantile(posthocDataWhite$`AUC test`, probs=0.975))

TrainTestAUC["AUC train", "Black"] <- sprintf("%.3f[%.3f %.3f]",
                                              quantile(posthocDataBlack$`AUC train`, probs=0.5),
                                              quantile(posthocDataBlack$`AUC train`, probs=0.025),
                                              quantile(posthocDataBlack$`AUC train`, probs=0.975))
TrainTestAUC["AUC test", "Black"] <- sprintf("%.3f[%.3f %.3f]",
                                             quantile(posthocDataBlack$`AUC test`, probs=0.5),
                                             quantile(posthocDataBlack$`AUC test`, probs=0.025),
                                             quantile(posthocDataBlack$`AUC test`, probs=0.975))

convertToTexTable(TrainTestAUC, "posthoc_AUC.tex",
                  caption="Post-hoc estimates of the AUC values for train and test set. Point estimate + 95\\% CI",
                  reflabel="posthoc_AUC",
                  rows.named=TRUE)

