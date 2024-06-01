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
valSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
data_white_val <- data_white[valSelectA,]
data_white_use <- data_white[-valSelectA,]
set.seed(3)
trainSelectA <- sample(x=1:nrow(data_white_use), size=round(0.75*nrow(data_white_use)))
data_white_train <- data_white_use[trainSelectA,]
data_white_test <- data_white_use[-trainSelectA,]


data_white_train_balanced <- oversample(data_white_use)


# black ethnicity data

data_black <- subset(data_regular, subset = Race=="Black")
set.seed(2)
valSelectB <- sample(x=1:nrow(data_black), size=round(0.2*nrow(data_black)))
data_black_val <- data_black[valSelectB,]
data_black_use <- data_black[-valSelectB,]
set.seed(3)
trainSelectB <- sample(x=1:nrow(data_black_use), size=round(0.75*nrow(data_black_use)))
data_black_train <- data_black_use[trainSelectB,]
data_black_test <- data_black_use[-trainSelectB,]


# tabulate cross-validated performance of all model recipes

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

if (!file.exists("trainWhite.tex")) {
  trainWhite <- tabulateValidation(effects, types, transformations, ethnicity="White", balancing, data_white_use)
  trainBlack <- tabulateValidation(effects1, types1, transformations1, ethnicity="Black", new.data=data_black_use)
  convertToTexTable(trainWhite, "trainWhite.tex", 
                    caption="Cross-validated prediction evaluation of the models on white ethnicity training data.", 
                    reflabel="trainWhite")
  convertToTexTable(trainBlack, "trainBlack.tex", 
                    caption="Cross-validated prediction evaluation of the models on black ethnicity training data.", 
                    reflabel="trainBlack")
}

# detection of metabolic outliers
# From the cross-validation, model recipes interactionOLSInvModelWhiteBalanced for white and 
# interactionLASSOInvModelBlack for black were selected as the most appropriate. 
# Train models on all data for retrieving the metabolic outliers

modelWhiteFull <- 
  trainLASSORidge(effect="interaction", type="LASSO", transformation="Inv", new.data=oversample(data_white))
modelBlackFull <- 
  trainLASSORidge(effect="interaction", type="LASSO", transformation="Inv", new.data=data_black)

allPredsWhite <- 1/predict(object=modelWhiteFull, 
                           newx=makeMatrix(data_white, includeInteraction=TRUE)$mat,
                           type="response")[,"s0"]
allPredsClassWhite <- c("predicted: Normal weight", "predicted: Overweight", 
                        "predicted: Obese") [1 + (allPredsWhite>25) + (allPredsWhite>30)]
allPredsClassWhite <- factor(allPredsClassWhite, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
allObsClassWhite <- data_white$ObesityClass
allPredictTableWhite <- table(allPredsClassWhite, allObsClassWhite)
convertToTexTable(allPredictTableWhite, "allPredictTableWhite.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for all white ethnicity patients.",
                  reflabel="allPredictTableWhite")

allPredsBlack <- 1/predict(object=modelBlackFull, 
                           newx=makeMatrix(data_black, includeInteraction=TRUE)$mat,
                           type="response")[,"s0"]
allPredsClassBlack <- c("predicted: Normal weight", "predicted: Overweight", 
                        "predicted: Obese") [1 + (allPredsBlack>25) + (allPredsBlack>30)]
allPredsClassBlack <- factor(allPredsClassBlack, levels=c("predicted: Normal weight", "predicted: Overweight", "predicted: Obese"))
allObsClassBlack <- data_black$ObesityClass
allPredictTableBlack <- table(allPredsClassBlack, allObsClassBlack)
convertToTexTable(allPredictTableBlack, "allPredictTableBlack.tex", rows.named=TRUE,
                  caption="Contingency table of observed and predicted BMI class for all black ethnicity patients.",
                  reflabel="allPredictTableBlack")

# make posters of model diagnostics of chosen models

finalWhite <- 
  trainLASSORidge(effect="interaction", type="LASSO", transformation="Inv", 
                  new.data=oversample(data_white_use))
finalBlack <- 
  trainLASSORidge(effect="interaction", type="LASSO", transformation="Inv", 
                  new.data=data_black_use)

interactionLASSOInvModelWhiteBalanced <- finalWhite
interactionLASSOInvModelBlack <- finalBlack
modelDiagnostics(effects="interaction", types="LASSO", transformations="Inv", balancing="Balanced", 
                 ethnicity="White", new.data=data_white_use)
modelDiagnostics(effects="interaction", types="LASSO", transformations="Inv", 
                 ethnicity="Black", new.data=data_black_use)




## =============================================================================
# Output results of metabolic outlier analysis

#source("metabolicAnalysis.R")


# calculate power and sample size of fitted metabolite beta coefficients

whitePredictions <- 1/predict(finalWhite, newx=makeMatrix(data_white_use, includeInteraction=TRUE)$mat)[,"s0"]
blackPredictions <- 1/predict(finalBlack, newx=makeMatrix(data_black_use, includeInteraction=TRUE)$mat)[,"s0"]

if (any(!file.exists(c("whiteSimulation.rds", "blackSimulation.rds")))) {
  whiteSimulation <- simulateAlternative(finalWhite, data_white_use, 1/whitePredictions,
                                         transformation="Inv", type="LASSO", nsim=1000, oversampled=TRUE)
  blackSimulation <- simulateAlternative(finalBlack, data_black_use, 1/blackPredictions,
                                         transformation="Inv", type="Ridge", nsim=1000)
  
  saveRDS(whiteSimulation, "whiteSimulation.rds")
  saveRDS(blackSimulation, "blackSimulation.rds")
}

whiteSimulation <- readRDS("whiteSimulation.rds")
blackSimulation <- readRDS("blackSimulation.rds")


# calculate power from the standard deviations of the simulated coefficients


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
                            x="Normal quantile", y="Standardized beta coefficient")

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
                            x="Normal quantile", y="Standardized beta coefficient")

ggsave(filename="simBetasWhiteQQplot.pdf", plotSimBetasWhite)
ggsave(filename="simBetasBlackQQplot.pdf", plotSimBetasBlack)


# use this assumption to perform a sample size calculation
sampleSizeWhite <- sampleSizeObeseMetabolites(betas=whiteSimulation$betas,
                                              sample_sizes=c(50, 60, 80, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000), 
                                              n_sample=nrow(data_white_use), alpha=0.05)
sampleSizeBlack <- sampleSizeObeseMetabolites(betas=blackSimulation$betas,
                                              sample_sizes=c(50, 60, 80, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000), 
                                              n_sample=nrow(data_black_use), alpha=0.05)

whitePower <- plotSampleSize(sampleSizeWhite, requested_power=0.75, at_samplesize=1000, plottitle="Power of predictors (White ethnicity)")
blackPower <- plotSampleSize(sampleSizeBlack, requested_power=0.75, at_samplesize=1000, plottitle="Power of predictors (Black ethnicity)")

ggsave(filename="whitePower.pdf", whitePower)
ggsave(filename="blackPower.pdf", blackPower)


# The analysis on metabolites related to obesity needs to be checked with the left-out validation data


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
