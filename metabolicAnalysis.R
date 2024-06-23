
scaledEffects <- function(data, effect, type, transformation, balancing="", boot.n=100) {
  
  # if balancing the data set was part of the model formulation, scaling data on an oversampled replicate was necessary for ...
  if (balancing=="") {data_ref <- data}
  if (balancing=="Balanced") {data_ref <- oversample(data)}
  
  # generate a scaled version of the data set 
  data_scaled <- data
  for (met in metabolites) {
    data_scaled[[met]] <- data_scaled[[met]]/sd(data_ref[[met]])
  }
  
  # standardize BMI on type of transformation
  if (transformation=="Log") {data_scaled[["BMI"]] <- data_scaled[["BMI"]]/sd(data_ref[["BMI"]])}
  if (transformation=="Inv") {
    stdev <- sd(exp(-data_ref[["BMI"]]))
    transformScaled <- exp(-data_scaled[["BMI"]])/stdev
    data_scaled[["BMI"]] <- -log(transformScaled)
  }
  dataMatrix <- makeMatrix(data_scaled, includeInteraction = effect=="interaction")$mat
  
  # for boot.n bootstrap replicates, recalculate the regression coefficients
  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Parallelized loop
  coeffs <- foreach(i=1:boot.n, .combine=bind_rows, .packages=c("dplyr", "glmnet"), 
                    .export=c("trainOLS", "trainRidgeLASSO", 
                              "metabolites", "riskLevel", "aic", "oversample", 
                              "SMOTE", "makeMatrix")) %dopar% {
    set.seed(i)
    # retrieve stratified bootstrap replicate
    w <- sample(1:nrow(data_scaled), size=nrow(data_scaled), replace=TRUE)
    data_rep <- data_scaled[w,]
    
    if (type=="OLS") {
      if (balancing=="") {
        model <- trainOLS(effect, type, transformation, data_rep)
      }
      if (balancing=="Balanced") {
        model <- trainOLS(effect, type, transformation, oversample(data_rep))
      }
      return(model$coefficients)
    }
    if (type=="LASSO"|type=="Ridge") {
      if (balancing=="") {
        new.x <- dataMatrix[w,]
        new.y <- data_scaled$BMI[w]
        model <- trainRidgeLASSO(effect, type, transformation, new.x=new.x, new.y=new.y)
      }
      if (balancing=="Balanced") {
        model <- trainRidgeLASSO(effect, type, transformation, oversample(data_rep))
      }
      return(model$beta[,"s0"])
    }
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  return(coeffs)
  
}

plotANOVA <- function(data, ethnicity) {
  groupLevels <- subset(data, subset=!is.na(predictionGroup), 
                             select=c(metabolites, "predictionGroup"))
  
  levelDiffsNormal <- matrix(nrow=length(metabolites), ncol=6)
  rownames(levelDiffsNormal) <- metabolites
  colnames(levelDiffsNormal) <- c("NO-NN", "NO-NN (lCI)", "NO-NN (uCI)",
                                       "OO-NO", "OO-NO (lCI)", "OO-NO (uCI)")
  
  levelDiffsObese <- matrix(nrow=length(metabolites), ncol=6)
  rownames(levelDiffsObese) <- metabolites
  colnames(levelDiffsObese) <- c("OO-ON", "OO-ON (lCI)", "OO-ON (uCI)",
                                      "NN-ON", "NN-ON (lCI)", "NN-ON (uCI)")
  
  legend <- c()
  for (met in metabolites) {
    diffs <- TukeyHSD(aov(scale(groupLevels[[met]])~groupLevels$predictionGroup))$`groupLevels$predictionGroup`
    levelDiffsNormal[met,] <- c(diffs["NO-NN", "diff"], diffs["NO-NN", "lwr"], 
                                     diffs["NO-NN", "upr"], diffs["OO-NO", "diff"], 
                                     diffs["OO-NO", "lwr"], diffs["OO-NO", "upr"])
    levelDiffsObese[met,] <- c(diffs["OO-ON", "diff"], diffs["OO-ON", "lwr"], 
                                    diffs["OO-ON", "upr"], -diffs["ON-NN", "diff"], 
                                    -diffs["ON-NN", "upr"], -diffs["ON-NN", "lwr"])
    if (abs(diffs["OO-NN", "diff"])>0.5) {
      legend <- c(legend, met)
    } else {
      legend <- c(legend, "other")
    }
  }
  
  levelDiffsNormal <- data.frame(met=rownames(levelDiffsNormal),
                                      legend = legend,
                                      as.data.frame(levelDiffsNormal),
                                      check.names=FALSE)
  levelDiffsObese <- data.frame(met=rownames(levelDiffsObese),
                                     legend=legend,
                                     as.data.frame(levelDiffsObese),
                                     check.names=FALSE)
  
  pN <- ggplot(data=levelDiffsNormal, 
               aes(x=`OO-NO`, y=`NO-NN`, label=met, col=met)) + 
          theme(legend.position = "none") +
          coord_fixed(ratio=1) +
          geom_abline(slope=0, intercept=0, lty="dashed") +
          geom_vline(xintercept=0, lty="dashed") +
          geom_pointrange(aes(ymin=`NO-NN (lCI)`, ymax=`NO-NN (uCI)`)) +
          geom_pointrange(aes(xmin=`OO-NO (lCI)`, xmax=`OO-NO (uCI)`)) +
          geom_text(vjust=0, nudge_y=0.05, size=2) +
          labs(title=paste0("Normal weight with obese metabolome (NO) (", ethnicity, " ethnicity)"), 
               x="Scaled log conc. diff. OO-NO", 
               y="Scaled log conc. diff. NO-NN")
  
  pO <- ggplot(data=levelDiffsObese, 
               aes(x=`NN-ON`, y=`OO-ON`, label=met, col=met)) + 
          theme(legend.position = "none") +
          coord_fixed(ratio=1) + 
          geom_abline(slope=0, intercept=0, lty="dashed") +
          geom_vline(xintercept=0, lty="dashed") +
          geom_pointrange(aes(ymin=`OO-ON (lCI)`, ymax=`OO-ON (uCI)`)) +
          geom_pointrange(aes(xmin=`NN-ON (lCI)`, xmax=`NN-ON (uCI)`)) +
          geom_text(vjust=0, nudge_y=0.05, size=2) +
          labs(title=paste0("Obese with normal metabolome (ON) (", ethnicity, " ethnicity)"), 
               x="Scaled log conc. diff. NN-ON", 
               y="Scaled log conc. diff. OO-ON")
  
  list(normal=pN, obese=pO)
}
