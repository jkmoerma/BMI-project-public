
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


