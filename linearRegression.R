aic <- function(model, oversampled=FALSE, amount=0) {
  n <- length(model$coefficients)+model$df.residual
  SSE <- sum(model$residuals**2)
  p <- length(model$coefficients)
  
  if (oversampled) {
    n <- n-amount
    SSE <- SSE*n/(n+amount)
  }
  
  n*log(SSE) - n*log(n) + 2*p
}

oversample <- function(data_train) {
  
  data_train$ID <- NULL
  data_train$Smoking <- as.factor(data_train$Smoking)
  
  overweight_data <- subset(data_train, subset = !ObesityClass=="Obese")
  overweight_data$ObesityClass <- factor(ifelse(overweight_data$ObesityClass == "Overweight", 
                                                "Overweight","Normal weight")) 
  counts <- table(overweight_data$ObesityClass)
  minority <- min(counts)
  majority <- max(counts)
  overweight_resampled <- SMOTE(ObesityClass ~ met_110, overweight_data, 
                                perc.over=100*(majority/minority-1), 
                                perc.under=100/(1-minority/majority))
  
  obese_data <- subset(data_train, subset = !ObesityClass=="Overweight")
  obese_data$ObesityClass <- factor(ifelse(obese_data$ObesityClass == "Obese", 
                                           "Obese","Normal weight")) 
  counts <- table(obese_data$ObesityClass)
  minority <- min(counts)
  majority <- max(counts)
  obese_resampled <- SMOTE(ObesityClass ~ met_110, obese_data, 
                           perc.over=100*(majority/minority-1), 
                           perc.under=100/(1-minority/majority))
  
  data_balanced <- bind_rows(subset(data_train, subset = ObesityClass=="Normal weight"),
                             subset(overweight_resampled, subset = ObesityClass=="Overweight"), 
                             subset(obese_resampled, subset = ObesityClass=="Obese"))
  data_balanced$Smoking <- as.numeric(data_balanced$Smoking)-1
  data_train$Smoking <- as.numeric(data_train$Smoking)-1
  
  # recalculate metabolite ratios ON LOG SCALE !!
  ratios <- metabolites[which(grepl(pattern="/", metabolites))]
  for (ratio in ratios) {
    mets <- strsplit(ratio, split="/")[[1]]
    data_balanced[[ratio]] <- with(data_balanced, eval(parse(text=paste(mets[1], "-", mets[2]))))
  }
  
  data_balanced
  
}

makeMatrix <- function(data, includeInteraction=FALSE) {
  irrelevant <- which(colnames(data)%in%c("Race", "ID", "ObesityClass", "BMI", "transBMI", "predicted"))
  RidgeLASSO_data <- as.matrix(data[,-irrelevant])
  vars <- colnames(data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (includeInteraction) {
    for (i in 1:(length(vars)-1)) {
      for (j in (i+1):length(vars)) {
        var1 <- vars[i]
        var2 <- vars[j]
        RidgeLASSO_data <- cbind(RidgeLASSO_data, RidgeLASSO_data[,var1]*RidgeLASSO_data[,var2])
        interactions <- c(interactions, paste(var1, var2, sep="*"))
      }
    }
    colnames(RidgeLASSO_data) <- c(vars, interactions)
  }
  return(list(mat=RidgeLASSO_data, interactions=interactions))
}


trainOLS <- function(effect, type, transformation, balancing="", new.data) {
  
  n0 <- nrow(new.data)
  if (balancing=="Balanced") {new.data <- oversample(new.data)}
  amount <- nrow(new.data) - n0
  
  new.data$transBMI <- new.data$BMI
  if (transformation=="Inv") {new.data$transBMI <- exp(-new.data$BMI)}
  
  interactionEffect <- FALSE
  if (effect=="interaction") {
    interactionEffect <- TRUE
  }
  
  # main effects, log(BMI), ethnicity="White"
  intercept_only <- lm(transBMI ~ 1, data=new.data)
  all <- lm(data=new.data,
            formula = eval(parse(text=paste0("transBMI~`", paste(c("Age", "Smoking", metabolites), 
                                                            collapse="`+`"), "`"))))
  model <- step(intercept_only, direction='both', scope=formula(all), trace=0,
                k=2*(n0+amount)/n0)
  
  if (interactionEffect) {
    # interaction effects, log(BMI), ethnicity="White"
    variables <- 
      colnames(model$model)[-which(colnames(model$model)=="transBMI")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=new.data,
                         formula=eval(parse(text=paste0("transBMI~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        AIC0 <- aic(model, oversampled=balancing=="Balanced", amount=amount)
        AIC1 <- aic(addedModel, oversampled=balancing=="Balanced", amount=amount)
        if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=new.data,
              formula = eval(parse(text=paste0("transBMI~`", paste0(c(variables, interactionSet), 
                                                               collapse="`+`"), "`"))))
    model <- step(model, direction='both', scope=formula(all), trace=0,
                  k=2*(n0+amount)/n0)
  }
  model
}

trainLASSORidge <- function(effect, type, transformation, new.data) {
  
  interactionEffect <- FALSE
  if (effect=="interaction") {interactionEffect <- TRUE}
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI", "transBMI"))
  LASSO_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (interactionEffect) {
    for (i in 1:(length(vars)-1)) {
      for (j in (i+1):length(vars)) {
        var1 <- vars[i]
        var2 <- vars[j]
        LASSO_data <- cbind(LASSO_data, LASSO_data[,var1]*LASSO_data[,var2])
        interactions <- c(interactions, paste(var1, var2, sep="*"))
      }
    }
    colnames(LASSO_data) <- c(vars, interactions)
  }
  
  outcomes <- new.data$BMI
  if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  
  # divide LASSO_data in 4 folds for cross-validation
  set.seed(4)
  foldid <- sample(1:4, size=nrow(LASSO_data), replace=TRUE)
  
  # stabilize probed lambda values with first validation + initiate IQRs
  w <- which(foldid==1)
  LASSOModel <- glmnet(x=LASSO_data[-w, c(vars, interactions)],
                       y=outcomes[-w],
                       alpha=alpha,
                       family="gaussian")
  LASSOPreds <- predict(object=LASSOModel, 
                        newx=LASSO_data[w, c(vars, interactions)],
                        type="response")
  LASSOOuts <- matrix(rep(outcomes[w], times=LASSOModel$dim[2]), 
                      byrow=FALSE, ncol=LASSOModel$dim[2])
  LASSORes <- LASSOPreds - LASSOOuts
  IQRs <- sapply(X=1:LASSOModel$dim[2], 
                 FUN=function(j) {quants <- quantile(LASSORes[,j], probs=c(0.25, 0.75));
                 quants[2] - quants[1]}
  )
  lambdas <- LASSOModel$lambda
  
  # iterate over other cross folds
  for (i in 2:4) {
    w <- which(foldid==i)
    LASSOModel <- glmnet(x=LASSO_data[-w, c(vars, interactions)],
                         y=outcomes[-w],
                         alpha=alpha,
                         lambda=lambdas,
                         family="gaussian")
    LASSOPreds <- predict(object=LASSOModel, 
                          newx=LASSO_data[w, c(vars, interactions)],
                          type="response")
    LASSOOuts <- matrix(rep(outcomes[w], times=LASSOModel$dim[2]), 
                        byrow=FALSE, ncol=LASSOModel$dim[2])
    LASSORes <- LASSOPreds - LASSOOuts
    IQRi <- sapply(X=1:LASSOModel$dim[2], 
                   FUN=function(j) {quants <- quantile(LASSORes[,j], probs=c(0.25, 0.75));
                   quants[2] - quants[1]}
                   )
    IQRs <- IQRs + IQRi
  }
  IQRs <- IQRs/4
  
  # choose lambda parameter tuned to the smallest inter-quartile range in validation set
  tuneIndex <- which.min(IQRs)
  lambda <- LASSOModel$lambda[tuneIndex]
  
  # return model trained on all the received data with tuned lambda parameter
  glmnet(x=LASSO_data[, c(vars, interactions)],
         y=outcomes,
         alpha=alpha,
         lambda=lambda,
         family="gaussian")
}

riskLevel <- function(observed, predicted, clinicalSignificance=2, lowRange=25, upRange=30) {
  levels <- vector(mode="character", length=length(observed))
  
  levels[which(observed>lowRange & observed<upRange & 
                 predicted>lowRange & predicted<upRange)] <- "B"
  
  levels[which(observed>lowRange & observed<upRange & predicted>upRange)] <- "C1"
  
  levels[which(observed>lowRange & observed<upRange & predicted<lowRange)] <- "C2"
  
  levels[which(predicted>lowRange & predicted<upRange & 
                 (observed<lowRange | observed>upRange))] <- "D"
  
  levels[which(observed<lowRange & predicted>upRange)] <- "E1"
  
  levels[which(predicted<lowRange & observed>upRange)] <- "E2"
  
  levels[which(abs(observed-predicted) < clinicalSignificance)] <- "A"
  levels[which(observed<lowRange & predicted<lowRange)] <- "A"
  levels[which(observed>upRange & predicted>upRange)] <- "A"
  
  levels
  
}

validateModelOLS <- function(effect, type, transformation, balancing, new.data) {
  
  stopifnot("Function is for Ridge and LASSO regression only" = type=="OLS")
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  outcome <- "BMI"
  if (transformation=="Inv") {outcome <- "exp(-BMI)"}
  
  # initiate empty vector for collecting prediction error classes
  errorClasses <- c()
  
  # split data in 4 folds
  set.seed(4)
  foldid <- sample(1:4, size=nrow(new.data), replace=TRUE)
  
  # for every fold:
  for (i in 1:4) {
    w <- which(foldid==i)
    # train on other folds
    data_train <- new.data[-w, ]
    n0 <- nrow(data_train)
    if (balancing=="Balanced") {data_train <- oversample(data_train)}
    amount <- nrow(data_train) - n0
    
    intercept_only <- lm(eval(parse(text=paste0(outcome, " ~ 1"))), data=data_train)
    all <- lm(data=data_train,
              formula = eval(parse(text=paste0(outcome, "~`", paste(c("Age", "Smoking", metabolites), 
                                                              collapse="`+`"), "`"))))
    mainModel <- step(intercept_only, direction='both', scope=formula(all), trace=0,
                      k=2*(n0+amount)/n0)
    
    if (effect=="main") {model <- mainModel}
    if (effect=="interaction") {
      variables <- 
        colnames(mainModel$model)[-which(colnames(mainModel$model)==outcome)]
      interactionSet <- c()
      for (i in 1:(length(variables)-1)) {
        for (j in (i+1):length(variables)) {
          term <- paste0(variables[i], "`*`", variables[j])
          addedModel <- lm(data=data_train,
                           formula=eval(parse(text=paste0(outcome, "~`", 
                                                          paste0(c(variables, term), 
                                                                 collapse="`+`"), "`"))))
          AIC0 <- aic(mainModel, oversampled=balancing=="Balanced", amount=amount)
          AIC1 <- aic(addedModel, oversampled=balancing=="Balanced", amount=amount)
          if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
        }
      }
      all <- lm(data=data_train,
                formula = eval(parse(text=paste0(outcome, "~`", paste0(c(variables, interactionSet), 
                                                                 collapse="`+`"), "`"))))
      model <- step(mainModel, direction='both', scope=formula(all), trace=0,
                    k=2*(n0+amount)/n0)
    }
    
    # predict left out fold
    foldPreds <- predict(object=model, newdata=new.data[w, ])
    
    # store error classes
    if (transformation=="Log") {
      obs <- exp(new.data[w, "BMI"])
      preds <- exp(foldPreds)
    }
    if (transformation=="Inv") {
      obs <- exp(new.data[w, "BMI"])
      preds <- 1/foldPreds
    }
    errorClasses <- c(errorClasses, riskLevel(obs, preds))
    
  }
  
  # return table of error classes
  errorClasses <- factor(errorClasses, levels=c("A", "B", "C1", "C2", "D", "E1", "E2"))
  table(errorClasses)
  
}

validateModelRidgeLASSO <- function(effect, type, transformation, balancing, new.data) {
  
  stopifnot("Function is for Ridge and LASSO regression only" = type=="Ridge"|type=="LASSO")
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  includeInteraction <- FALSE
  if (effect=="interaction") {includeInteraction <- TRUE}
  
  outcomes <- new.data$BMI
  if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  
  # initiate empty vector for collecting prediction error classes
  errorClasses <- c()
  
  # split data in 4 folds
  set.seed(4)
  foldid <- sample(1:4, size=nrow(new.data), replace=TRUE)
  
  # for every fold:
  for (i in 1:4) {
    w <- which(foldid==i)
    # train on other folds
    data_train <- new.data[-w,]
    if (balancing=="Balanced") {data_train <- oversample(data_train)}
    RidgeLASSOtrain <- makeMatrix(data_train, includeInteraction)
    if (transformation=="Log") {data_train$transBMI <- data_train$BMI}
    if (transformation=="Inv") {data_train$transBMI <- exp(-data_train$BMI)}
    interactions <- RidgeLASSOtrain$interactions
    RidgeLASSOModel <- glmnet(x=RidgeLASSOtrain$mat,
                              y=data_train$transBMI,
                              alpha=alpha,
                              family="gaussian")
    # predict on other folds to choose the optimal lambda parameter
    RidgeLASSOPreds <- predict(object=RidgeLASSOModel, 
                               newx=RidgeLASSOtrain$mat,
                               type="response")
    
    RidgeLASSOOuts <- matrix(rep(data_train$transBMI, times=RidgeLASSOModel$dim[2]), 
                             byrow=FALSE, ncol=RidgeLASSOModel$dim[2])
    RidgeLASSORes <- RidgeLASSOPreds - RidgeLASSOOuts
    IQR <- sapply(X=1:RidgeLASSOModel$dim[2], 
                  FUN=function(j) {quants <- quantile(RidgeLASSORes[,j], probs=c(0.25, 0.75));
                  quants[2] - quants[1]}
    )
    
    # choose lambda parameter tuned to the smallest inter-quartile range in other folds
    tuneIndex <- which.min(IQR)
    lambda <- RidgeLASSOModel$lambda[tuneIndex]
    
    # train model with chosen lambda parameter all the received data with tuned lambda parameter
    pruneModel <- glmnet(x=RidgeLASSOtrain$mat,
                         y=data_train$transBMI,
                         alpha=alpha,
                         lambda=lambda,
                         family="gaussian")
    
    # predict left out fold
    RidgeLASSOtest <- makeMatrix(new.data[w,], includeInteraction)
    foldPreds <- predict(object=pruneModel, 
                         newx=RidgeLASSOtest$mat,
                         type="response")[,"s0"]
    
    # store error classes
    if (transformation=="Log") {
      obs <- exp(outcomes[w])
      preds <- exp(foldPreds)
    }
    if (transformation=="Inv") {
      obs <- 1/outcomes[w]
      preds <- 1/foldPreds
    }
    errorClasses <- c(errorClasses, riskLevel(obs, preds))
    
  }
  
  # return table of error classes
  errorClasses <- factor(errorClasses, levels=c("A", "B", "C1", "C2", "D", "E1", "E2"))
  return(table(errorClasses))
  
}

tabulateValidation <- function(effects, types, transformations, ethnicity, 
                                balancing=NULL, new.data) {
  if (is.null(ethnicities)) {
    ethnicities <- rep("", times=length(effects))
  } else {
    ethnicities <- rep(ethnicity, times=length(effects))
  }
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Parallelized loop
  validations <- foreach(i=1:length(effects), .combine=bind_rows, .packages=c("dplyr", "glmnet"), 
                         .export=c("validateModelOLS", "validateModelRidgeLASSO", 
                                   "metabolites", "riskLevel", "aic", "oversample", 
                                   "SMOTE", "makeMatrix")) %dopar% {
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    
    # Predict values of the given data with the requested model
    if (type == "OLS") {
      validation <- validateModelOLS(effect, type, transformation, balance, new.data)
    }
    if (type == "Ridge" || type == "LASSO") {
      validation <- validateModelRidgeLASSO(effect, type, transformation, balance, new.data)
    }
    
    return(validation)
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  formatTypes <- types
  formatTypes[which(formatTypes=="Ridge")] <- "ridge"
  formatTransformations <- tolower(transformations)
  formatTransformations[which(formatTransformations=="inv")] <- "1/x"
  
  return(bind_cols("effect"=effects, "type"=formatTypes, "transf."=formatTransformations, 
                   "balancing"=balancing, validations))
}

tabulatePredictionEvaluation <- function(effects, types, transformations, ethnicity,
                                         balancing=NULL, new.data) {
  
  if (is.null(ethnicities)) {ethnicities <- rep("", times=length(effects))}
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (any(effects=="interaction" & (types=="Ridge"|types=="LASSO"))) {
    for (i in 1:(length(vars)-1)) {
      for (j in (i+1):length(vars)) {
        var1 <- vars[i]
        var2 <- vars[j]
        ridge_data <- cbind(ridge_data, ridge_data[,var1]*ridge_data[,var2])
        interactions <- c(interactions, paste(var1, var2, sep="*"))
      }
    }
    colnames(ridge_data) <- c(vars, interactions)
  }
  
  IQR <- vector(mode="numeric", length=length(effects))
  errorTable <- matrix(nrow=length(effects), ncol=6, 0)
  colnames(errorTable) <- c("A", "B", "C1", "C2", "D", "E1", "E2")
  
  for (i in 1:length(effects)) {
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    model <- eval(parse(text=paste0(effect, type, transformation, "Model", 
                                    ethnicity, balance)))
    
    # predict values of the given data with the requested model
    if (type=="OLS") {
      valuePreds <- predict(model, newdata=new.data)
    }
    if (type=="Ridge"|type=="LASSO") {
      if (effect=="main") {
        valuePreds <- 
          predict(model, newx=ridge_data[, c(metabolites, "Age", "Smoking")])[, model$dim[2]]
      }
      if (effect=="interaction") {
        valuePreds <- 
          predict(newx=ridge_data[, c(metabolites, interactions, "Age", "Smoking")],
                  object=model)[, model$dim[2]]
      }
    }
    
    if (transformation=="Log") {predictions <- exp(valuePreds)}
    if (transformation=="Inv") {predictions <- 1/valuePreds}
    
    errorLevels <- riskLevel(observed=exp(new.data$BMI), predicted=predictions)
    w <- which(new.data$ObesityClass== "Overweight")
    QRs <- quantile(exp(new.data$BMI)[w] - predictions[w], probs=c(0.25, 0.75))
    IQR[i] <- QRs[2] - QRs[1]
    errorTable[i, colnames(errorTable)] <- table(errorLevels)[colnames(errorTable)]
    
  }
  errorTable[which(is.na(errorTable))] <- 0
  
  formatTypes <- types
  formatTypes[which(formatTypes=="Ridge")] <- "ridge"
  formatTransformations <- tolower(transformations)
  formatTransformations[which(formatTransformations=="inv")] <- "1/x"
  
  cbind(effect=effects, type=formatTypes, transformation=formatTransformations, 
        balancing, errorTable, "IQR(res)"=round(IQR, digits=3))
  
}

modelDiagnostics <- function(effects, types, transformations, balancing=NULL, 
                             ethnicity, model=NULL, new.data) {
  
  pdf(file=paste0("modelDiagnostics", ethnicity, ".pdf"), width=8.27, height=11.69)
  
  ethnicities <- rep(ethnicity, times=length(effects))
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI", "transBMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  
  init_ridge_data <- makeMatrix(new.data, includeInteraction=TRUE)
  ridge_data <- init_ridge_data$mat
  interactions <- init_ridge_data$interactions
  
  
  for (i in 1:length(effects)) {
    
    layout(matrix(1:2, ncol=1))
    
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    modeltitle <- paste0(effect, type, transformation, "Model", ethnicity, balance)
    
    if (is.null(model)) {
      model <- eval(parse(text=modeltitle))
    }
    
    if (is.null(new.data$predicted)) {
      # predict values of the given data with the requested model
      if (type=="OLS") {
        df_model <- model$df
        valuePreds <- predict(model, newdata=new.data)
      }
      if (type=="Ridge"|type=="LASSO") {
        valuePreds <- 
          predict(newx=ridge_data[, c(vars, interactions)],
                  object=model)[, "s0"]
      }
    } else {
      if (transformation=="Log") {valuePreds <- log(new.data$predicted)}
      if (transformation=="Inv") {valuePreds <- 1/new.data$predicted}
    }
    
    if (transformation=="Log") {valueOuts <- new.data$BMI
                                SSTOT <- sum((new.data$BMI-mean(new.data$BMI))**2)
                                xlabel <- "predicted log(BMI)"}
    if (transformation=="Inv") {valueOuts <- exp(-new.data$BMI)
                                SSTOT <- sum((exp(-new.data$BMI)-mean(exp(-new.data$BMI)))**2)
                                xlabel <- "predicted 1/BMI"}
    
    valueResiduals <- valueOuts - valuePreds
    residuals <- data.frame(prediction=valuePreds, residual=valueResiduals, 
                            residualSquare=(valueResiduals)**2)
    
    
    # plot residuals as function of predicted transformed BMI to assess linearity
    
    plot(x=valuePreds, y=valueResiduals, main=paste0(modeltitle, ": linearity assessment"), 
         xlab=xlabel, ylab="residual")
    abline(h=0, lty="dashed")
    linearity <- loess(formula = residual~prediction, data=residuals)
    linearitySmooth <- predict(linearity, se=TRUE, 
                               newdata=data.frame(prediction=sort(valuePreds)))
    lines(sort(valuePreds), linearitySmooth$fit, col="red")
    lines(sort(valuePreds), linearitySmooth$fit - 1.96*linearitySmooth$se.fit, 
          col="red", lty="dashed")
    lines(sort(valuePreds), linearitySmooth$fit + 1.96*linearitySmooth$se.fit, 
          col="red", lty="dashed")
    
    linearityFit <- lm(formula=residual~1, data=residuals)
    linearitySSE <- sum(linearityFit$residuals**2)
    linearityDOF <- linearityFit$df.residual
    
    smoothSSE <- sum(linearity$residuals**2)
    smoothDOF <- linearity$n - linearity$enp
    
    Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))
    Rsquare <- 1 - linearitySSE/SSTOT
    
    pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)
    text(x=sum(range(valuePreds))/2, y=max(valueResiduals), pos=1, col="red",
         labels=sprintf("linearity: p = %.4f \n R-square = %.2f", pval, Rsquare))
    
    
    # plot squared residuals as function of predicted transformed BMI to assess homoscedasticity
    
    plot(x=range(valuePreds), y=c(1e-5, max((valueResiduals)**2)), 
         col="white", main=paste0(modeltitle, ": homoscedasticity assessment"), 
         xlab=xlabel, ylab="squared residual", log="y")
    points(x=valuePreds, y=(valueResiduals)**2)
    abline(h=mean((valueResiduals)**2), lty="dashed", col="blue")
    homoscedacity <- loess(formula = residualSquare~prediction, data=residuals)
    homoscedacitySmooth <- predict(homoscedacity, se=TRUE, 
                                   newdata=data.frame(prediction=sort(valuePreds)))
    lines(sort(valuePreds), homoscedacitySmooth$fit, col="red")
    lines(sort(valuePreds), homoscedacitySmooth$fit - 1.96*homoscedacitySmooth$se.fit, 
          col="red", lty="dashed")
    lines(sort(valuePreds), homoscedacitySmooth$fit + 1.96*homoscedacitySmooth$se.fit, 
          col="red", lty="dashed")
    
    assess_homoscedasticity <- lm(residualSquare~prediction, data=residuals)
    
    pval <- summary(assess_homoscedasticity)$coefficients["prediction", "Pr(>|t|)"]
    Qpreds <- quantile(valuePreds, probs=c(0.25, 0.75))
    sigma2 <- predict(homoscedacity, Qpreds)
    factorIncrease <- sigma2[2]/sigma2[1]
    
    text(x=sum(range(valuePreds))/2, y=max(valueResiduals**2), pos=1, col="red",
         labels=sprintf("homoscedacity: p = %.4f\nQ1(pred) vs Q3(pred) rel. var. incr. = %.2f", 
                        pval, factorIncrease))
    
  }
  
  dev.off()
  
}

plotPredictions <- function(effects, types, transformations, balancing=NULL, 
                            ethnicity, chosen, new.data) {
  
  pdf(file=paste0("plotPrediction", ethnicity, ".pdf"), width=8.27, height=5.845)
  
  ethnicities <- rep(ethnicity, times=length(effects))
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (any(effects=="interaction" & (types=="Ridge"|types=="LASSO"))) {
    for (i in 1:(length(vars)-1)) {
      for (j in (i+1):length(vars)) {
        var1 <- vars[i]
        var2 <- vars[j]
        ridge_data <- cbind(ridge_data, ridge_data[,var1]*ridge_data[,var2])
        interactions <- c(interactions, paste(var1, var2, sep="*"))
      }
    }
    colnames(ridge_data) <- c(vars, interactions)
  }
  
  riskCs <- list()
  riskEs <- list()
  allPredicted <- c()
  allObserved <- c()
  
  for (i in 1:length(effects)) {
    
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    
    title <- paste0(effect, type, transformation, "Model", ethnicity, balance)
    model <- eval(parse(text=title))
    
    if (i==chosen) {legendlabel <- title}
    
    # predict values of the given data with the requested model
    if (type=="OLS") {
      df_model <- model$df
      preds <- predict(model, newdata=new.data)
    }
    if (type=="Ridge"|type=="LASSO") {
      df_model <- model$dof
      if (effect=="main") {
        preds <- 
          predict(model, newx=ridge_data[, c(metabolites, "Age", "Smoking")])[, model$dim[2]]
      }
      if (effect=="interaction") {
        preds <- 
          predict(newx=ridge_data[, c(metabolites, interactions, "Age", "Smoking")],
                  object=model)[, model$dim[2]]
      }
    }
    
    valueOuts <- exp(new.data$BMI)
    xlabel <- "predicted BMI"
    
    if (transformation=="Log") {valuePreds <- exp(preds)}
    if (transformation=="Inv") {valuePreds <- 1/preds}
    
    allPredicted <- c(allPredicted, valuePreds)
    allObserved <- c(allObserved, valueOuts)
    
    riskCs[[title]] <- riskLevel(valueOuts, valuePreds)%in%c("C1", "C2")
    riskEs[[title]] <- riskLevel(valueOuts, valuePreds)%in%c("E1", "E2")
    
    
  }
  
  # make plot of prediction results
  
  modelColor <- rep("black", times=length(effects))
  modelColor[chosen] <- "red"
  color <- rep(modelColor, each=nrow(new.data))
  
  riskC <- rep(FALSE, times=nrow(new.data))
  riskE <- rep(FALSE, times=nrow(new.data))
  for (title in names(riskCs)) {
    riskC <- riskC | riskCs[[title]]
    riskE <- riskE | riskEs[[title]]
  }
  riskC <- which(riskC)
  riskE <- which(riskE)
  
  subsetAllC <- rep(0:(length(effects)-1), each=length(riskC))*nrow(new.data) + rep(riskC, times=length(effects))
  subsetAllE <- rep(0:(length(effects)-1), each=length(riskE))*nrow(new.data) + rep(riskE, times=length(effects))
  
  plot(x=allPredicted, y=allObserved, xlab="predicted BMI", ylab="observed BMI", 
       main=paste0("validation, ethnicity: ", ethnicity), pch=".", col=color)
  abline(a=-2, b=1, col="green", lty="dashed")
  abline(a=2, b=1, col="green", lty="dashed")
  abline(h=c(25,30), col="green", lty="dashed")
  abline(v=c(25,30), col="green", lty="dashed")
  
  text(x=allPredicted[subsetAllE], y=allObserved[subsetAllE], cex=0.5,
       labels=rep(new.data$ID[riskE], times=length(effects)), col=color[subsetAllE])
  text(x=allPredicted[subsetAllC], y=allObserved[subsetAllC], cex=0.5,
       labels=rep(new.data$ID[riskC], times=length(effects)), col=color[subsetAllC])
  
  legend("topleft", legend=legendlabel, col="red", pch=".")
  
  dev.off()
  
}


