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

trainOLS <- function(filenames) {
  
  if (any(c("mainOLSLogModelWhite.rds", "interactionOLSLogModelWhite.rds")%in%filenames)) {
  
    # main effects, log(BMI), ethnicity="White"
    intercept_only <- lm(BMI ~ 1, data=data_white_train)
    all <- lm(data=data_white_train,
              formula = eval(parse(text=paste0("BMI~`", paste(c("Age", "Smoking", metabolites), 
                                                              collapse="`+`"), "`"))))
    mainOLSLogModelWhite <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSLogModelWhite, "mainOLSLogModelWhite.rds")
    
    # interaction effects, log(BMI), ethnicity="White"
    variables <- 
      colnames(mainOLSLogModelWhite$model)[-which(colnames(mainOLSLogModelWhite$model)=="BMI")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_white_train,
                         formula=eval(parse(text=paste0("BMI~`", 
                                                   paste0(c(variables, term), 
                                                          collapse="`+`"), "`"))))
        if (aic(mainOLSLogModelWhite)>aic(addedModel)) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_white_train,
              formula = eval(parse(text=paste0("BMI~`", paste0(c(variables, interactionSet), 
                                                               collapse="`+`"), "`"))))
    interactionOLSLogModelWhite <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSLogModelWhite, "interactionOLSLogModelWhite.rds")
  
  }
  
  if (any(c("mainOLSInvModelWhite.rds", "interactionOLSInvModelWhite.rds")%in%filenames)) {
    
    # main effects, 1/BMI, ethnicity="White"
    intercept_only <- lm(exp(-BMI) ~ 1, data=data_white_train)
    all <- lm(data=data_white_train,
              formula = eval(parse(text=paste0("exp(-BMI)~`", paste(c("Age", "Smoking", metabolites), 
                                                                    collapse="`+`"), "`"))))
    mainOLSInvModelWhite <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSInvModelWhite, "mainOLSInvModelWhite.rds")
    
    # interaction effects, 1/BMI, ethnicity="White"
    variables <- 
      colnames(mainOLSInvModelWhite$model)[-which(colnames(mainOLSInvModelWhite$model)=="exp(-BMI)")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_white_train,
                         formula=eval(parse(text=paste0("exp(-BMI)~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        if (aic(mainOLSInvModelWhite)>aic(addedModel)) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_white_train,
              formula = eval(parse(text=paste0("exp(-BMI)~`", paste0(c(variables, interactionSet), 
                                                               collapse="`+`"), "`"))))
    interactionOLSInvModelWhite <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSInvModelWhite, "interactionOLSInvModelWhite.rds")
    
  }
  
  
  
  if (any(c("mainOLSLogModelWhiteBalanced.rds", "interactionOLSLogModelWhiteBalanced.rds")%in%filenames)) {
  
    # main effects, log(BMI), ethnicity="White", balanced
    intercept_only <- lm(BMI ~ 1, data=data_white_train_balanced)
    all <- lm(data=data_white_train_balanced,
              formula = eval(parse(text=paste0("BMI~`", paste(c("Age", "Smoking", metabolites), 
                                                              collapse="`+`"), "`"))))
    mainOLSLogModelWhiteBalanced <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSLogModelWhiteBalanced, "mainOLSLogModelWhiteBalanced.rds")
    
    # interaction effects, log(BMI), ethnicity="White", balanced
    variables <- 
      colnames(mainOLSLogModelWhiteBalanced$model)[-which(colnames(mainOLSLogModelWhiteBalanced$model)=="BMI")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_white_train_balanced,
                         formula=eval(parse(text=paste0("BMI~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        AIC0 <- aic(mainOLSLogModelWhiteBalanced, oversampled=TRUE, 
                    amount=nrow(data_white_train_balanced)-nrow(data_white_train))
        AIC1 <- aic(addedModel, oversampled=TRUE, 
                    amount=nrow(data_white_train_balanced)-nrow(data_white_train))
        if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_white_train_balanced,
              formula = eval(parse(text=paste0("BMI~`", paste0(c(variables, interactionSet), 
                                                               collapse="`+`"), "`"))))
    interactionOLSLogModelWhiteBalanced <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSLogModelWhiteBalanced, "interactionOLSLogModelWhiteBalanced.rds")
    
  }
  
  if (any(c("mainOLSInvModelWhiteBalanced.rds", "interactionOLSInvModelWhiteBalanced.rds")%in%filenames)) {
    
    # main effects, 1/BMI, ethnicity="White", balanced
    intercept_only <- lm(exp(-BMI) ~ 1, data=data_white_train_balanced)
    all <- 
      lm(data=data_white_train_balanced,
         formula=eval(parse(text=paste0("exp(-BMI)~`", paste(c("Age", "Smoking", metabolites), 
                                                             collapse="`+`"), "`"))))
    mainOLSInvModelWhiteBalanced <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSInvModelWhiteBalanced, "mainOLSInvModelWhiteBalanced.rds")
    
    # interaction effects, 1/BMI, ethnicity="White", balanced
    variables <- 
      colnames(mainOLSInvModelWhiteBalanced$model)[-which(colnames(mainOLSInvModelWhiteBalanced$model)=="exp(-BMI)")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_white_train_balanced,
                         formula=eval(parse(text=paste0("exp(-BMI)~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        AIC0 <- aic(mainOLSInvModelWhiteBalanced, oversampled=TRUE, 
                    amount=nrow(data_white_train_balanced)-nrow(data_white_train))
        AIC1 <- aic(addedModel, oversampled=TRUE, 
                    amount=nrow(data_white_train_balanced)-nrow(data_white_train))
        if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_white_train_balanced,
              formula = eval(parse(text=paste0("exp(-BMI)~`", paste0(c(variables, interactionSet), 
                                                                     collapse="`+`"), "`"))))
    interactionOLSInvModelWhiteBalanced <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSInvModelWhiteBalanced, "interactionOLSInvModelWhiteBalanced.rds")
    
  }
  
  
  
  if (any(c("mainOLSLogModelBlack.rds", "interactionOLSLogModelBlack.rds")%in%filenames)) {
    
    # main effects, log(BMI), ethnicity="Black"
    intercept_only <- lm(BMI ~ 1, data=data_black_train)
    all <- lm(data=data_black_train,
              formula = eval(parse(text=paste0("BMI~`", paste(c("Age", "Smoking", metabolites), 
                                                              collapse="`+`"), "`"))))
    mainOLSLogModelBlack <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSLogModelBlack, "mainOLSLogModelBlack.rds")
    
    # interaction effects, log(BMI), ethnicity="Black"
    variables <- 
      colnames(mainOLSLogModelBlack$model)[-which(colnames(mainOLSLogModelBlack$model)=="BMI")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_black_train,
                         formula=eval(parse(text=paste0("BMI~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        if (aic(mainOLSLogModelBlack)>aic(addedModel)) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_black_train,
              formula = eval(parse(text=paste0("BMI~`", paste0(c(variables, interactionSet), 
                                                               collapse="`+`"), "`"))))
    interactionOLSLogModelBlack <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSLogModelBlack, "interactionOLSLogModelBlack.rds")
    
  }
  
  if (any(c("mainOLSInvModelBlack.rds", "interactionOLSInvModelBlack.rds")%in%filenames)) {
  
    # main effects, 1/BMI, ethnicity="Black"
    intercept_only <- lm(exp(-BMI) ~ 1, data=data_black_train)
    all <- lm(data=data_black_train,
              formula=eval(parse(text=paste0("exp(-BMI)~`", paste(c("Age", "Smoking", metabolites), 
                                                                    collapse="`+`"), "`"))))
    mainOLSInvModelBlack <- step(intercept_only, direction='both', scope=formula(all), trace=0)
    saveRDS(mainOLSInvModelBlack, "mainOLSInvModelBlack.rds")
    
    # interaction effects, 1/BMI, ethnicity="Black"
    variables <- 
      colnames(mainOLSInvModelBlack$model)[-which(colnames(mainOLSInvModelBlack$model)=="exp(-BMI)")]
    interactionSet <- c()
    for (i in 1:(length(variables)-1)) {
      for (j in (i+1):length(variables)) {
        term <- paste0(variables[i], "`*`", variables[j])
        addedModel <- lm(data=data_black_train,
                         formula=eval(parse(text=paste0("exp(-BMI)~`", 
                                                        paste0(c(variables, term), 
                                                               collapse="`+`"), "`"))))
        if (aic(mainOLSInvModelBlack)>aic(addedModel)) {interactionSet <- c(interactionSet, term)}
      }
    }
    all <- lm(data=data_black_train,
              formula=eval(parse(text=paste0("exp(-BMI)~`", paste0(c(variables, interactionSet), 
                                                                   collapse="`+`"), "`"))))
    interactionOLSInvModelBlack <- step(all, direction='both', scope=formula(all), trace=0)
    saveRDS(interactionOLSInvModelBlack, "interactionOLSInvModelBlack.rds")
  
  }
  
}

trainRidge <- function(new.data, transformation="Log", interactionEffect=FALSE) {
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ObesityClass", "BMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (interactionEffect) {
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
  
  outcomes <- new.data$BMI
  if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  
  # divide ridge_data in 4 folds for cross-validation
  set.seed(4)
  foldid <- sample(1:4, size=nrow(ridge_data), replace=TRUE)
  
  # stabilize probed lambda values with first validation + initiate IQRs
  w <- which(foldid==1)
  ridgeModel <- glmnet(x=ridge_data[-w, c(metabolites, interactions, "Age", "Smoking")],
                       y=outcomes[-w],
                       alpha=0,
                       lambda.min.ratio=0.000001,
                       family="gaussian")
  ridgePreds <- predict(object=ridgeModel, 
                        newx=ridge_data[w, c(metabolites, interactions, "Age", "Smoking")],
                        type="response")
  ridgeOuts <- matrix(rep(outcomes[w], times=ridgeModel$dim[2]), 
                      byrow=FALSE, ncol=ridgeModel$dim[2])
  ridgeRes <- ridgePreds - ridgeOuts
  IQRs <- sapply(X=1:ridgeModel$dim[2], 
                 FUN=function(j) {quants <- quantile(ridgeRes[,j], probs=c(0.25, 0.75));
                 quants[2] - quants[1]}
  )
  lambdas <- ridgeModel$lambda
  
  # iterate over other cross folds
  for (i in 2:4) {
    w <- which(foldid==i)
    ridgeModel <- glmnet(x=ridge_data[-w, c(metabolites, interactions, "Age", "Smoking")],
                         y=outcomes[-w],
                         alpha=0,
                         lambda=lambdas,
                         family="gaussian")
    ridgePreds <- predict(object=ridgeModel, 
                          newx=ridge_data[w, c(metabolites, interactions, "Age", "Smoking")],
                          type="response")
    ridgeOuts <- matrix(rep(outcomes[w], times=ridgeModel$dim[2]), 
                        byrow=FALSE, ncol=ridgeModel$dim[2])
    ridgeRes <- ridgePreds - ridgeOuts
    IQRi <- sapply(X=1:ridgeModel$dim[2], 
                   FUN=function(j) {quants <- quantile(ridgeRes[,j], probs=c(0.25, 0.75));
                   quants[2] - quants[1]}
    )
    IQRs <- IQRs + IQRi
  }
  IQRs <- IQRs/4
  
  # choose lambda parameter tuned to the smallest inter-quartile range in validation set
  tuneIndex <- which.min(IQRs)
  lambda <- ridgeModel$lambda[tuneIndex]
  
  # return model trained on all the received data with tuned lambda parameter
  glmnet(x=ridge_data[, c(metabolites, interactions, "Age", "Smoking")],
         y=outcomes,
         alpha=0,
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
  
  levels[which((predicted<lowRange & observed>upRange) | 
                 (observed<lowRange & predicted>upRange))] <- "E"
  
  levels[which(abs(observed-predicted) < clinicalSignificance)] <- "A"
  levels[which(observed<lowRange & predicted<lowRange)] <- "A"
  levels[which(observed>upRange & predicted>upRange)] <- "A"
  
  levels
  
}

tabulatePredictionEvaluation <- function(effects, types, transformations, ethnicity,
                                         balancing=NULL, new.data) {
  
  if (is.null(ethnicities)) {ethnicities <- rep("", times=length(effects))}
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ObesityClass", "BMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (any(effects=="interaction" & types=="Ridge")) {
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
  colnames(errorTable) <- c("A", "B", "C1", "C2", "D", "E")
  
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
    if (type=="Ridge") {
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
    QRs <- quantile(exp(new.data$BMI) - predictions, probs=c(0.25, 0.75))
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
                             ethnicity, new.data) {
  
  pdf(file=paste0("modelDiagnostics", ethnicity, ".pdf"), width=8.27, height=11.69)
  
  ethnicities <- rep(ethnicity, times=length(effects))
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ObesityClass", "BMI"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars <- colnames(new.data)[-irrelevant]
  metabolites <- vars[-which(vars %in% c("Smoking", "Age"))]
  interactions <- c()
  if (any(effects=="interaction" & types=="Ridge")) {
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
  
  
  for (i in 1:length(effects)) {
    
    layout(matrix(1:2, ncol=1))
    
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    
    title <- paste0(effect, type, transformation, "Model", ethnicity, balance)
    model <- eval(parse(text=title))
    
    # predict values of the given data with the requested model
    if (type=="OLS") {
      df_model <- model$df
      valuePreds <- predict(model, newdata=new.data)
    }
    if (type=="Ridge") {
      df_model <- model$dof
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
    
    plot(x=valuePreds, y=valueResiduals, main=paste0(title, ": linearity assessment"), 
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
         col="white", main=paste0(title, ": homoscedasticity assessment"), 
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




