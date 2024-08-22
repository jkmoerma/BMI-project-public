#' Computes the Akaike information criterion of an ordinary least square regression model. 
#' The AIC was corrected for synthetically oversampled data, for which the sample size is unrepresentable
#'
#' @param model An ordinary least square regression model
#' @param oversampled A logical value indicating whether the data set used to train the regression model was oversampled or not.
#' @param amount The number of synthetically generated observations added to the oversampled data set
#' @return A numeric AIC value calculated as $n log(SSE) - n log(n) + 2p$
#' @examples
#' 
aic <- function(model, oversampled=FALSE, amount=0) {
  n <- length(model$coefficients)+model$df.residual
  SSE <- sum(model$residuals**2)
  p <- length(model$coefficients)
  
  if (oversampled) {
    n <- n-amount   # substract the amount of synthetical observations from the sample size
    SSE <- SSE*n/(n+amount)   # rescale the SSE to the original sample size
  }
  
  n*log(SSE) - n*log(n) + 2*p
}

#' Oversamples a given data set to balance out the "Overweight" and "Obese" obesity classes with respect to the "Normal weight" using SMOTE.
#' The data set must contain the two columns called "ObesityClass" and "met_110",  and some columns of the globally specified "metabolites", "BMI" and/or "Age".
#' 
#' @param data_train A data frame. One of the colums must be called "ObesityClass", consisting of values "Normal weight", "Overweight" and "Obese". Some columns must be 
#' @return A SMOTE oversampled version of the data set as a data frame. Patient ID and other irrelevant columns are not included
#' @examples
#' library(DMwR)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' 
#' metabolites <- c("met1", "met2", "met_110")
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=20),
#'                           Age=rnorm(n=20, mean=30, sd=5),
#'                           BMI=runif(n=20, min=18, max=25),
#'                           met1=rnorm(n=20, mean=0.5, sd=0.05),
#'                           met2=rnorm(n=20, mean=5, sd=1),
#'                           met_110=rnorm(n=20, mean=100, sd=30))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=10),
#'                           Age=rnorm(n=10, mean=30, sd=5),
#'                           BMI=runif(n=10, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.05),
#'                           met2=rnorm(n=10, mean=4, sd=1),
#'                           met_110=rnorm(n=10, mean=100, sd=30))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=5),
#'                           Age=rnorm(n=5, mean=30, sd=5),
#'                           BMI=runif(n=5, min=30, max=40),
#'                           met1=rnorm(n=5, mean=1, sd=0.05),
#'                           met2=rnorm(n=5, mean=3, sd=1),
#'                           met_110=rnorm(n=5, mean=100, sd=30))
#'                           
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' oversampled_data <- oversample(data_train)
#' 
#' table(oversampled_data$ObesityClass)
#' boxplot(BMI~ObesityClass, oversampled_data)
#' boxplot(BMI~ObesityClass, data_train)
#' 
oversample <- function(data_train) {
  
  stopifnot("met_110 is used for SMOTE interpolation and is not included" = "met_110"%in%names(data_train))
  
  irrelevant <- which(! names(data_train) %in% c(metabolites, "BMI", "Age", "ObesityClass"))
  
  for (name in names(data_train)[irrelevant]) {
    data_train[[name]] <- NULL
  }
  
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
  
  # recalculate metabolite ratios ON LOG SCALE !!
  ratios <- metabolites[which(grepl(pattern="/", metabolites))]
  for (ratio in ratios) {
    mets <- strsplit(ratio, split="/")[[1]]
    data_balanced[[ratio]] <- with(data_balanced, eval(parse(text=paste(mets[1], "-", mets[2]))))
  }
  
  data_balanced
  
}

#' This function makes a matrix of a data frame with metabolite measurements and Age. It is usually used for training and predicting with the glmnet-function.
#' Columns for interaction terms can be included if requested
#' 
#' @param data A data frame of metabolite values and "Age".
#' @param includeInteraction Logical value for optionally including interaction terms if requested
#' @return A list containing a matrix (mat) and the interaction terms included in the matrtix (interactions)
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' 
#' metabolites <- c("met1", "met2", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=20),
#'                           Age=rnorm(n=20, mean=30, sd=5),
#'                           BMI=runif(n=20, min=18, max=25),
#'                           met1=rnorm(n=20, mean=0.5, sd=0.05),
#'                           met2=rnorm(n=20, mean=5, sd=1),
#'                           met3=rnorm(n=20, mean=100, sd=30))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=10),
#'                           Age=rnorm(n=10, mean=30, sd=5),
#'                           BMI=runif(n=10, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.05),
#'                           met2=rnorm(n=10, mean=4, sd=1),
#'                           met3=rnorm(n=10, mean=100, sd=30))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=5),
#'                           Age=rnorm(n=5, mean=30, sd=5),
#'                           BMI=runif(n=5, min=30, max=40),
#'                           met1=rnorm(n=5, mean=1, sd=0.05),
#'                           met2=rnorm(n=5, mean=3, sd=1),
#'                           met3=rnorm(n=5, mean=100, sd=30))
#'                           
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' 
#' makeMatrix(data_train, includeInteraction = TRUE)
#' 
makeMatrix <- function(data, includeInteraction=FALSE) {
  relevant <- which(colnames(data) %in% c("Age", metabolites))
  RidgeLASSO_data <- as.matrix(data[, relevant])
  vars0 <- colnames(data)[relevant]
  interactions <- c()
  
  if (includeInteraction) {
    n <- ncol(RidgeLASSO_data)
    num_interactions <- n * (n - 1) / 2
    
    # Pre-allocate a matrix for interactions
    interaction_matrix <- matrix(NA, nrow = nrow(RidgeLASSO_data), ncol = num_interactions)
    
    interaction_index <- 1
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        interaction_matrix[, interaction_index] <- RidgeLASSO_data[,i] * RidgeLASSO_data[,j]
        interactions <- c(interactions, paste(vars0[i], vars0[j], sep="*"))
        interaction_index <- interaction_index + 1
      }
    }
    
    # Combine original data with interactions
    RidgeLASSO_data <- cbind(RidgeLASSO_data, interaction_matrix)
    colnames(RidgeLASSO_data) <- c(vars0, interactions)
  }
  
  return(list(mat = RidgeLASSO_data, interactions = interactions))
}

#' Train a stepwise OLS model from the data with a specification about main/interaction effects and a BMI transformation.
#' 
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "OLS" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An OLS regression model
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' 
#' metabolites <- c("met1", "met2", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=20),
#'                           Age=rnorm(n=20, mean=30, sd=5),
#'                           BMI=runif(n=20, min=18, max=25),
#'                           met1=rnorm(n=20, mean=0.5, sd=0.05),
#'                           met2=rnorm(n=20, mean=5, sd=1),
#'                           met_110=rnorm(n=20, mean=100, sd=30))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=10),
#'                           Age=rnorm(n=10, mean=30, sd=5),
#'                           BMI=runif(n=10, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.05),
#'                           met2=rnorm(n=10, mean=4, sd=1),
#'                           met_110=rnorm(n=10, mean=100, sd=30))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=5),
#'                           Age=rnorm(n=5, mean=30, sd=5),
#'                           BMI=runif(n=5, min=30, max=40),
#'                           met1=rnorm(n=5, mean=1, sd=0.05),
#'                           met2=rnorm(n=5, mean=3, sd=1),
#'                           met_110=rnorm(n=5, mean=100, sd=30))
#'                           
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' data_train$BMI <- log(data_train$BMI)
#' trainOLS(effect="main", type="OLS", transformation="Log", balancing="balanced", new.data=data_train)
#' 
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
            formula = eval(parse(text=paste0("transBMI~`", paste(c("Age", metabolites), 
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

#' Hyper parameter tuning of a penalty lambda parameter for LASSO and ridge regression. A 4-fold cross validation is performed with the received data.
#' For every fold, a model is trained on the 3 other folds with 100 probed lambda parameters to predict the values of the left-out fold.
#' Per fold, the AUC is calculated for normal weight versus obese patients for every probed lambda.
#' The AUC values are averaged over the 4 folds for each probed lambda.
#' The lambda with the highest AUC fold average is chosen
#' 
#' @param data_matrix Measured metabolite values and possibly in matrix format
#' @param alpha A parameter passed to glmnet. Set to 0 for ridge regression and to 1 for LASSO regression
#' @param new.y A vector having the same length as the number of rows of data_matrix. It represents the log(BMI) values according the data in data_matrix
#' @param effect "main" or "interaction" dependent on which effects are wished in the model
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param returnAll Logical. If FALSE, return the tuned lambda parameter, if TRUE return all probed lambdas and probed AUCs as a data frame
#' @return 
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(glmnet)
#' library(pROC)
#' 
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=5, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=4, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                           
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' data_train$BMI <- log(data_train$BMI)
#' data_train$met_coll <- data_train$met1 + 7.2*data_train$met2
#' 
#' data_train <- oversample(data_train)
#' AUCs <- tuneLambda(data_matrix=makeMatrix(data_train)$mat, 
#'                    new.y = data_train$BMI, alpha=0,
#'                    effect="main", transformation="Log",
#'                    returnAll=TRUE)
#' bestLambda <- AUCs$lambda[which.max(AUCs$AUC)]
#' plot(AUC~log(lambda), AUCs)
#' abline(v=log(bestLambda), col="blue")
#' 
tuneLambda <- function(data_matrix, alpha, new.y, effect,
                          transformation, returnAll=FALSE) {
  
  # divide data_matrix in 10 folds for cross-validation
  set.seed(15)
  tuneSelect <- sample(x=1:nrow(data_matrix), size=round(0.25*nrow(data_matrix)))
  
  # define outcomes for regression
  if (transformation=="Log") {outcomes <- new.y}
  if (transformation=="Inv") {outcomes <- exp(-new.y)}
  
  ## calculate a start set of lambda values
  
  # the highest lambda parameter meaningful for the model
  model0 <- glmnet(x=data_matrix[-tuneSelect, ],
                   y=outcomes[-tuneSelect],
                   alpha=alpha,
                   nlambda=3,
                   lambda.min.ratio=0.1,
                   family="gaussian")
  max.lambda <- max(model0$lambda)
  
  # a lower limit lambda parameter defined by where the percentage deviance explained is almost 1
  min.lambda.series <- max.lambda * 10**(-(0:10))
  model1 <- glmnet(x=data_matrix[-tuneSelect, ],
                   y=outcomes[-tuneSelect],
                   alpha=alpha,
                   lambda=min.lambda.series,
                   family="gaussian")
  minIndex <- min(which(1-model1$dev.ratio/max(model1$dev.ratio)<1e-2)) + 2
  min.lambda <- model1$lambda[minIndex]
  
  # a sequence of desired lambda values from the model
  lambda_sequence <- exp(seq(from=log(max.lambda), to=log(min.lambda), length=100))
  
  # define folds for lambda cross validation
  set.seed(7)
  foldid <- sample(1:4, size=nrow(data_matrix), replace=TRUE) 
  
  AUCs <- 0
  for (i in 1:4) {
    w <- which(foldid==i)
    
    # train models with probed lambda values
    model <- glmnet(x=data_matrix[-w, ],
                    y=outcomes[-w],
                    alpha=alpha,
                    lambda=lambda_sequence,
                    family="gaussian")
    preds <- predict(object=model, 
                     newx=data_matrix[w, ],
                     type="response")
    BMIs <- exp(new.y)[w]
    
    # calculate AUC from predictions of the left-out data at the different lambda values probed
    AUCs_i <- sapply(X=1:model$dim[2], 
                   FUN=function(j) {
                     w_normal <- which(BMIs < 25)
                     w_obese <- which(BMIs > 30)
                     roc_j <- roc(controls=preds[w_normal,j], 
                                  cases=preds[w_obese,j])
                     auc(roc_j)
                   }
    )
    AUCs <- AUCs + AUCs_i
  }
  AUCs <- AUCs/4
 
  
  # return lambda (returnAll=FALSE) or whole evaluation (returnAll=TRUE)
  if (returnAll) {
    return(data.frame(lambda=model$lambda, AUC=AUCs))
  } else {
    # choose lambda parameter tuned to the highest cross-validated AUC
    tuneIndex <- which.max(AUCs)
    lambda <- model$lambda[tuneIndex]
    return(lambda)
  }
  
}

#' Train a Ridge or LASSO model from the data with a specification about main/interaction effects and a BMI transformation.
#' 
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return A glmnet regression model with 1 lambda parameter
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(glmnet)
#' library(pROC)
#' 
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=5, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=4, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                         
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' data_train$BMI <- log(data_train$BMI)
#' data_train$met_coll <- data_train$met1 + 7.2*data_train$met2 + rnorm(n=nrow(data_train), mean=0, sd=0.1)
#' 
#' data_train <- oversample(data_train)
#' trainRidgeLASSO(effect="main", type="Ridge", transformation="Log", new.data=data_train)
#' 
trainRidgeLASSO <- function(effect, type, transformation,
                            new.data=NULL, new.x=NULL, new.y=NULL, lambda=NULL) {
  
  if (!is.null(new.x)) {
    stopifnot("Outcomes must be provided separately if data matrix is used" = 
                !is.null(new.y))
  }
  
  interactionEffect <- FALSE
  if (effect=="interaction") {interactionEffect <- TRUE}
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  if (is.null(new.x)) {
    generateX <- makeMatrix(new.data, includeInteraction = effect=="interaction")
    regression_data <- generateX$mat
    interactions <- generateX$interactions
    new.y <- new.data$BMI
    outcomes <- new.data$BMI
    if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  } else {
    regression_data <- new.x
    outcomes <- new.y
    if (transformation=="Inv") {outcomes <- exp(-new.y)}
  }
  
  # if not provided, tune a lambda parameter for the final model
  if (is.null(lambda)) {
    lambda <- tuneLambda(regression_data, alpha, new.y, effect,
                         transformation, returnAll=FALSE)
  }
  
  # return model trained on all the received data with tuned lambda parameter
  glmnet(x=regression_data,
         y=outcomes,
         alpha=alpha,
         lambda=lambda,
         family="gaussian")
}

#' Use 4-fold cross validation inside new.data to score a regression formulation with type OLS. 
#' Calculate the AUC of every fold building a model from the other 3 folds.
#' Return the average AUC over the 4 folds.
#' 
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "OLS" for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#'
validateModelOLS <- function(effect, type, transformation, balancing, new.data) {
  
  stopifnot("Function is for Ridge and LASSO regression only" = type=="OLS")
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  outcome <- "BMI"
  if (transformation=="Inv") {outcome <- "exp(-BMI)"}
  
  # initiate data frame collecting data + predictions
  AUC <- 0
  
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
              formula = eval(parse(text=paste0(outcome, "~`", paste(c("Age",  metabolites), 
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
    w_normal <- which(obs < 25)
    w_obese <- which(obs > 30)
    roc_j <- roc(controls=preds[w_normal], cases=preds[w_obese])
    
    AUC <- AUC + auc(roc_j)
  }
  
  return(round(AUC/4, digits=3))
  
}

#' Use 4-fold cross validation inside new.data to score a regression formulation with type Ridge or LASSO. 
#' Calculate the AUC of every fold building a model from the other 3 folds.
#' Return the average AUC over the 4 folds.
#' 
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#'
validateModelRidgeLASSO <- function(effect, type, transformation, balancing, 
                                    new.data) {
  
  stopifnot("Function is for Ridge and LASSO regression only" = type=="Ridge"|type=="LASSO")
  
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  includeInteraction <- FALSE
  if (effect=="interaction") {includeInteraction <- TRUE}
  
  outcomes <- new.data$BMI
  if (transformation=="Inv") {outcomes <- exp(-new.data$BMI)}
  
  # initiate pooled AUC
  AUC <- 0
  
  # split data in 4 folds
  set.seed(4)
  foldid <- sample(1:4, size=nrow(new.data), replace=TRUE)
  
  # for every fold:
  for (i in 1:4) {
    w <- which(foldid==i)
    
    # train on other folds (no lambda specified, therefore tuned)
    pruneModel <- trainRidgeLASSO(effect, type, transformation, 
                                  new.data=new.data[-w,])
    
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
    w_normal <- which(obs < 25)
    w_obese <- which(obs > 30)
    roc_j <- roc(controls=preds[w_normal], cases=preds[w_obese])
    
    AUC <- AUC + auc(roc_j)
    
  }
  
  return(round(AUC/4, digits=3))
  
}

#' Use 4-fold cross validation inside new.data to score the provided regression formulations 
#' For every regression formulation, calculate the AUC of every fold building a model from the other 3 folds.
#' For every regression formulation, return the average AUC over the 4 folds.
#' 
#' @param effects a vector of "main" and/or "interaction" effects of every formulation
#' @param types a vector of regression types "OLS", "Ridge" and/or "LASSO" of every formulation
#' @param transformation a vector of BMI transformations "Log" and/or "Inv" of every formulation
#' @param balancing defaults to NULL, or a vector of "" and/or "Balanced" reflecting balancing requests of every formulation
#' @param new.data Data to build a regression model from
#' @return An averaged AUC value representing a score of the model formulation
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(glmnet)
#' library(pROC)
#' library(doParallel)
#' 
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=5, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=4, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                         
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' data_train$BMI <- log(data_train$BMI)
#' data_train$met_coll <- data_train$met1 + 7.2*data_train$met2 + rnorm(n=nrow(data_train), mean=0, sd=0.1)
#' 
#' formulations <- expand.grid(effects=c("main", "interaction"), 
#'                             types=c("OLS", "Ridge", "LASSO"), 
#'                             transformations=c("Log", "Inv"),
#'                             balancing=c("", "Balanced"),
#'                             stringsAsFactors=FALSE)
#' 
#' formulationScores <- tabulateValidation(formulations$effects, formulations$types, 
#'                                         formulations$transformations,
#'                                         balancing=formulations$balancing, data_train)
#' 
#' View(formulationScores)
#' 
tabulateValidation <- function(effects, types, transformations, 
                                balancing=NULL, new.data) {
  
  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Parallelized loop
  validations <- foreach(i=1:length(effects), .combine=c, 
                         .packages=c("dplyr", "glmnet"), 
                         .export=c("validateModelOLS", "validateModelRidgeLASSO", 
                                   "metabolites", "aic", "oversample", 
                                   "SMOTE", "makeMatrix", "roc", "auc", 
                                   "trainRidgeLASSO", "tuneLambda")) %dopar% {
    effect <- effects[i]
    type <- types[i]
    transformation <- transformations[i]
    balance <- balancing[i]
    
    # Predict values of the given data with the requested model
    if (type == "OLS") {
      validation <- validateModelOLS(effect, type, transformation, balance, new.data)
    }
    if (type == "Ridge" || type == "LASSO") {
      validation <- validateModelRidgeLASSO(effect, type, transformation, balance,
                                            new.data)
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
                   "balancing"=balancing, "AUC.cv"=validations))
}

#' Stores a pdf file called paste0("modelDiagnostics", ethnicity, ".pdf") with plots of the residuals and squared residuals as a function of the predicted BMI.
#' The functional form of the predictors significantly deviating from linearity are included.
#' All p-values are calculated from F-tests between an loess-estimator and an intercept estimator of the residuals
#' 
#' @param effects a vector of "main" and/or "interaction" effects of formulations the model diagnostics need to be examined of
#' @param types a vector of regression types "OLS", "Ridge" and/or "LASSO" of formulations the model diagnostics need to be examined of
#' @param transformation a vector of BMI transformations "Log" and/or "Inv" of formulations the model diagnostics need to be examined of
#' @param balancing defaults to NULL, or a vector of "" and/or "Balanced" reflecting balancing requests of formulations the model diagnostics need to be examined of
#' @param ethnicity the ethnicity of patients in new.data
#' @param new.data Data to build a regression model from
#' @return 
#' 
modelDiagnostics <- function(effects, types, transformations, balancing=NULL, 
                             ethnicity, model=NULL, new.data) {
  
  pdf(file=paste0("modelDiagnostics", ethnicity, ".pdf"), width=8.27, height=11.69)
  
  ethnicities <- rep(ethnicity, times=length(effects))
  if (is.null(balancing)) {balancing <- rep("", times=length(effects))}
  
  irrelevant <- which(colnames(new.data)%in%c("Race", "ID", "ObesityClass", "BMI", "transBMI", "Smoking"))
  ridge_data <- as.matrix(new.data[,-irrelevant])
  vars0 <- colnames(new.data)[-irrelevant]
  
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
        if (effect=="main") {
          valuePreds <- 
            predict(newx=ridge_data[, c(vars0)],
                    object=model)[, "s0"]
        } else {
          valuePreds <- 
            predict(newx=ridge_data[, c(vars0, interactions)],
                    object=model)[, "s0"]
        }
        
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
    
    
    # make second page plotting the covariates with a non-linear functional form
    
    # checking the functional form of the variables entering the model
    nonlinears <- c()
    for (met in c("Age", metabolites)) {
      residuals <- data.frame(feature=new.data[[met]], residual=valueResiduals, 
                              residualSquare=(valueResiduals)**2)
      
      linearity <- loess(formula = residual~feature, data=residuals)
      linearitySmooth <- predict(linearity, se=TRUE, 
                                 newdata=data.frame(feature=sort(new.data[[met]])))
      
      linearityFit <- lm(formula=residual~feature, data=residuals)
      linearitySSE <- sum(linearityFit$residuals**2)
      linearityDOF <- linearityFit$df.residual
      
      smoothSSE <- sum(linearity$residuals**2)
      smoothDOF <- linearity$n - linearity$enp
      
      Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))
      
      pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)
      
      if (pval<0.05) {nonlinears <- c(nonlinears, met)}
      
    }
    
    if (length(nonlinears) > 0) {
      
      layout(matrix(1:6, ncol=2, byrow=TRUE))
      
      variables <- c("Age", metabolites)
      
      reduced_vars <- variables[-which(variables %in% nonlinears)]
      reduced_model <- lm(formula = eval(parse(text=paste0("BMI~`", paste(reduced_vars, collapse="`+`"), "`"))),
                          data=new.data)
      for (met in nonlinears) {
        
        residuals <- data.frame(feature=new.data[[met]], residual=valueResiduals, 
                                residualSquare=(valueResiduals)**2)
        
        plot(x=new.data[[met]], y=valueResiduals, main=paste0("model residuals vs. ", met), 
             xlab=met, ylab="residual")
        abline(h=0, lty="dashed")
        linearity <- loess(formula = residual~feature, data=residuals)
        linearitySmooth <- predict(linearity, se=TRUE, 
                                   newdata=data.frame(feature=sort(new.data[[met]])))
        lines(sort(new.data[[met]]), linearitySmooth$fit, col="red")
        lines(sort(new.data[[met]]), linearitySmooth$fit - 1.96*linearitySmooth$se.fit, 
              col="red", lty="dashed")
        lines(sort(new.data[[met]]), linearitySmooth$fit + 1.96*linearitySmooth$se.fit, 
              col="red", lty="dashed")
        
        linearityFit <- lm(formula=residual~feature, data=residuals)
        linearitySSE <- sum(linearityFit$residuals**2)
        linearityDOF <- linearityFit$df.residual
        
        smoothSSE <- sum(linearity$residuals**2)
        smoothDOF <- linearity$n - linearity$enp
        
        Fstat <- (linearitySSE-smoothSSE)*smoothDOF/(smoothSSE*(linearityDOF-smoothDOF))
        
        pval <- 1-pf(q=Fstat, df1=linearityDOF-smoothDOF, df2=smoothDOF)
        text(x=sum(range(new.data[[met]]))/2, y=max(valueResiduals), pos=1,
             labels=sprintf("linearity: p = %.2e", pval), col="red")
        
        modelVar <- lm(formula = eval(parse(text=paste(met,"~",paste(reduced_vars, collapse="+")))),
                       data=new.data)
        
        plot(x=modelVar$residuals, y=reduced_model$residuals, 
             xlab=paste0(met, "-pred(", met, ")"), ylab="residual", 
             main=paste0("Desired functional form ", met))
        abline(h=0, lty="dashed")
        residuals <- data.frame(feature=modelVar$residuals, residual=reduced_model$residuals)
        linearity <- loess(formula = residual~feature, data=residuals)
        linearitySmooth <- predict(linearity, se=TRUE,
                                   newdata=data.frame(feature=sort(modelVar$residuals)))
        lines(sort(modelVar$residuals), linearitySmooth$fit, col="red")
        lines(sort(modelVar$residuals), linearitySmooth$fit - 1.96*linearitySmooth$se.fit, 
              col="red", lty="dashed")
        lines(sort(modelVar$residuals), linearitySmooth$fit + 1.96*linearitySmooth$se.fit, 
              col="red", lty="dashed")
        
      }
      
    }
    
  }
  
  dev.off()
  
}


if (FALSE) {
  
  # check of lambda parameter tuning for ridge and LASSO regression
  
  formulations <- expand.grid(effects=c("main", "interaction"), 
                              types=c("Ridge", "LASSO"), 
                              transformations=c("Log", "Inv"),
                              stringsAsFactors=FALSE)
  balancing <- rep("Balanced", times=nrow(formulations))
  
  
  for (i in 1:nrow(formulations)) {
    formulation <- formulations[i, ]
    
    # specify alpha for the type of regression
    if (type=="Ridge") {alpha <- 0}
    if (type=="LASSO") {alpha <- 1}
    
    lambdaWhite <- 
      tuneLambda(data_matrix=makeMatrix(data_white_train_balanced, 
                                        includeInteraction=formulation$effect=="interaction")$mat, 
                 alpha= alpha, 
                 new.y=data_white_train_balanced$BMI, 
                 transformation=formulation$transformation, 
                 effect=formulation$effect, returnAll=TRUE)
    bestLambdaWhite <- lambdaWhite$lambda[which.max(lambdaWhite$AUC)]
    p_tuneLambdaWhite <- ggplot(lambdaWhite, aes(x=lambda, y=AUC)) +
      geom_point() + scale_x_log10() +
      geom_vline(xintercept=bestLambdaWhite, col="blue") +
      labs(title=paste0("Demonstration of parameter tuning: ", formulation$type, " regression"),
           subtitle=paste0("Ethnicity: white, effects: ", formulation$effect, 
                           ", BMI transformation: ", formulation$transformation))
    print(p_tuneLambdaWhite)
  }
  
  for (i in 1:nrow(formulations)) {
    formulation <- formulations[i, ]
    
    # specify alpha for the type of regression
    if (type=="Ridge") {alpha <- 0}
    if (type=="LASSO") {alpha <- 1}
    
    lambdaBlack <- 
      tuneLambda(data_matrix=makeMatrix(data_black_train, 
                                        includeInteraction=formulation$effect=="interaction")$mat, 
                 alpha= alpha, 
                 new.y=data_black_train$BMI, 
                 transformation=formulation$transformation, 
                 effect=formulation$effect, returnAll=TRUE)
    bestLambdaBlack <- lambdaBlack$lambda[which.max(lambdaBlack$AUC)]
    p_tuneLambdaBlack <- ggplot(lambdaBlack, aes(x=lambda, y=AUC)) +
      geom_point() + scale_x_log10() +
      geom_vline(xintercept=bestLambdaBlack, col="blue") +
      labs(title=paste0("Demonstration of parameter tuning: ", formulation$type, " regression"),
           subtitle=paste0("Ethnicity: black, effects: ", formulation$effect, 
                           ", BMI transformation: ", formulation$transformation))
    print(p_tuneLambdaBlack)
  }
  
}



