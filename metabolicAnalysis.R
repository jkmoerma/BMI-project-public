trainAllOLS <- function(new.data, oversampled=FALSE) {
  
  if (oversampled) {
    n0 <- nrow(new.data)
    new.data <- oversample(new.data)
    amount <- nrow(new.data) - n0
  }
  
  intercept_only <- lm(exp(-BMI) ~ 1, data=new.data)
  all <- lm(data=new.data,
            formula = eval(parse(text=paste0("exp(-BMI)~`", paste(c("Age", "Smoking", metabolites), 
                                                                  collapse="`+`"), "`"))))
  mainOLSModel <- step(intercept_only, direction='both', scope=formula(all), trace=0)
  
  # interaction effects, 1/BMI, ethnicity="White"
  variables <- 
    colnames(mainOLSModel$model)[-which(colnames(mainOLSModel$model)=="exp(-BMI)")]
  interactionSet <- c()
  for (i in 1:(length(variables)-1)) {
    for (j in (i+1):length(variables)) {
      term <- paste0(variables[i], "`*`", variables[j])
      addedModel <- lm(data=new.data,
                       formula=eval(parse(text=paste0("exp(-BMI)~`", 
                                                      paste0(c(variables, term), 
                                                             collapse="`+`"), "`"))))
      AIC0 <- aic(mainOLSModel, oversampled=oversampled, amount=amount)
      AIC1 <- aic(addedModel, oversampled=oversampled, amount=amount)
      if (AIC1 < AIC0) {interactionSet <- c(interactionSet, term)}
    }
  }
  all <- lm(data=new.data,
            formula = eval(parse(text=paste0("exp(-BMI)~`", paste0(c(variables, interactionSet), 
                                                                   collapse="`+`"), "`"))))
  step(all, direction='both', scope=formula(all), trace=0)
}

simulateAlternative <- function(model, data, nsim=1000, oversampled=FALSE) {
  
  if (!oversampled) {
    vars <- names(model$coefficients)[-which(names(model$coefficients)=="(Intercept)")]
    formula <- paste0("exp(-BMI)~", paste0(vars, collapse="+"), "")
    formula <- gsub(pattern=":", replacement="*", formula, fixed=TRUE)
    preds <- model$fitted.values
    resid <- model$residuals
  }
  
  if (oversampled) {
    vars <- names(model$coefficients)[-which(names(model$coefficients)=="(Intercept)")]
    formula <- paste0("exp(-BMI)~", paste0(vars, collapse="+"), "")
    formula <- gsub(pattern=":", replacement="*", formula, fixed=TRUE)
    preds <- predict(model, newdata=data)
    resid <- exp(-data$BMI) - preds
  }
  
  sigma <- sd(resid)
  
  if (shapiro.test(resid)$p.value < 1e-4) warning("residuals are not normally distributed")
  
  pvals <- matrix(ncol=length(vars), nrow=nsim)
  colnames(pvals) <- vars
  
  for (i in 1:nsim) {
    data_new <- data
    data_new$BMI <- -log(preds+rnorm(n=nrow(data), mean=0, sd=sigma))
    
    if (oversampled) {data_new <- oversample(data_new)}
    
    repModel <- lm(formula=eval(parse(text=formula)), data=data_new)
    pvals[i, vars] <- summary(repModel)$coefficients[vars,"Pr(>|t|)"]
    
  }
  
  as.data.frame(pvals)
  
}


calculatePower <- function(vec, alpha=0.05) {
  mean(vec<alpha)
}

