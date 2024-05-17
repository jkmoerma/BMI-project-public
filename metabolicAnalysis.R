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
    set.seed(i)
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

plotPower <- function(simulation, plottitle, power_threshold=0.5,
                      alphas=c(0.0001, 0.001, 0.01, 0.05, 0.1)) {
  power <- c()
  for (alpha in alphas) {
    power <- bind_rows(power, c(alpha=alpha, 
                                unlist(lapply(X=simulation, FUN=calculatePower, 
                                              alpha=alpha))))
  }
  
  plotdata <- pivot_longer(power, cols=colnames(simulation), 
                           names_to="met", values_to="power")
  line_legend <- data.frame(met=colnames(simulation), 
                            power0=unlist(as.vector(power[1,-1]))) %>% 
                   mutate(leg=1+(power0<power_threshold))
  line_legend$lab <- "other"
  line_legend$lab[which(line_legend$leg==1)] <- 
    line_legend$met[which(line_legend$leg==1)]
  line_legend$power0 <- NULL
  line_legend$leg <- NULL
  
  plotdata <- merge(plotdata, line_legend, ID="met")
  ggplot(data=plotdata, 
         aes(x=alpha, y=power, by=met, logscale="x")) + 
    geom_line(aes(color=lab)) +
    scale_x_continuous(trans='log10')+
    labs(title=plottitle)
}

#plotMetabolicBands <- function() {}

multipleANOVA <- function(data, predictions, band_lower, band_upper, nsim=1000) {
  w <- which(predictions>band_lower & predictions<band_upper)
  use_data <- data[w,]
  
  # calculate lowest p-value of sample
  p_sample <- 1
  for (met in metabolites) {
    sampleModel <- lm(formula=eval(parse(text=paste0("exp(-BMI)~`", met, "`"))), 
                      data=use_data)
    colname <- met
    if (grepl(pattern="/", x=met)) {colname <- paste0("`", met, "`")}
    p <- 1-pf(summary(sampleModel)$fstatistic, df1=1, df2=sampleModel$df.residual)
    p_sample <- min(p_sample, p)
  }
  
  
  # constitute a simulation null distribution
  p0 <- vector(mode="numeric", length=nsim)
  for (i in 1:nsim) {
    p_i <- 1
    set.seed(i)
    null_data <- use_data
    null_data$BMI <- sample(use_data$BMI, size=nrow(use_data), replace=FALSE)
    for (met in metabolites) {
      nullModel <- lm(formula=eval(parse(text=paste0("exp(-BMI)~`", met, "`"))), 
                      data=null_data)
      colname <- met
      if (grepl(pattern="/", x=met)) {colname <- paste0("`", met, "`")}
      p <- 1-pf(summary(nullModel)$fstatistic, df1=1, df2=nullModel$df.residual)
      p_i <- min(p_i, p)
    }
    p0[i] <- p_i
  }
  
  # return a corrected p-value for the sample
  mean(p0 < p_sample)
}




