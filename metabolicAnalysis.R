trainAllOLS <- function(new.data, modelstring) {
  
  model <- eval(parse(text=modelstring))
  
  if (grepl(x=modelstring, pattern="OLS")) {type <- "OLS"}
  if (grepl(x=modelstring, pattern="LASSO")) {type <- "LASSO"}
  if (grepl(x=modelstring, pattern="Ridge")) {type <- "Ridge"}
  
  if (grepl(x=modelstring, pattern="Log")) {outcome <- "BMI"}
  if (grepl(x=modelstring, pattern="Inv")) {outcome <- "exp(-BMI)"}
  
  if (grepl(x=modelstring, pattern="Balanced")) {new.data <- oversample(new.data)}
  
  if (type=="OLS") {
    variables <- 
      names(model$coefficients)[-which(names(model$coefficients)=="(Intercept)")]
    return(lm(data=new.data,
              formula = eval(parse(text=paste0(outcome, "~", 
                                               paste(variables, collapse="+"))))))
  }
  
  if (type=="LASSO" | type=="Ridge") {
    w <- which(model$beta[,"s0"]==0)
    if (length(w)==0) {
      variables <- names(model$beta[,"s0"])
    } else {
      variables <- names(model$beta[-w,"s0"])
    }
    ratioL <- which(grepl(x=variables, pattern="/"))
    variables[ratioL] <- paste0("`", variables[ratioL], "`")
    alias_model <- lm(data=new.data,
                      formula = eval(parse(text=paste0(outcome, "~", 
                                                       paste(variables, collapse="+")))))
    a <- which(variables %in% colnames(alias(alias_model)$Complete))
    b <- which(variables %in% rownames(alias(alias_model)$Complete))
    
    vars0 <- variables[-c(a,b)]
    if (length(vars0)==0) {
      startModel <- lm(data=new.data,
                       formula = eval(parse(text=paste0(outcome, "~", "1"))))
    } else {
      startModel <- lm(data=new.data,
                       formula = eval(parse(text=paste0(outcome, "~", 
                                                        paste(vars0, collapse="+")))))
    }
    
    
    # "eat" residual variance of the model alias structure by alias structure
    for (i in nrow(alias(alias_model)$Complete)) {
      
      varsA <- rownames(alias(alias_model)$Complete)[i]
      varsB <- colnames(which(alias(alias_model)$Complete[i,]!=0))
      
      aliasModel <- lm(data=new.data,
                       formula = eval(parse(text=paste0(outcome, "~", 
                                                        paste(c(vars0, varsA, varsB), collapse="+")))))
      startModel <- step(startModel, direction='forward', scope=formula(aliasModel), trace=0)
      vars0 <- colnames(startModel$model)[-which(colnames(startModel$model)==outcome)]
    }
    
    # stepwise regression will end as soon as the worst predictive aliased variables remain
    return(startModel)
  }
  
  
}

simulateAlternative <- function(model, data, preds, transformation, type,
                                interactionEffect=FALSE, nsim=1000, oversampled=FALSE) {
  
  if (transformation=="Log") {data$transBMI <- data$BMI}
  if (transformation=="Inv") {data$transBMI <- exp(-data$BMI)}
  
  # extract predictive variables from the model
  vars0 <- names(model$beta[which(model$beta!=0),])
  
  # retrieve residuals 
  resid <- data$transBMI - preds
  sigma <- sd(resid)
  if (shapiro.test(resid)$p.value < 1e-4) warning("residuals are not normally distributed")
  
  # make best subselection of variables in which there is no alias structure
  # do this by building a stepwise linear regression model with the variables
  formula0 <- paste0("transBMI~", paste0(vars0, collapse="+"), "")
  formula0 <- gsub(pattern=":", replacement="*", formula0, fixed=TRUE)
  
  intercept_only <- lm(transBMI ~ 1, data=data)
  all <- lm(data=data, formula=eval(parse(text=formula0)))
  model <- step(intercept_only, direction='both', scope=formula(all), trace=0)
  
  
  # extract variables that are selected in the stepwise regression model
  vars <- names(model$coefficients)[-which(names(model$coefficients)=="(Intercept)")]
  
  formula <- paste0("transBMI~", paste0(vars, collapse="+"), "")
  formula <- gsub(pattern=":", replacement="*", formula, fixed=TRUE)
  
  betas <- matrix(ncol=length(vars), nrow=nsim)
  colnames(betas) <- vars
  
  pvals <- matrix(ncol=length(vars), nrow=nsim)
  colnames(pvals) <- vars
  
  for (i in 1:nsim) {
    
    set.seed(i)
    data_new <- data
    data_new$transBMI <- preds+rnorm(n=nrow(data), mean=0, sd=sigma)
    
    if (oversampled) {data_new <- oversample(data_new)}
    
    repModel <- lm(formula=eval(parse(text=formula)), data=data_new)
    betas[i, vars] <- summary(repModel)$coefficients[vars,"Estimate"]
    pvals[i, vars] <- summary(repModel)$coefficients[vars,"Pr(>|t|)"]
    
  }
  
  list(
    betas = as.data.frame(betas),
    pvals = as.data.frame(pvals)
  )
  
  
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


subsetMetabolicBands <- function(data, predictions, band_lower, band_upper) {
  w <- which(predictions>band_lower & predictions<band_upper)
  use_data <- data[w,]
  use_data
}


multipleCorrTest <- function(use_data, nsim0=1000) {
  
  # calculate most extreme t-statistic of sample
  t_sample <- 0
  for (met in metabolites) {
    r <- cor.test(x=use_data[[met]], y=use_data$BMI, alternative="two.sided", 
                  method="spearman", exact=FALSE)$estimate
    # adapt test to hypothetical sample size
    t <- r*sqrt((nrow(use_data)-2)/(1-r**2))
    t_sample <- max(t_sample, abs(t))
  }
  
  
  # constitute a simulation null distribution of lowest t-statistics
  simNull <- function(use_data) {
    #set.seed(seed)
    w <- sample(1:nrow(use_data), size=nrow(use_data), replace=FALSE)
    null_data <- use_data
    null_data$BMI <- use_data$BMI[w]
    t_0_i <- 0
    for (met in metabolites) {
      r <- cor.test(x=null_data[[met]], y=null_data$BMI, alternative="two.sided", 
                    method="spearman", exact=FALSE)$estimate
      t <- r*sqrt((nrow(use_data)-2)/(1-r**2))
      t_0_i <- max(t_0, abs(t))
    }
    t_0_i
  }
  
  t_0 <- replicate(nsim0, simNull(use_data))
  
  # Return a corrected p-value for the sample if null hypothesis is tested
  mean(t_0 > t_sample)
  
}


powerCorrTest <- function(use_data, n, nsim0=1000, nsimP=1000, alpha=0.05) {
  
  # If power of sample effects is requested, simulate alternative as sample point estimates
  # since metabolites and according residuals can be correlated, sample replicates 
  # are constituted as the predicted metabolite values given the BMI, injected with
  # a joint series of residuals as found from one patient in the sample.
  predicted_mets <- matrix(ncol=length(metabolites), nrow=nrow(use_data))
  colnames(predicted_mets) <- metabolites
  joint_residuals <- matrix(ncol=length(metabolites), nrow=nrow(use_data))
  colnames(joint_residuals) <- metabolites
  for (met in metabolites) {
    reg <- lm(formula=eval(parse(text=paste(met, "~ BMI"))), data=use_data)
    predicted_mets[,met] <- reg$fitted.values
    joint_residuals[,met] <- reg$residuals
  }
  
  p <- vector(mode="numeric", length=nsimP)
  for (i in 1:nsimP) {
    # simulate a sample under alternative hypothesis
    v <- sample(1:nrow(use_data), size=n, replace=TRUE)
    w <- sample(1:nrow(use_data), size=n, replace=TRUE)
    pred <- predicted_mets[v,]
    res <- joint_residuals[w,]
    alt_data <- as.data.frame(pred+res)
    alt_data$BMI <- use_data$BMI[v]
    
    # calculate p-value under the alternative
    p_i <- multipleCorrTest(alt_data, nsim0)
    
    # store
    p[i] <- p_i
  }
  
  mean(p<alpha)
  
}




