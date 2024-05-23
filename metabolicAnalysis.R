trainAllOLS <- function(new.data, modelstring) {
  
  model <- eval(parse(text=modelstring))
  
  if (grepl(x=modelstring, pattern="OLS")) {type <- "OLS"}
  if (grepl(x=modelstring, pattern="LASSO")) {type <- "LASSO"}
  if (grepl(x=modelstring, pattern="Ridge")) {type <- "Ridge"}
  
  stopifnot("Ridge regression models are not suitable for this function" = type!="Ridge")
  
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
  
  if (type=="LASSO") {
    w <- which(model$beta[,"s0"]==0)
    variables <- names(model$beta[-w,"s0"])
    ratioL <- which(grepl(x=variables, pattern="/"))
    variables[ratioL] <- paste0("`", variables[ratioL], "`")
    alias_model <- lm(data=new.data,
                      formula = eval(parse(text=paste0(outcome, "~", 
                                                       paste(variables, collapse="+")))))
    a <- which(alias(alias_model)$Complete[1,]!=0)
    b <- which(variables %in% rownames(alias(alias_model)$Complete))
    
    vars0 <- variables[-c(a,b)]
    varsA <- variables[c(a,b)]
    
    startModel <- lm(data=new.data,
                     formula = eval(parse(text=paste0(outcome, "~", 
                                                      paste(vars0, collapse="+")))))
    # stepwise regression will end as soon as the worst predictive aliased variables remain
    return(step(startModel, direction='forward', scope=formula(alias_model), trace=0))
  }
  
  
}

simulateAlternative <- function(model, data, transformation, nsim=1000, oversampled=FALSE) {
  
  if (transformation=="Log") {data$transBMI <- data$BMI}
  if (transformation=="Inv") {data$transBMI <- exp(-data$BMI)}
  
  vars <- names(model$coefficients)[-which(names(model$coefficients)=="(Intercept)")]
  formula <- paste0("transBMI~", paste0(vars, collapse="+"), "")
  formula <- gsub(pattern=":", replacement="*", formula, fixed=TRUE)
  
  if (!oversampled) {
    preds <- model$fitted.values
    resid <- model$residuals
  }
  
  if (oversampled) {
    preds <- predict(model, newdata=data)
    resid <- data$transBMI - preds
  }
  
  sigma <- sd(resid)
  
  if (shapiro.test(resid)$p.value < 1e-4) warning("residuals are not normally distributed")
  
  pvals <- matrix(ncol=length(vars), nrow=nsim)
  colnames(pvals) <- vars
  
  for (i in 1:nsim) {
    set.seed(i)
    data_new <- data
    data_new$transBMI <- preds+rnorm(n=nrow(data), mean=0, sd=sigma)
    
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

multipleCorrTest <- function(data, predictions, band_lower, band_upper, 
                             nsim0=1000, new.n=NULL, returnPower=FALSE, nsimP=1000) {
  w <- which(predictions>band_lower & predictions<band_upper)
  use_data <- data[w,]
  
  # if no sample size is specified, perform hypothesis test under given sample size
  if (is.null(new.n)) {new.n <- 1}
  
  # calculate most extreme t-statistic of sample
  t_sample <- 0
  for (met in metabolites) {
    r <- cor.test(x=use_data[[met]], y=use_data$BMI, alternative="two.sided", 
                  method="spearman", exact=FALSE)$estimate
    # adapt test to hypothetical sample size
    t <- r*sqrt((new.n*nrow(use_data)-2)/(1-r**2))
    t_sample <- max(t_sample, abs(t))
  }
  
  # If power of sample effects is requested, simulate alternative as sample point estimates
  # since metabolites and according residuals can be correlated, sample replicates 
  # are constituted as the predicted metabolite values given the BMI, injected with
  # a joint series of residuals as found from one patient in the sample.
  if (returnPower) {
    predicted_mets <- matrix(ncol=length(metabolites), nrow=nrow(use_data))
    colnames(predicted_mets) <- metabolites
    joint_residuals <- matrix(ncol=length(metabolites), nrow=nrow(use_data))
    colnames(joint_residuals) <- metabolites
    for (met in metabolites) {
      reg <- lm(formula=eval(parse(text=paste(met, "~ BMI"))), data=use_data)
      predicted_mets[,met] <- reg$fitted.values
      joint_residuals[,met] <- reg$residuals
    }
    
    simAlt <- function(new.n, use_data, predicted_mets, joint_residuals) {
      #set.seed(seed+0.5)
      t_alt_i <- vector(mode="numeric", length=length(new.n))
      names(t_alt_i) <- as.character(new.n)
      for (n in new.n) {
        v <- rep(1:nrow(use_data), times=n)
        w <- sample(1:nrow(use_data), size=n*nrow(use_data), replace=TRUE)
        pred <- predicted_mets[v,]
        res <- joint_residuals[w,]
        alt_data <- as.data.frame(pred+res)
        alt_data$BMI <- use_data$BMI[v]
        t_n <- 0
        for (met in metabolites) {
          r <- cor.test(x=alt_data[[met]], y=alt_data$BMI, alternative="two.sided", 
                        method="spearman", exact=FALSE)$estimate
          t <- r*sqrt((n*nrow(use_data)-2)/(1-r**2))
          t_n <- max(t_n, abs(t))
        }
        t_alt_i[as.character(n)] <- t_n
      }
      t_alt_i
    }
    
    t_alt <- replicate(nsimP, simAlt(new.n, use_data, predicted_mets, joint_residuals))
    if (length(new.n)==1) {
      dim(t_alt) <- c(length(new.n), length(t_alt))
      rownames(t_alt) <- as.character(new.n)
    }
    
    
  }
  
  # constitute a simulation null distribution of lowest t-statistics
  simNull <- function(new.n, use_data) {
    #set.seed(seed)
    t_0_i <- vector(mode="numeric", length=length(new.n))
    names(t_0_i) <- as.character(new.n)
    for (n in new.n) {
      v <- rep(1:nrow(use_data), times=n)
      w <- sample(1:(nrow(use_data)*n), size=n*nrow(use_data), replace=FALSE)
      null_data <- use_data[v,]
      null_data$BMI <- rep(use_data$BMI, times=n)[w]
      t_n <- 0
      for (met in metabolites) {
        r <- cor.test(x=null_data[[met]], y=null_data$BMI, alternative="two.sided", 
                      method="spearman", exact=FALSE)$estimate
        t <- r*sqrt((n*nrow(use_data)-2)/(1-r**2))
        t_n <- max(t_n, abs(t))
      }
      t_0_i[as.character(n)] <- t_n
    }
    t_0_i
  }
  
  t_0 <- replicate(nsim0, simNull(new.n, use_data))
  if (length(new.n)==1) {
    dim(t_0) <- c(length(new.n), length(t_0))
    rownames(t_0) <- as.character(new.n)
  }
  
  # Return a corrected p-value for the sample if null hypothesis is tested
  if (!returnPower) {return(mean(t_0 > t_sample))}
  
  # return series of corrected p-values from the alternative hypothesis simulation
  # if power of the alternative hypothesis gets tested.
  if (returnPower) {
    returns <- list()
    for (n in as.character(new.n)) {
      returns[[n]] <- sapply(1:ncol(t_alt), function(i) mean(t_0[n,] < t_alt[n,i]))
    }
    return(as.data.frame(returns, check.names=FALSE))
  }
}




