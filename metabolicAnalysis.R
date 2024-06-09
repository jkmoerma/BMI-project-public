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
  
  data$predicted <- NULL
  
  if (transformation=="Log") {data$transBMI <- data$BMI}
  if (transformation=="Inv") {data$transBMI <- exp(-data$BMI)}
  
  # extract predictive variables from the model
  if (type=="Ridge"|type=="LASSO") {vars0 <- names(model$beta[-which(model$beta[,"s0"]==0),])}
  if (type=="OLS") {vars0 <- names(model$coefficients)[-1]}
  
  # retrieve residuals 
  resid <- data$transBMI - preds
  sigma <- sd(resid)
  if (shapiro.test(resid)$p.value < 1e-4) warning(sprintf("residuals are not normally distributed, Shapiro-Wilk p-val: %.2e", shapiro.test(resid)$p.value))
  
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
    
    t_oversampled <- summary(repModel)$coefficients[vars,"t value"]
    t2_corrected <- t_oversampled**2 * nrow(data)/nrow(data_new)
    pvals[i, vars] <- 1-pchisq(t2_corrected, df=1)
    
  }
  
  betas <- cbind("simulation"=1:nsim, betas)
  pvals <- cbind("simulation"=1:nsim, pvals)
  
  betas <- as.data.frame(betas)
  pvals <- as.data.frame(pvals)
  
  list(betas=betas, pvals=pvals)
  
}


calculatePower <- function(vec, alpha=0.05) {
  mean(vec<alpha)
}

plotPower <- function(simulationPvals, plottitle, power_threshold=0.5,
                      alphas=c(0.0001, 0.001, 0.01, 0.05, 0.1)) {
  power <- c()
  for (alpha in alphas) {
    power <- bind_rows(power, c(alpha=alpha, 
                                unlist(lapply(X=simulationPvals, FUN=calculatePower, 
                                              alpha=alpha))))
  }
  
  plotdata <- pivot_longer(power, cols=colnames(simulationPvals), 
                           names_to="met", values_to="power")
  line_legend <- data.frame(met=colnames(simulationPvals), 
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
    scale_x_continuous(trans='log10') +
    labs(title=plottitle)
}

#' Calculate power metabolite beta-coefficients related to obesity for different sample sizes
#'
#' @param betas A data frame containing beta coefficients from a simulation under the alternative hypothesis
#' @param sample_sizes A vector containing the sample sizes under which the power for every metabolite beta-coefficient needs to be calculated
#' @param n_sample Integer number representing the sample size of the sample used in the simulation experiment
#' @param alpha Number representing the significance level to reject the null hypothesis
#' @return A data frame with the probed sample sizes and the power of metabolite beta-coefficients under given significance level
#' @examples
#' 
sampleSizeObeseMetabolites <- function(betas, sample_sizes, n_sample, alpha) {
  powers <- matrix(ncol=length(colnames(betas))-1, nrow=length(sample_sizes))
  rownames(powers) <- as.character(sample_sizes)
  colnames(powers) <- colnames(betas)[-1]
  Tc <- qnorm(p=c(alpha/2, 1-alpha/2))
  for (met in colnames(betas)[-1]) {
    coeffs <- betas[[met]]
    average <- mean(coeffs)
    stdev <- sd(coeffs)
    for (n in sample_sizes) {
      TA <- sqrt(n)*coeffs/(stdev*sqrt(n_sample))
      power <- (length(which(TA<Tc[1])) + length(which(TA>Tc[2])))/length(TA)
      powers[as.character(n), met] <- power
    }
  }
  powers <- as.data.frame(powers)
  powers$n <- sample_sizes
  powers
}

plotSampleSize <- function(powers, requested_power, at_samplesize, plottitle) {
  pivotPowers <- pivot_longer(powers, cols=colnames(powers)[-which(colnames(powers)=="n")],
                              names_to="met", values_to="power")
  line_legend <- subset(pivotPowers, subset= n==at_samplesize)
  line_legend$lab <- "other"
  line_legend$lab[which(line_legend$power>requested_power)] <- 
    line_legend$met[which(line_legend$power>requested_power)]
  line_legend$power <- NULL
  line_legend$n <- NULL
  
  pivotPowers <- merge(pivotPowers, line_legend, ID="met")
  ggplot(pivotPowers, aes(x=n, y=power, by=met)) + 
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
      t_0_i <- max(t_0_i, abs(t))
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




if (FALSE) {
  
  # Verification statistical power computation alternative hypothesis
  
  # This code checks the calculation of p-values from replicated simulations of 
  # the alternative hypothesis.
  # The Wald p-values stored from the simulation experiment were compared with 
  # the p-values computed from the simulation repeatedly fitted beta-coefficients
  
  
  # The analysis for black ethnicity patients was never oversampled.
  # This first test section inspects p-value calculation when no oversampling 
  # nor increasing sample size was applied.
  
  blackSimulation <- simulateAlternative(modelBlackTrain, data_black_train, 1/blackPredictions,
                                         transformation="Inv", type="Ridge", nsim=20)
  
  simBetasBlack <- pivot_longer(blackSimulation$betas, cols=colnames(blackSimulation$betas)[-1], 
                                names_to="met", values_to="beta")
  standardize <- "(beta)/sd(beta)"
  pWald <- paste0("1-pchisq((", standardize, ")**2, df=1)")
  simBetasBlack <- simBetasBlack %>% group_by(met) %>% 
    mutate(standardized=eval(parse(text=standardize)),
           "Replicate p-value"=eval(parse(text=pWald)))
  simPvalsBlack <- pivot_longer(blackSimulation$pvals, cols=colnames(blackSimulation$betas)[-1], 
                                names_to="met", values_to="Wald p-value")
  comparePvalsBlack <- merge(simBetasBlack, simPvalsBlack, key=c("simulation", "met"))
  
  ggplot(comparePvalsBlack, aes(x=`Replicate p-value`, y=`Wald p-value`)) + 
    geom_point(aes(col=met)) + 
    geom_abline(slope=1, intercept=0, lty="dashed", col="gray") +
    scale_x_continuous(trans='log10')+
    scale_y_continuous(trans='log10')+
    labs(title="Concordance of replicated p-val. with Wald p-val.",
         subtitle="No changed sample size, no balancing")
  
  
  # The analysis for white ethnicity patients was oversampled.
  # This second test section inspects p-value calculation when oversampling 
  # without increasing sample size was applied.
  
  whiteSimulation <- simulateAlternative(modelWhiteTrain, data_white_train, 1/whitePredictions,
                                         transformation="Inv", type="LASSO", nsim=20, oversampled=TRUE)
  
  simBetasWhite <- pivot_longer(whiteSimulation$betas, cols=colnames(whiteSimulation$betas)[-1], 
                                names_to="met", values_to="beta")
  standardize <- "(beta)/sd(beta)"
  pWald <- paste0("1-pchisq((", standardize, ")**2, df=1)")
  simBetasWhite <- simBetasWhite %>% group_by(met) %>% 
    mutate(standardized=eval(parse(text=standardize)),
           "Replicate p-value"=eval(parse(text=pWald)))
  simPvalsWhite <- pivot_longer(whiteSimulation$pvals, cols=colnames(whiteSimulation$betas)[-1], 
                                names_to="met", values_to="Wald p-value")
  comparePvalsWhite <- merge(simBetasWhite, simPvalsWhite, key=c("simulation", "met"))
  
  selectLegend <- comparePvalsWhite %>% group_by(met) %>% summarise(select=min(`Wald p-value`)<1e-4)
  selectLegend$legend <- "other"
  selectLegend$legend[which(selectLegend$select)] <- selectLegend$met[which(selectLegend$select)]
  selectLegend$select <- NULL
  comparePvalsWhite <- merge(comparePvalsWhite, selectLegend, key="met")
  
  ggplot(comparePvalsWhite, aes(x=`Replicate p-value`, y=`Wald p-value`)) + 
    geom_point(aes(col=legend)) + 
    geom_abline(slope=1, intercept=0, lty="dashed", col="gray") +
    scale_x_continuous(trans='log10')+
    scale_y_continuous(trans='log10', limits=c(1e-16, 1))+
    labs(title="Concordance of replicated p-val. with Wald p-val.",
         subtitle="No changed sample size, oversampling in simulation")
  
  
  # The analysis for black ethnicity patients was not oversampled.
  # This third test section inspects p-value estimation for no oversampling 
  # but increased sample size.
  # Distribution of p-values are compared for increasing the sample size and 
  
  blackSimulation2 <- simulateAlternative(modelBlackTrain, 
                                          bind_rows(data_black_train, data_black_train), 
                                          rep(1/blackPredictions, times=2),
                                          transformation="Inv", type="Ridge", nsim=1000)
  
  blackSimulation <- simulateAlternative(modelBlackTrain, data_black_train, 1/blackPredictions,
                                         transformation="Inv", type="Ridge", nsim=1000)
  
  pvals1 <- sampleSizeObeseMetabolites(blackSimulation2$betas, 2*nrow(data_black_train), 
                             2*nrow(data_black_train), alpha=0.05)
  pvals2 <- sampleSizeObeseMetabolites(blackSimulation$betas, 2*nrow(data_black_train), 
                             nrow(data_black_train), alpha=0.05)
  View(bind_rows(pvals1, pvals2))
  
}

