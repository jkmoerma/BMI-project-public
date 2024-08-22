
#' Calculates regression coefficients of a model formulation and bootstrap replicates of the data set. 
#' The predictors are scaled to unit variance in order to compare the regression coefficients of different predictors equally
#' 
#' @param data Data to draw bootstrap replicates from
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @param balancing Defaults to "", not oversampling the data. Specify "Balanced" if oversampling is needed.
#' @param boot.n the amount of bootstrap replicates to be drawn
#' @return A data tibble with regression coefficients for every replicate.
#'
scaledEffects <- function(data, effect, type, transformation, balancing="", 
                          boot.n=100) {
  
  # if balancing the data set was part of the model formulation, scaling data on an oversampled replicate was necessary for ...
  if (balancing=="") {data_ref <- data}
  if (balancing=="Balanced") {data_ref <- oversample(data)}
  
  # specify alpha for the type of regression
  if (type=="Ridge") {alpha <- 0}
  if (type=="LASSO") {alpha <- 1}
  
  # generate a scaled version of the data set 
  # do not standardize BMI on type of transformation
  data_scaled <- data
  for (met in c("Age", metabolites)) {
    data_scaled[[met]] <- data_scaled[[met]]/sd(data_ref[[met]])
  }
  
  if (balancing=="") {
    tune_outcome <- data_scaled$BMI
    dataMatrix <- makeMatrix(data_scaled, includeInteraction = effect=="interaction")$mat
  }
  if (balancing=="Balanced") {
    tune_data <- oversample(data_scaled)
    tune_outcome <- tune_data$BMI
    dataMatrix <- makeMatrix(tune_data, includeInteraction = effect=="interaction")$mat
  }
  
  lambda <- tuneLambda(dataMatrix, alpha = alpha, new.y=tune_outcome, 
                       effect, transformation, returnAll=FALSE)
  
  # for boot.n bootstrap replicates, recalculate the regression coefficients
  # Register parallel backend
  numCores <- detectCores() - 1  # Use one less than the total number of cores
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  # Parallelized loop
  coeffs <- foreach(i=1:boot.n, .combine=bind_rows, .packages=c("dplyr", "glmnet"), 
                    .export=c("trainOLS", "trainRidgeLASSO", 
                              "metabolites", "aic", "oversample", 
                              "SMOTE", "makeMatrix", "roc", "auc")) %dopar% {
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
        model <- trainRidgeLASSO(effect, type, transformation, new.x=new.x, 
                                 new.y=new.y, lambda=lambda)
      }
      if (balancing=="Balanced") {
        new.x <- makeMatrix(data_rep, includeInteraction = effect=="interaction")$mat
        new.y <- data_rep$BMI
        model <- trainRidgeLASSO(effect, type, transformation, new.x=new.x, 
                                 new.y=new.y, lambda=lambda)
      }
      return(model$beta[,"s0"])
    }
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # return coefficients
  return(coeffs)
  
}

#' Visualize the regression coefficients + 95 pc CI of the bootstrap replicates.
#' The groups of correlated metabolites were specified from the data inspection and are not customizable in this version of the function.
#' 
#' @param bootstrapCoeffsWhite Bootstrap replicated coefficients of the white patients data
#' @param bootstrapCoeffsBlack Bootstrap replicated coefficients of the black patients data
#' @param effect "main" or "interaction" dependent on which effects are wished
#' @param type Regression type, must be "Ridge" or "LASSO for this function
#' @param transformation "Log" or "Inv", dependent on whether regression on log(BMI) or 1/BMI is required
#' @return If model only uses main effects, a ggplot2 plot is returned of the main effects, if interaction effects used, a list of two plots for main effects and 50 most important interaction effects.
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
#' data_normal_white <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=5, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight_white <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=4, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese_white <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                         
#' data_train_white <- bind_rows(data_normal_white, data_overweight_white, data_obese_white)
#' data_train_white$BMI <- log(data_train_white$BMI)
#' data_train_white$met_coll <- data_train_white$met1 + 7.2*data_train_white$met2 + rnorm(n=nrow(data_train_white), mean=0, sd=0.1)
#' 
#' data_normal_black <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=7, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight_black <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=5, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese_black <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                         
#' data_train_black <- bind_rows(data_normal_black, data_overweight_black, data_obese_black)
#' data_train_black$BMI <- log(data_train_black$BMI)
#' data_train_black$met_coll <- data_train_black$met1 + 7.2*data_train_black$met2 + rnorm(n=nrow(data_train_black), mean=0, sd=0.1)
#' 
#' 
#' bootstrapCoeffsWhite <- scaledEffects(data_train_white, effect="main", 
#'                                       type="LASSO", transformation="Log", 
#'                                       balancing="Balanced", boot.n=1000)
#' bootstrapCoeffsWhite <- pivot_longer(bootstrapCoeffsWhite, 
#'                                      cols=colnames(bootstrapCoeffsWhite),
#'                                      names_to="met", values_to="coeff")
#' bootstrapCoeffsBlack <- scaledEffects(data_train_black, effect="main", 
#'                                       type="LASSO", transformation="Log", 
#'                                       balancing="Balanced", boot.n=1000)
#' bootstrapCoeffsBlack <- pivot_longer(bootstrapCoeffsBlack, 
#'                                      cols=colnames(bootstrapCoeffsBlack),
#'                                      names_to="met", values_to="coeff")
#' 
#' plotScaledEffects(bootstrapCoeffsWhite, bootstrapCoeffsBlack, 
#'                   effect="main", type="LASSO", transformation="Log")
#' 
#' 
plotScaledEffects <- function(bootstrapCoeffsWhite, bootstrapCoeffsBlack, 
                              effect, type, transformation) {
  
  # generate coefficient estimates for main effects in white ethnicity
  CIcoeffsWhite <- bootstrapCoeffsWhite %>% 
    group_by(met) %>% summarise("est."=quantile(x=coeff, probs=0.5),
                                "95pc lCI"=quantile(x=coeff, probs=0.025),
                                "95pc uCI"=quantile(x=coeff, probs=0.975))
  CIcoeffsWhite["Race"] <- "White"
  
  # generate coefficient estimates for main effects in black ethnicity
  CIcoeffsBlack <- bootstrapCoeffsBlack %>% 
    group_by(met) %>% summarise("est."=quantile(x=coeff, probs=0.5),
                                "95pc lCI"=quantile(x=coeff, probs=0.025),
                                "95pc uCI"=quantile(x=coeff, probs=0.975))
  CIcoeffsBlack["Race"] <- "Black"
  
  # check arguments given to function
  stopifnot("beta coefficients must originate from same modeling formulation" = 
              nrow(CIcoeffsWhite)==nrow(CIcoeffsBlack))
  if (effect=="main") {
    stopifnot("coefficients originate from model with interaction effects while effect is specified as 'main'" = 
                !any(grepl(pattern="*", x=CIcoeffsWhite$met, fixed=TRUE)))
  }
  if (effect=="interaction") {
    stopifnot("coefficients originate from model with main effects while effect is specified as 'interaction'" = 
                any(grepl(pattern="*", x=CIcoeffsWhite$met, fixed=TRUE)))
  }
  
  # separate main and interaction regression coefficients
  CIcoeffsWhiteMain <- subset(CIcoeffsWhite, subset = met%in%c("Age", metabolites))
  CIcoeffsBlackMain <- subset(CIcoeffsBlack, subset = met%in%c("Age", metabolites))
  if (effect=="interaction") {
    CIcoeffsWhiteInt <- subset(CIcoeffsWhite, subset = ! met%in%c("Age", metabolites))
    CIcoeffsBlackInt <- subset(CIcoeffsBlack, subset = ! met%in%c("Age", metabolites))
  }
  
  # rearrange beta-coefficients for the main effects plot
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[order(abs(CIcoeffsWhiteMain$est.), decreasing=TRUE),]
  CIcoeffsBlackMain <- CIcoeffsBlackMain[match(CIcoeffsWhiteMain$met, CIcoeffsBlackMain$met), ]
  w <- order(abs(CIcoeffsWhiteMain$est.)+abs(CIcoeffsBlackMain$est.), 
             decreasing=TRUE)   # from large effects to small effects
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[w,]
  CIcoeffsBlackMain <- CIcoeffsBlackMain[w,]
  cluster1 <- c("met_013", "met_028", "met_031")
  cluster2 <- c("met_066", "met_071", "met_132", "met_133", "met_134")
  cluster3 <- c("met_003", "met_005", "met_011", "met_015", "met_018", "met_019", 
                "met_020", "met_027", "met_029", "met_034", "met_037", "met_049", 
                "met_050", "met_064", "met_073", "met_084", "met_114")
  cluster4 <- c("met_059", "met_060")
  cluster012 <- c("met_012", "met_012/met_038", "met_012/met_041", "met_078/met_012", 
                  "met_010/met_012", "met_012/met_026")
  cluster017 <- c("met_017", "met_010/met_017", "met_078/met_017")
  cluster026 <- c("met_026", "met_026/met_038", "met_026/met_041", "met_026/met_047", 
                "met_026/met_093")
  cluster032 <- c("met_032", "met_038/met_032", "met_041/met_032")
  cluster047 <- c("met_047", "met_010/met_047", "met_038/met_047", "met_041/met_047", 
                  "met_078/met_047")
  cluster093 <- c("met_093", "met_010/met_093", "met_038/met_093", "met_041/met_093", 
                  "met_078/met_093")
  
  
  w <- match(c("Age", cluster1, cluster2, cluster3, cluster4, cluster012, cluster017, 
               cluster026, cluster032, cluster047, cluster093), 
             CIcoeffsWhiteMain$met) # clusters first
  w <- w[which(!is.na(w))]
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[c(w, (1:nrow(CIcoeffsWhiteMain))[-w]),]
  CIcoeffsWhiteMain$order <- 1:nrow(CIcoeffsWhiteMain)
  CIcoeffsBlackMain <- CIcoeffsBlackMain %>% 
    mutate(order=match(met, CIcoeffsWhiteMain$met))
  
  CIcoeffsMain <- bind_rows(CIcoeffsWhiteMain, CIcoeffsBlackMain)
  CIcoeffsMain$Race <- factor(CIcoeffsMain$Race, levels=c("White", "Black"))
  
  assignCluster <- function(met) {
    cluster <- ifelse(met %in% cluster1, "A", 
                        ifelse(met %in% cluster2, "B", 
                        ifelse(met %in% cluster3, "C", 
                        ifelse(met %in% cluster4, "D", 
                        ifelse(met %in% cluster012, "met_012", 
                        ifelse(met %in% cluster017, "met_017", 
                        ifelse(met %in% cluster026, "met_026", 
                        ifelse(met %in% cluster032, "met_032", 
                        ifelse(met %in% cluster047, "met_047", 
                        ifelse(met %in% cluster093, "met_093", 
                               "not clustered"))))))))))
    return(cluster)
  }
  CIcoeffsMain <- CIcoeffsMain %>% mutate(cluster=assignCluster(met))
  
  # store plot of main effects
  title <- "Standardized coefficients + 95% CI"
  
  p_CIcoeffsMain <- ggplot(CIcoeffsMain,aes(x=reorder(met,-order), y=est., 
                                            ymin=`95pc lCI`, ymax=`95pc uCI`, 
                                            by=Race, col=cluster)) +
    geom_pointrange(stat="identity", aes(pch=Race), 
                    position=position_dodge2(width=0.5, padding=0.5)) +
    #geom_errorbar(stat="identity",
    #              position=position_dodge2(width=0.25, padding=0.5)) +
    geom_abline(slope=0, intercept=0, lty="dashed") +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab("") + ylab("corrected effect") +
    coord_flip() +
    labs(title=title, subtitle=paste0("Effect: ", effect, ", type: ", type, 
                                      ", transformation: ", transformation))
  
  # rearrange beta-coefficients for the interaction effects plot
  if (effect=="interaction") {
    CIcoeffsWhiteInt <- CIcoeffsWhiteInt[order(abs(CIcoeffsWhiteInt$est.), decreasing=TRUE),]
    CIcoeffsBlackInt <- CIcoeffsBlackInt[match(CIcoeffsWhiteInt$met, CIcoeffsBlackInt$met), ]
    w <- order(abs(CIcoeffsWhiteInt$est.)+abs(CIcoeffsBlackInt$est.), 
               decreasing=TRUE)   # from large effects to small effects
    CIcoeffsWhiteInt <- CIcoeffsWhiteInt[w,]
    CIcoeffsBlackInt <- CIcoeffsBlackInt[w,]
    
    CIcoeffsWhiteInt$order <- 1:nrow(CIcoeffsWhiteInt)
    CIcoeffsBlackInt$order <- 1:nrow(CIcoeffsBlackInt)
    
    CIcoeffsInt <- bind_rows(CIcoeffsWhiteInt, CIcoeffsBlackInt)
    CIcoeffsInt <- subset(CIcoeffsInt, subset = order<=50)
    CIcoeffsInt$Race <- factor(CIcoeffsInt$Race, levels=c("White", "Black"))
    
    clusterInteraction <- function(metabolite) {
      mets <- strsplit(metabolite, "*", fixed=TRUE)[[1]]
      ifelse(any(mets=="Age"), return("Age"), 
         ifelse(any(mets=="met_068"), return("met_068"), 
            ifelse(any(mets=="met_084"), return("met_084"), 
               {c1 <- ifelse(mets[1] %in% cluster1, "A", 
                            ifelse(mets[1] %in% cluster2, "B", 
                                   ifelse(mets[1] %in% cluster3, "C", 
                                          "X")));
               c2 <- ifelse(mets[2] %in% cluster1, "A", 
                            ifelse(mets[2] %in% cluster2, "B", 
                                   ifelse(mets[2] %in% cluster3, "C", "X")));
               intercluster <- paste0(c(c1, c2)[order(c(c1, c2))], collapse="*");
               return(intercluster)})))
      
    }
    
    # assign category to interaction effects
    CIcoeffsInt <- CIcoeffsInt %>% rowwise() %>% mutate(cluster=clusterInteraction(met))
    
    # assign order of interaction effects within assigned category
    CIcoeffsInt <- CIcoeffsInt %>% group_by(cluster, Race) %>% 
      mutate(groupOrder=sort.int(order, decreasing=TRUE, index.return=TRUE)$ix)
    
    # calculate size and cumulative size of categories
    cumGroupSize <- CIcoeffsInt %>% group_by(cluster) %>%
      summarise(sizes=n()) %>% mutate(cumSizes=cumsum(sizes))
    
    # assign order of metabolites on plot
    CIcoeffsInt <- merge(CIcoeffsInt, cumGroupSize, by="cluster") %>%
      mutate(order=groupOrder+cumSizes)
  }
  
  # store plot of interaction effects
  if (effect=="interaction") {
    p_CIcoeffsInt <- ggplot(CIcoeffsInt, aes(x=reorder(met,-order), y=est., 
                                              ymin=`95pc lCI`, ymax=`95pc uCI`, 
                                              by=Race, col=cluster)) +
      geom_pointrange(stat="identity", aes(pch=Race)) +
      geom_abline(slope=0, intercept=0, lty="dashed") +
      theme(axis.text.x = element_text(angle = 90)) +
      xlab("") + ylab("corrected effect") +
      coord_flip() +
      labs(title=title, subtitle=paste0("Effect: ", effect, ", type: ", type, 
                                        ", transformation: ", transformation))
  }
  
  # return plot(s)
  if (effect=="main") {
    return(p_CIcoeffsMain)
  }
  if (effect=="interaction") {
    return(list(main=p_CIcoeffsMain, interaction=p_CIcoeffsInt))
  }
  
}


