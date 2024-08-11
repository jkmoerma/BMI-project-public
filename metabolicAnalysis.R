
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
  
  # calculate correlations +
  # translate coefficients to a measure of relatedness to obesity
  effects <- coeffs
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
                  "met_026/met_093", "met_012/met_026")
  cluster032 <- c("met_032", "met_038/met_032", "met_041/met_032")
  cluster047 <- c("met_047", "met_010/met_047", "met_038/met_047", "met_041/met_047", 
                  "met_078/met_047")
  cluster093 <- c("met_093", "met_010/met_093", "met_038/met_093", "met_041/met_093", 
                  "met_078/met_093")
  
  summedEffects <- function(coeffs, cluster) {
    effects <- coeffs
    if (length(cluster)==2) {
      rho <- cor(data_ref[[cluster[1]]], data_ref[[cluster[2]]])
      corr <- data.frame(rowname = cluster,
                         met1 = c(1, rho),
                         met2 = c(rho, 1))
      names(corr) <- c("rowname", cluster)
    }
    if (length(cluster)>2) {
      corr <- cor_mat(data_ref[,cluster])
    }
    for (i in 1:length(cluster)) {
      met <- corr$rowname[i]
      partial_effects <- lapply(X=cluster, 
                                FUN = function(met) corr[[i,met]]*coeffs[[met]])
      sum_effects <- vector(mode="numeric", length=nrow(coeffs))
      for (j in 1:length(cluster)) {
        sum_effects <- sum_effects + partial_effects[[j]]
      }
      effects[[met]] <- sum_effects
    }
    return(effects)
  }
  
  effects <- summedEffects(effects, cluster1)
  effects <- summedEffects(effects, cluster2)
  effects <- summedEffects(effects, cluster3)
  effects <- summedEffects(effects, cluster4)
  effects <- summedEffects(effects, cluster012)
  effects <- summedEffects(effects, cluster017)
  effects <- summedEffects(effects, cluster026)
  effects <- summedEffects(effects, cluster032)
  effects <- summedEffects(effects, cluster047)
  effects <- summedEffects(effects, cluster093)
  
  # return effects
  return(effects)
  
}

plotScaledEffects <- function(bootstrapCoeffsWhite, bootstrapCoeffsBlack, 
                              bootstrapCoeffsAll, effect, type, transformation) {
  
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
  
  # generate coefficient estimates for main effects unstratified
  CIcoeffsAll <- bootstrapCoeffsAll %>% 
    group_by(met) %>% summarise("est."=quantile(x=coeff, probs=0.5),
                                "95pc lCI"=quantile(x=coeff, probs=0.025),
                                "95pc uCI"=quantile(x=coeff, probs=0.975))
  CIcoeffsAll["Race"] <- "All"
  
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
  CIcoeffsAllMain <- subset(CIcoeffsAll, subset = met%in%c("Age", metabolites))
  if (effect=="interaction") {
    CIcoeffsWhiteInt <- subset(CIcoeffsWhite, subset = ! met%in%c("Age", metabolites))
    CIcoeffsBlackInt <- subset(CIcoeffsBlack, subset = ! met%in%c("Age", metabolites))
    CIcoeffsAllInt <- subset(CIcoeffsAll, subset = ! met%in%c("Age", metabolites))
  }
  
  # rearrange beta-coefficients for the main effects plot
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[order(abs(CIcoeffsWhiteMain$est.), decreasing=TRUE),]
  CIcoeffsBlackMain <- CIcoeffsBlackMain[match(CIcoeffsWhiteMain$met, CIcoeffsBlackMain$met), ]
  CIcoeffsAllMain <- CIcoeffsAllMain[match(CIcoeffsWhiteMain$met, CIcoeffsAllMain$met), ]
  w <- order(abs(CIcoeffsWhiteMain$est.)+abs(CIcoeffsBlackMain$est.), 
             decreasing=TRUE)   # from large effects to small effects
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[w,]
  CIcoeffsBlackMain <- CIcoeffsBlackMain[w,]
  CIcoeffsAllMain <- CIcoeffsAllMain[w,]
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
  CIcoeffsWhiteMain <- CIcoeffsWhiteMain[c(w, (1:nrow(CIcoeffsWhiteMain))[-w]),]
  CIcoeffsWhiteMain$order <- 1:nrow(CIcoeffsWhiteMain)
  CIcoeffsBlackMain <- CIcoeffsBlackMain %>% 
    mutate(order=match(met, CIcoeffsWhiteMain$met))
  CIcoeffsAllMain <- CIcoeffsAllMain %>% 
    mutate(order=match(met, CIcoeffsWhiteMain$met))
  
  CIcoeffsMain <- bind_rows(CIcoeffsWhiteMain, CIcoeffsBlackMain, CIcoeffsAllMain)
  #CIcoeffsMain <- subset(CIcoeffsMain, subset = order<=50)
  
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
  title <- "Correlation-corrected effects + 95% CI"
  
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

plotANOVA <- function(data, ethnicity) {
  for (met in metabolites) {
    data[[met]] <- scale(data[[met]])
  }
  groupLevels <- subset(data, subset=!is.na(predictionGroup), 
                             select=c(metabolites, "predictionGroup"))
  
  levelDiffsNormal <- matrix(nrow=length(metabolites), ncol=6)
  rownames(levelDiffsNormal) <- metabolites
  colnames(levelDiffsNormal) <- c("NO-NN", "NO-NN (lCI)", "NO-NN (uCI)",
                                       "OO-NO", "OO-NO (lCI)", "OO-NO (uCI)")
  
  levelDiffsObese <- matrix(nrow=length(metabolites), ncol=6)
  rownames(levelDiffsObese) <- metabolites
  colnames(levelDiffsObese) <- c("OO-ON", "OO-ON (lCI)", "OO-ON (uCI)",
                                      "NN-ON", "NN-ON (lCI)", "NN-ON (uCI)")
  
  legend <- c()
  for (met in metabolites) {
    diffs <- TukeyHSD(aov(groupLevels[[met]]~groupLevels$predictionGroup))$`groupLevels$predictionGroup`
    levelDiffsNormal[met,] <- c(diffs["NO-NN", "diff"], diffs["NO-NN", "lwr"], 
                                     diffs["NO-NN", "upr"], diffs["OO-NO", "diff"], 
                                     diffs["OO-NO", "lwr"], diffs["OO-NO", "upr"])
    levelDiffsObese[met,] <- c(diffs["OO-ON", "diff"], diffs["OO-ON", "lwr"], 
                                    diffs["OO-ON", "upr"], -diffs["ON-NN", "diff"], 
                                    -diffs["ON-NN", "upr"], -diffs["ON-NN", "lwr"])
    if (abs(diffs["OO-NN", "diff"])>0.5) {
      legend <- c(legend, met)
    } else {
      legend <- c(legend, "other")
    }
  }
  
  levelDiffsNormal <- data.frame(met=rownames(levelDiffsNormal),
                                      legend = legend,
                                      as.data.frame(levelDiffsNormal),
                                      check.names=FALSE)
  levelDiffsObese <- data.frame(met=rownames(levelDiffsObese),
                                     legend=legend,
                                     as.data.frame(levelDiffsObese),
                                     check.names=FALSE)
  
  pN <- ggplot(data=levelDiffsNormal, 
               aes(x=`OO-NO`, y=`NO-NN`, label=met, col=met)) + 
          theme(legend.position = "none") +
          coord_fixed(ratio=1) +
          geom_abline(slope=0, intercept=0, lty="dashed") +
          geom_abline(slope=-1, intercept=0, lty="dashed") +
          geom_vline(xintercept=0, lty="dashed") +
          geom_pointrange(aes(ymin=`NO-NN (lCI)`, ymax=`NO-NN (uCI)`), fatten=1, lwd=0.1) +
          geom_pointrange(aes(xmin=`OO-NO (lCI)`, xmax=`OO-NO (uCI)`), fatten=1, lwd=0.1) +
          geom_text(vjust=0, nudge_y=0.05, size=2) +
          labs(title=paste0("NO (", ethnicity, " ethnicity)"), 
               x="Scaled log conc. diff. OO-NO", 
               y="Scaled log conc. diff. NO-NN")
  
  pO <- ggplot(data=levelDiffsObese, 
               aes(x=`NN-ON`, y=`OO-ON`, label=met, col=met)) + 
          theme(legend.position = "none") +
          coord_fixed(ratio=1) + 
          geom_abline(slope=0, intercept=0, lty="dashed") +
          geom_abline(slope=1, intercept=0, lty="dashed") +
          geom_vline(xintercept=0, lty="dashed") +
          geom_pointrange(aes(ymin=`OO-ON (lCI)`, ymax=`OO-ON (uCI)`), fatten=1, lwd=0.1) +
          geom_pointrange(aes(xmin=`NN-ON (lCI)`, xmax=`NN-ON (uCI)`), fatten=1, lwd=0.1) +
          geom_text(vjust=0, nudge_y=0.05, size=2) +
          labs(title=paste0("ON (", ethnicity, " ethnicity)"), 
               x="Scaled log conc. diff. NN-ON", 
               y="Scaled log conc. diff. OO-ON")
  
  list(normal=pN, obese=pO)
}
