#' Make a table of patient counts for every ethnicity inside every obesity class
#'
#' @param df A data frame containing columns "ObesityClass" and "Race".
#' @return A matrix with all the observed ethnicities as columns, the different obesity classes and the ethnicity total as rows and the counts as entries.
#' @examples
#' 
tableBiometrics <- function(df) {
  count_table <- df %>% group_by(ObesityClass, Race) %>% summarise(n())
  biometrics <- spread(count_table, key = Race, value = `n()`)
  
  count_all <- df %>% group_by(Race) %>% summarise(n())
  biometrics1 <- spread(count_all, key = Race, value = `n()`)
  
  answer_tibble <- bind_rows(biometrics, biometrics1)
  answer <- as.matrix(answer_tibble[,-1])
  rownames(answer) <- c(as.character(answer_tibble$ObesityClass[-nrow(answer_tibble)]), "total")
  answer
}

#' Make a table of patient counts of smokers and non-smokers for every ethnicity
#'
#' @param df A data frame containing columns "Smoking status" and "Race".
#' @return A matrix with all the observed ethnicities as columns, counted "Smoking" and "Non-smoking" patients as entries.
#' @examples
#' 
tableSmokingStatus <- function(df) {
  count_table <- df %>% group_by(`Smoking status`, Race) %>% summarise(n())
  SmokingStatus <- spread(count_table, key = Race, value = `n()`, fill=0)
  answer <- as.matrix(SmokingStatus[,-1])
  answer <- cbind(answer, "all"=margin.table(answer, margin=1))
  rownames(answer) <- c("Non-smoking", "Smoking")[1+(SmokingStatus$`Smoking status`)]
  answer
}


#' Assess log-normality of metabolite distributions + tagging outliers
#'
#' @param df A data frame containing columns "ID", "Race", "BMI" and metabolite levels.
#' @return A list object containing a data frame with outliers for every metabolite, a count of outliers + frequency of common outliers with other metabolites, a table of metabolites which were not log-normally distributed + p-value of the Shapiro-Wilk test
#' @examples
#' 
logNormalityAndOutliers <- function(df) {
  
  outliers <- data.frame(ID=df$ID, Race=df$Race)
  logNormality <- vector(mode="numeric", length=length(metabolites)+1)
  names(logNormality) <- c("BMI", metabolites)
  
  for (met in c("BMI", metabolites)) {
    
    # for standardizing log concentrations of metabolites
    standardize <- paste0("(log(", met, ")-mean(log(", met, 
                               "), na.rm=TRUE))/sd(log(", met, "), na.rm=TRUE)")
    
    # for calculating theoretical quantiles under log-normality
    Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")
    
    # for retrieving the smooth distribution value of all observations
    quantSmooth <- "smoother <- loess(formula=standardized~normQuantile, span=0.4, 
                                     degree=2, family='symmetric');
                    predict(smoother, normQuantile)"
    
    # applying the above defined formulas on the data
    investigation <- df %>% group_by(Race) %>% 
      mutate(standardized=eval(parse(text=standardize)),
             normQuantile=eval(parse(text=Zscore)),
             quantSmooth=eval(parse(text=quantSmooth)),
             outliers=(standardized-quantSmooth)*normQuantile/abs(normQuantile)>0.7)
    
    # store detected outliers
    outliers <- merge(outliers, investigation[,c("ID", "outliers")], by="ID")
    outliers[[met]] <- outliers$outliers
    outliers$outliers <- NULL
    
    # test normality
    logNormality[met] <- shapiro.test(investigation$standardized[which(!outliers[[met]])])$p.value
    
  }
  
  nonNormals <- logNormality[which(logNormality<1e-4)]
  totable <- cbind("met."=names(nonNormals), "p-val."=signif(nonNormals, digits=2))
  
  # create table to investigate dependencies in outliers between metabolites
  countOutliers <- data.frame(met=metabolites, 
                              count=vector(mode="numeric", length=length(metabolites)),
                              commons=vector(mode="character", length=length(metabolites)),
                              BMIperc=vector(mode="character", length=length(metabolites)))
  BMIpercentages <- df %>% group_by(Race) %>% mutate(percent=percent_rank(BMI)*100)
  BMIpercentages <- BMIpercentages$percent
  for (i in 1:length(metabolites)) {
    met1 <- countOutliers$met[i]
    countOutliers$count[i] <- sum(outliers[[met1]])
    if (sum(outliers[[met1]], na.rm=TRUE)>0) {
      countOutliers$BMIperc[i] <- paste(paste0(round(BMIpercentages, digits=1)[which(outliers[[met1]])], "\\%"), collapse=" ")
    }
    commonSet <- c()
    for (j in (1:length(metabolites))[-i]) {
      met2 <- countOutliers$met[j]
      if (length(which(outliers[[met1]]&outliers[[met2]]))>0) {
        commonSet <- c(commonSet, paste0(met2,":", length(which(outliers[[met1]]&outliers[[met2]]))))
      }
    }
    if (length(commonSet) > 0) {
      countOutliers$commons[i] <- paste0(commonSet, collapse=" ")
    }
  }
  
  list(outliers=outliers, countOutliers=countOutliers, nonNormals=totable)
  
}

#' Makes QQ-plots of a subselection of metabolites, stored in the file "QQplots.pdf"
#'
#' @param df A data frame containing columns "ID", "BMI", "Race" and metabolite levels.
#' @param mets A subselection of metabolites to make a QQ-plot of
#' @return 
#' @examples
#' 
makeQQplot <- function(df, mets, pvals) {
  
  outliers <- data.frame(ID=df$ID, Race=df$Race)
  logNormality <- vector(mode="numeric", length=length(metabolites)+1)
  names(logNormality) <- c("BMI", metabolites)
  colors <- c("red","blue","green","violet","black")
  chars <- c("W", "B", "S", "E", "M")
  
  pdf(file="QQplots.pdf", width=8.27, height=11.69)
  layout(matrix(1:6, byrow=TRUE, ncol=2))
  par(mai=c(0.6, 0.6, 0.4, 0.2), mex=1)
  for (i in 1:length(mets)) {
    
    met <- mets[i]
    pval <- pvals[i]
    
    # for standardizing log concentrations of metabolites
    standardize <- paste0("(log(", met, ")-mean(log(", met, 
                          "), na.rm=TRUE))/sd(log(", met, "), na.rm=TRUE)")
    
    # for calculating theoretical quantiles under log-normality
    Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")
    
    # for retrieving the smooth distribution value of all observations
    quantSmooth <- "smoother <- loess(formula=standardized~normQuantile, span=0.4, 
                                     degree=2, family='symmetric');
                    predict(smoother, normQuantile)"
    
    # applying the above defined formulas on the data
    investigation <- df %>% group_by(Race) %>% 
      mutate(standardized=eval(parse(text=standardize)),
             normQuantile=eval(parse(text=Zscore)),
             quantSmooth=eval(parse(text=quantSmooth)),
             outliers=(standardized-quantSmooth)*normQuantile/abs(normQuantile)>0.7)
    
    # plot QQ-plot of metabolite distribution
    normRangeX <- range(investigation$normQuantile, na.rm=TRUE)
    normRangeY <- range(investigation$standardized, na.rm=TRUE)
    plot(x=normRangeX, y=normRangeY, col="white", main=met,
         xlab="theoretical quantiles", ylab="standardized log concentrations")
    abline(a=0, b=1, lty="dashed", col="darkgrey")
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      normPlot <- subset(investigation, subset = Race==race)
      plotchar <- rep(chars[j], times=nrow(normPlot))
      w <- which(normPlot$outliers)
      plotchar[w] <- "."
      points(x=normPlot$normQuantile, y=normPlot$standardized, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$normQuantile[w], y=normPlot$standardized[w], col=colors[j], 
             labels=normPlot$ID[w], cex=0.7)
      }
    }
    
    # result normality test
    text(x=normRangeX[1], y=normRangeY[2], col="red", adj=c(0,1),
         labels=paste0("Shap.-Wilk p-val.: %.1e", pval))
         #labels=sprintf("Shap.-Wilk p-val.: %.1e", pval))
    
  }
  dev.off()
  
  
}



#' Stores a plot to demonstrate outlier detection with metabolites and ethnicities of interest, stored in the file "outlierDemonstration.pdf"
#'
#' @param df A data frame containing columns "ID", "BMI", "Race" and metabolite levels.
#' @param mets Metabolite names of interest
#' @param ethnicity Ethnicities according to the metabolites of interest
#' @return 
#' @examples
#' 
demonstrateOutlierDetection <- function(df, mets, ethnicity) {
  
  stopifnot("arguments 'mets' and 'ethnicity' must have the same length" = 
              length(mets)==length(ethnicity))
  
  outliers <- data.frame(ID=df$ID, Race=df$Race)
  colors <- c("red","blue","green","violet","black")
  fillRegular <- c("indianred","cadetblue","lightgreen","lightpink","grey")
  chars <- c("W", "B", "S", "E", "M")
  if (length(mets)%%2 == 0) {
    height <- 11.69 * length(mets)/6
    plot_pos <- 1:length(mets)
  }
  if (length(mets)%%2 == 1) {
    height <- 11.69 * (length(mets)+1)/6
    plot_pos <- c(1:length(mets), 0)
  }
  pdf(file="outlierDemonstration.pdf", width=8.27, height=height)
  layout(matrix(plot_pos, byrow=TRUE, ncol=2))
  par(mai=c(0.6, 0.6, 0.4, 0.2), mex=0.5)
  for (i in 1:length(mets)) {
    met <- mets[i]
    # for standardizing log concentrations of metabolites
    standardize <- paste0("(log(", met, ")-mean(log(", met, 
                          "), na.rm=TRUE))/sd(log(", met, "), na.rm=TRUE)")
    
    # for calculating theoretical quantiles under log-normality
    Zscore <- paste0("qqnorm(", standardize, ", plot.it=FALSE)$x")
    
    # for retrieving the smooth distribution value of all observations
    quantSmooth <- "smoother <- loess(formula=standardized~normQuantile, span=0.4, 
                                     degree=2, family='symmetric');
                    predict(smoother, normQuantile)"
    
    # applying the above defined formulas on the data
    investigation <- df %>% group_by(Race) %>% 
      mutate(standardized=eval(parse(text=standardize)),
             normQuantile=eval(parse(text=Zscore)),
             quantSmooth=eval(parse(text=quantSmooth)),
             outliers=(standardized-quantSmooth)*normQuantile/abs(normQuantile)>0.7)
    
    # plot QQ-plot of metabolite distribution
    normRangeX <- range(investigation$normQuantile, na.rm=TRUE)
    normRangeY <- range(investigation$standardized, na.rm=TRUE)
    plot(x=normRangeX, y=normRangeY, col="white", main=met,
         xlab="theoretical quantiles", ylab="standardized log concentrations")
    abline(a=0, b=1, lty="dashed", col="darkgrey")
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      normPlot <- subset(investigation, subset = Race==race)
      plotchar <- rep(chars[j], times=nrow(normPlot))
      w <- which(normPlot$outliers)
      plotchar[w] <- "."
      if (race == ethnicity[i]) {
        lines(x=sort(normPlot$normQuantile), y=sort(normPlot$quantSmooth), col=colors[j])
        polygon(x=c(sort(normPlot$normQuantile, decreasing=FALSE), sort(normPlot$normQuantile, decreasing=TRUE)), 
                y=c(sort(normPlot$quantSmooth, decreasing=FALSE)-0.7, 
                    sort(normPlot$quantSmooth, decreasing=TRUE)+0.7), 
                col=fillRegular[j], border="white")
      }
      points(x=normPlot$normQuantile, y=normPlot$standardized, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$normQuantile[w], y=normPlot$standardized[w], col=colors[j], 
             labels=normPlot$ID[w], cex=0.7)
      }
    }
    
  }
  dev.off()
}

#' This function was written to impute and explore missing data. It uses Multivariate Imputation Chained Equations to perform this.
#'
#' @param df A data frame as given by the owner of the data
#' @param makePlot logical, if TRUE, a plot will be generated and stored in file "boxplotMissing.pdf"
#' @return A similar data frame with MICE imputed missing values for met_002, met_010 and met_068.
#' 
boxplotMissing <- function(df, makePlot=TRUE) {
  
  # one missing value of met_010 will be imputed with the ethnicity stratified mean
  # metabolite ratios with met_010 will be recalculated
  met10colnames <- names(df)[which(grepl(names(df), pattern="met_010"))]
  stratified_means_met010 <- df %>% group_by(Race) %>% summarise(exp(mean(log(met_010), na.rm=TRUE)))
  df_complete <- merge(df, stratified_means_met010, by="Race")
  w <- which(is.na(df_complete$met_010))
  df_complete$met_010[w] <- df_complete$`exp(mean(log(met_010), na.rm = TRUE))`[w]
  df_complete$`exp(mean(log(met_010), na.rm = TRUE))` <- NULL
  for (name in met10colnames) {
    df_complete[[name]][w] <- with(df_complete[w,], eval(parse(text=name)))
  }
  
  data_model <- df_complete
  data_model$Smoking <- df_complete$`Smoking status`
  data_model$`Smoking status` <- NULL
  data_model$`Maternal Age` <- NULL
  data_model$AgeGroup <- NULL
  data_model$ObesityClass <- NULL
  for (met in c("BMI", metabolites)) {
    data_model[[met]] <- log(df_complete[[met]])
  }
  
  # impute met_002 and met_068 stratified on ethnicity
  log_all <- c()
  for (race in unique(data_model$Race)) {
    data_race <- subset(data_model, subset=Race==race)
    if (nrow(data_race) > length(metabolites)+10) {
      imputed_data <- mice(data_race, method="norm.predict")
    } else {
      imputed_data <- mice(data_race, method="pmm")
    }
    
    log_race <- complete(imputed_data) %>% 
                  mutate(met_002C=met_002) %>% 
                  mutate(met_068C=met_068)
    log_race <- subset(log_race, select=c("ID", "met_002C", "met_068C"))
    log_all <- bind_rows(log_all, log_race)
    
  }
  
  df_complete <- merge(df_complete, log_all, by="ID")
  df_complete$met_002 <- exp(df_complete$met_002C)
  df_complete$met_068 <- exp(df_complete$met_068C)
  
  if (makePlot) {
  
    missing02 <- wilcox.test(formula = value ~ missing, 
                             data=data.frame(value = df_complete$met_002, 
                                             missing = df_complete$ID %in% df$ID[which(is.na(df$met_002))]))
    missing68 <- wilcox.test(formula = value ~ missing, 
                             data=data.frame(value=df_complete$met_068, 
                                             missing = df_complete$ID %in% df$ID[which(is.na(df$met_068))]))
    
    pdf(file="boxplotMissing.pdf", width=8.27, height=5.845)
    
    layout(matrix(1:2, ncol=2))
    
    title <- sprintf("met_002 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f", 
                     missing02$p.value, 1-missing02$statistic/(38*1597))
    boxplot(formula = value ~ missing, 
            data=data.frame(value=df_complete$met_002, 
                            missing = df_complete$ID %in% df$ID[which(is.na(df$met_002))]), 
            xlab="", ylab="log(met_002)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))
    
    title <- sprintf("met_068 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f", 
                     missing68$p.value, 1-missing68$statistic/(49*1586))
    boxplot(formula = value ~ missing, 
            data=data.frame(value=df_complete$met_068, 
                            missing = df_complete$ID %in% df$ID[which(is.na(df$met_068))]), 
            xlab="", ylab="log(met_068)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))
    
    dev.off()
  
  }
  
  w <- match(df$ID, df_complete$ID)
  df_complete[w,]
  
}

#' This function splits the data set in outliers and regular observations
#'
#' @param df A data frame as given by the owner of the data with imputed values, it must contain a column with unique patient identifiers, called "ID".
#' @param outliers A data frame of logical values indicating whether values of a metabolite were flagged as outliers, it must contain a column with unique patient identifiers, called "ID".
#' @param mets The metabolites used to filter out outliers
#' @return A list containing the data frames of "regulars" and "outliers"
#' @examples 
#' df <- data.frame(ID=c("a", "b", "c", "d", "e"),
#'                  met1=c(2, 1, 10, 4, 3),
#'                  met2=c(20, 11, 12, 8, 9)
#'                  met3=c(0.5, 6, 0.1, 0.4, 1.1))
#' outliers <- data.frame(ID=c("a", "b", "c", "d", "e"),
#'                        met1=c(FALSE, FALSE, TRUE, FALSE, FALSE),
#'                        met2=c(TRUE, FALSE, FALSE, FALSE, FALSE)
#'                        met3=c(FALSE, TRUE, FALSE, FALSE, FALSE))
#' filterOutliers(df, outliers, mets=c("met1", "met2"))
#' filterOutliers(df, outliers, mets=c("met1", "met3"))
#' 
filterOutliers <- function(df, outliers, mets) {
  subsetOutlier <- with(outliers, eval(parse(text=paste(mets, collapse="|"))))
  outlierIDs <- outliers$ID[which(subsetOutlier)]
  w <- which(df$ID %in% outlierIDs)
  list(regulars=df[-w,], outliers=df[w,])
}





