#' Make a table of patient counts for every ethnicity inside every obesity class
#'
#' @param df A data frame containing columns "ObesityClass" and "Race".
#' @return A matrix with all the observed ethnicities as columns, the different obesity classes and the ethnicity total as rows and the counts as entries.
#' @examples
#' add(1, 1)
#' add(10, 1)
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


tableQuartiles <- function(df) {
  #df %>% group_by(Race) %>% summarise(IQR(met_002, na.rm=TRUE), IQR(met_003))
  
  quartileTable <- matrix(ncol=length(unique(df$Race))+1, nrow=length(metabolites))
  rownames(quartileTable) <- metabolites
  colnames(quartileTable) <- c(unique(df$Race), "ANOVA p-val.")
  for (met in metabolites) {
    quartileTable[met, "ANOVA p-val."] <- 
      anova_test(formula=eval(expr=parse(text=paste0("log(", met, ") ~ Race*BMI"))), 
                 data=df, type=3)$p[1]
    for (race in unique(df$Race)) {
      df_s <- subset(df, subset = Race==race)
      quartileTable[met, race] <- sprintf("%.02e [%.02e %.02e]", 
                                          log(quantile(p=0.5, df_s[[met]], na.rm=TRUE)),
                                          log(quantile(p=0.25, df_s[[met]], na.rm=TRUE)),
                                          log(quantile(p=0.75, df_s[[met]], na.rm=TRUE)))
    }
  }
}

#' Assess log-normality of metabolite distributions + tagging outliers
#'
#' @param df A data frame containing columns "ID", "BMI", "Race" and metabolite levels.
#' @return 
#' @examples
#' add(1, 1)
#' add(10, 1)
logNormalityAndOutliers <- function(df) {
  
  #df %>% group_by(Race) %>% 
  #  mutate(standardized=(log(met_002)-mean(log(met_002), na.rm=TRUE))/sd(log(met_002), na.rm=TRUE))
  
  outliers <- data.frame(ID=df$ID, Race=df$Race)
  logNormality <- vector(mode="numeric", length=length(metabolites)+1)
  names(logNormality) <- c("BMI", metabolites)
  colors <- c("red","blue","green","violet","black")
  chars <- c("W", "B", "S", "E", "M")
  pdf(file="QQplots.pdf", width=8.27, height=11.69)
  layout(matrix(1:6, byrow=TRUE, ncol=2))
  par(mai=c(0.6, 0.6, 0.4, 0.2), mex=1)
  for (met in c("BMI", metabolites)) {
    normalizedQuantiles <- list(x=c(), y=c(), ethnicity=c(), ID=c(), quantSmooth=c())
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      df_s <- subset(df, subset = Race==race)
      numbers <- log(df_s[[met]])
      standardNumbers <- (numbers-mean(numbers, na.rm=TRUE))/sd(numbers, na.rm=TRUE)
      quantileplot <- qqnorm(standardNumbers, plot.it=FALSE)
      normalizedQuantiles$x <- c(normalizedQuantiles$x, quantileplot$x)
      normalizedQuantiles$y <- c(normalizedQuantiles$y, quantileplot$y)
      normalizedQuantiles$ethnicity <- c(normalizedQuantiles$ethnicity, 
                                         rep(race, times=length(quantileplot$x)))
      normalizedQuantiles$ID <- c(normalizedQuantiles$ID, df_s$ID)
      
      smoother <- loess(formula=y~x, span=0.4, degree=2, family="symmetric",
                        data=data.frame(x=quantileplot$x, y=quantileplot$y))
      normalizedQuantiles$quantSmooth <- c(normalizedQuantiles$quantSmooth, 
                                           predict(smoother, quantileplot$x))
    }
    normalizedQuantiles <- as.data.frame(normalizedQuantiles)
    normalizedQuantiles$outliers <- with(normalizedQuantiles, (y-quantSmooth)*x/abs(x)>0.7)
    normRangeX <- range(normalizedQuantiles$x, na.rm=TRUE)
    normRangeY <- range(normalizedQuantiles$y, na.rm=TRUE)
    plot(x=normRangeX, y=normRangeY, col="white", main=met,
         xlab="theoretical quantiles", ylab="standardized log concentrations")
    abline(a=0, b=1, lty="dashed", col="darkgrey")
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      normPlot <- subset(normalizedQuantiles, subset = ethnicity==race)
      plotchar <- rep(chars[j], times=nrow(normPlot))
      w <- which(normPlot$outliers)
      plotchar[w] <- "."
      points(x=normPlot$x, y=normPlot$y, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$x[w], y=normPlot$y[w], col=colors[j], 
             labels=normPlot$ID[w], cex=0.7)
      }
    }
    
    # flag outliers
    outlierID <- subset(normalizedQuantiles, subset = outliers)$ID
    outliers[[met]] <- rep(FALSE, times=nrow(df))
    outliers[[met]][which(outliers$ID%in%outlierID)] <- TRUE
    
    # test normality
    logNormality[met] <- shapiro.test(normalizedQuantiles$y[which(!outliers[[met]])])$p.value
    text(x=normRangeX[1], y=normRangeY[2], col="red", adj=c(0,1),
         labels=sprintf("Shap.-Wilk p-val.: %.1e", logNormality[met]))
    
  }
  dev.off()
  
  nonNormals <- logNormality[which(logNormality<1e-4)]
  totable <- cbind("met."=names(nonNormals), "p-val."=signif(nonNormals, digits=2))
  
  countOutliers <- data.frame(met=metabolites, 
                              count=vector(mode="numeric", length=length(metabolites)),
                              commons=vector(mode="character", length=length(metabolites)),
                              BMIperc=vector(mode="character", length=length(metabolites)))
  BMIpercentages <- vector(mode="numeric", length=nrow(df))
  for (race in unique(df$Race)) {
    w <- which(df$Race==race)
    BMIpercentages[w] <- percent_rank(df$BMI[w])*100
  }
  for (i in 1:length(metabolites)) {
    met1 <- countOutliers$met[i]
    countOutliers$count[i] <- sum(outliers[[met1]])
    if (sum(outliers[[met1]])>0) {
      countOutliers$BMIperc[i] <- paste(paste0(round(BMIpercentages, digits=1)[which(outliers[[met1]])], "\\%"), collapse=" ")
    }
    commonSet <- c()
    for (j in (1:length(metabolites))[-i]) {
      met2 <- countOutliers$met[j]
      if (any(outliers[[met1]]&outliers[[met2]])) {
        commonSet <- c(commonSet, paste0(met2,":",sum(outliers[[met1]]&outliers[[met2]])))
      }
    }
    if (length(commonSet) > 0) {
      countOutliers$commons[i] <- paste0(commonSet, collapse=" ")
    }
  }
  
  list(outliers=outliers, countOutliers=countOutliers, nonNormals=totable)
  
}

# give a demonstration about ho outliers are detected
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
    normalizedQuantiles <- list(x=c(), y=c(), ethnicity=c(), ID=c(), quantSmooth=c())
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      df_s <- subset(df, subset = Race==race)
      numbers <- log(df_s[[met]])
      standardNumbers <- (numbers-mean(numbers, na.rm=TRUE))/sd(numbers, na.rm=TRUE)
      quantileplot <- qqnorm(standardNumbers, plot.it=FALSE)
      normalizedQuantiles$x <- c(normalizedQuantiles$x, quantileplot$x)
      normalizedQuantiles$y <- c(normalizedQuantiles$y, quantileplot$y)
      normalizedQuantiles$ethnicity <- c(normalizedQuantiles$ethnicity, 
                                         rep(race, times=length(quantileplot$x)))
      normalizedQuantiles$ID <- c(normalizedQuantiles$ID, df_s$ID)
      
      smoother <- loess(formula=y~x, span=0.4, degree=2, family="symmetric",
                        data=data.frame(x=quantileplot$x, y=quantileplot$y))
      normalizedQuantiles$quantSmooth <- c(normalizedQuantiles$quantSmooth, 
                                           predict(smoother, quantileplot$x))
    }
    normalizedQuantiles <- as.data.frame(normalizedQuantiles)
    normalizedQuantiles$outliers <- with(normalizedQuantiles, (y-quantSmooth)*x/abs(x)>0.7)
    normRangeX <- range(normalizedQuantiles$x, na.rm=TRUE)
    normRangeY <- range(normalizedQuantiles$y, na.rm=TRUE)
    plot(x=normRangeX, y=normRangeY, col="white", main=met,
         xlab="theoretical quantiles", ylab="standardized log concentrations")
    abline(a=0, b=1, lty="dashed", col="darkgrey")
    for (j in 1:length(ethnicities)) {
      race <- ethnicities[j]
      normPlot <- subset(normalizedQuantiles, subset = ethnicity==race)
      plotchar <- rep(chars[j], times=nrow(normPlot))
      w <- which(normPlot$outliers)
      plotchar[w] <- "."
      if (race == ethnicity[i]) {
        lines(x=sort(normPlot$x), y=sort(normPlot$quantSmooth), col=colors[j])
        polygon(x=c(sort(normPlot$x, decreasing=FALSE), sort(normPlot$x, decreasing=TRUE)), 
                y=c(sort(normPlot$quantSmooth, decreasing=FALSE)-0.7, 
                    sort(normPlot$quantSmooth, decreasing=TRUE)+0.7), 
                col=fillRegular[j], border="white")
      }
      points(x=normPlot$x, y=normPlot$y, col=colors[j], pch=plotchar)
      if (length(w)>0) {
        text(x=normPlot$x[w], y=normPlot$y[w], col=colors[j], 
             labels=normPlot$ID[w], cex=0.7)
      }
    }
    
  }
  dev.off()
}

# calculating measures for association between metabolites
tableCorrelations <- function(df) {
  
  pearsonCors <- matrix(nrow=length(metabolites), ncol=length(metabolites))
  rownames(pearsonCors) <- metabolites
  colnames(pearsonCors) <- metabolites
  pearsonCorsTable <- matrix(nrow=length(metabolites), ncol=length(metabolites))
  rownames(pearsonCorsTable) <- metabolites
  colnames(pearsonCorsTable) <- metabolites
  pearsonSignif <- matrix(nrow=length(metabolites), ncol=length(metabolites))
  rownames(pearsonSignif) <- metabolites
  colnames(pearsonSignif) <- metabolites
  for (i in 1:(length(metabolites)-1)) {
    for (j in (i+1):length(metabolites)) {
      met1 <- metabolites[i]
      met2 <- metabolites[j]
      corr <- cor.test(formula = eval(expr=parse(text=paste0("~ log(", met1, ") + log(", met2, ")"))), 
                       data=df, method="pearson", 
                       subset=eval(expr=parse(text=paste0("!is.na(", met1, ") & !is.na(", met2, ")"))))
      pearsonCors[met1,met2] <- 
        sprintf("%.3f [%.3f %.3f]", corr$estimate, corr$conf.int[1], corr$conf.int[2]) 
      pearsonCorsTable[met1,met2] <- corr$estimate
      pearsonSignif[met1,met2] <- corr$p.value
      pearsonCors[met2,met1] <- 
        sprintf("%.3f [%.3f %.3f]", corr$estimate, corr$conf.int[1], corr$conf.int[2]) 
      pearsonCorsTable[met2,met1] <- corr$estimate
      pearsonSignif[met2,met1] <- corr$p.value
      
    }
  }
  
  selectCors <- which(abs(pearsonCorsTable)>0.8)
  colNumbers <- selectCors%%length(metabolites)
  colNumbers[which(colNumbers==0)] <- length(metabolites)
  rowNumbers <- 1 + floor((selectCors-1)/length(metabolites))
  met1 <- metabolites[rowNumbers]
  met2 <- metabolites[colNumbers]
  
  write.table(x=round(pearsonCorsTable, digits=3), file="correlations.csv", sep=",")
  
  listCors <- vector(mode="numeric", length=length(met1))
  for (i in 1:length(met1)) {
    listCors[i] <- pearsonCorsTable[met2[i], met1[i]]
  }
  
  cbind(met1, met2, "Pearson Cor."=round(listCors, digits=3))

}

# missing data
boxplotMissing <- function(df, makePlot=TRUE) {
  
  # one missing value of met_010 will be imputed with the ethnicity stratified mean
  # metabolite ratios with met_010 will be recalculated
  stratified_means_met010 <- df %>% group_by(Race) %>% summarise(exp(mean(log(met_010))))
  df <- merge(df, stratified_means_met010, by="Race")
  met10colnames <- names(df)[which(grepl(names(df), pattern="met_010"))]
  w <- which(is.na(df$met_010))
  df$met_010[w] <- df$`exp(mean(log(met_010)))`[w]
  df$`exp(mean(log(met_010)))` <- NULL
  for (name in met10colnames) {
    df[[name]][w] <- with(df[w,], eval(parse(text=name)))
  }
  
  data_model <- as.data.frame(log(makeX(train=df[,c("BMI", metabolites)], 
                                        na.impute=FALSE)),
                              check.names=FALSE)
  data_model <- data.frame(Race=relevel(as.factor(df$Race), ref="White"), 
                           #Smoking=relevel(as.factor(df$`Smoking status`), ref="FALSE"), 
                           Smoking=as.numeric(df$`Smoking status`), 
                           Age=(df$`Maternal Age`-mean(df$`Maternal Age`))/sd(df$`Maternal Age`),
                           ObesityClass=relevel(as.factor(df$ObesityClass), ref="Normal weight"),
                           data_model, check.names=FALSE)
  colnames(data_model) <- c("Race", "Smoking", "Age", "ObesityClass", "BMI", metabolites)
  
  # impute met_002 and met_068 stratified on ethnicity
  for (race in unique(df$Race)) {
    # subset first ...
    log_complete <- complete(mice(data_model, method="lasso.norm"))
  }
  
  df_complete <- 
  
  if (makePlot) {
  
    missing02 <- wilcox.test(formula = value ~ missing, 
                             data=data.frame(value=df_complete$met_002, missing=is.na(df$met_002)))
    missing68 <- wilcox.test(formula = value ~ missing, 
                             data=data.frame(value=df_complete$met_068, missing=is.na(df$met_068)))
    
    pdf(file="boxplotMissing.pdf", width=8.27, height=5.845)
    
    layout(matrix(1:2, ncol=2))
    
    title <- sprintf("met_002 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f", 
                     missing02$p.value, 1-missing02$statistic/(38*1597))
    boxplot(formula = value ~ missing, 
            data=data.frame(value=log(df_complete$met_002), missing=is.na(df$met_002)), 
            xlab="", ylab="log(met_002)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))
    
    title <- sprintf("met_068 \n Wilcox. p-val. %.3f \n P(met[miss] > met) = %.3f", 
                     missing68$p.value, 1-missing68$statistic/(49*1586))
    boxplot(formula = value ~ missing, 
            data=data.frame(value=log(df_complete$met_068), missing=is.na(df$met_068)), 
            xlab="", ylab="log(met_068)", main=title, xaxt="n")
    axis(side=1, at=1:2, labels=c("measured", "missing"))
    
    dev.off()
  
  }
  
  df_complete
}

# split data on usual and unusual observations
filterOutliers <- function(df, outliers, mets) {
  subsetOutlier <- with(outliers, eval(parse(text=paste(mets, collapse="|"))))
  list(regulars=df[which(!subsetOutlier),], outliers=df[which(subsetOutlier),])
}





