calculateFraction <- function(n1, n2, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:2), size=n1+n2, prob=c(n1, n2)/(n1+n2), replace=TRUE))
    repCounts[2]/(repCounts[1]+repCounts[2])
  })
  sprintf("%.2f[%.2f-%.2f]", quantile(fracs, prob=0.5), 
          quantile(fracs, prob=0.025),quantile(fracs, prob=0.975))
}

countClasses <- function(data, race, strat) {
  classCounts <- data %>% 
    group_by(ObesityClass, Reclassified) %>% summarise(n=n())
  classCounts <- pivot_wider(classCounts, names_from="Reclassified", 
                             values_from="n")
  classCounts <- classCounts %>% 
    mutate(Race=race) %>%
    mutate(model=strat) %>%
    mutate(ObesityClass=as.character(ObesityClass)) %>%
    mutate(`pred.: Normal`=Normal) %>%
    mutate(`pred.: Obese`=Obese) %>%
    mutate("fraction pred. obese"=calculateFraction(Normal, Obese))
  classCounts$Normal <- NULL
  classCounts$Obese <- NULL
  classCounts
}

calculateFractionNumeric <- function(n1, n2, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:2), size=n1+n2, prob=c(n1, n2)/(n1+n2), replace=TRUE))
    repCounts[2]/(repCounts[1]+repCounts[2])
  })
  quantile(fracs, prob=c(0.025, 0.5, 0.975))
}

countClassesNumeric <- function(data, race, strat) {
  classCounts <- data %>% 
    group_by(ObesityClass, Reclassified) %>% summarise(n=n())
  classCounts <- pivot_wider(classCounts, names_from="Reclassified", 
                             values_from="n")
  classCounts <- classCounts %>% 
    mutate(Race=race) %>%
    mutate(model=strat) %>%
    mutate(ObesityClass=as.character(ObesityClass)) %>%
    mutate(`pred.: Normal`=Normal) %>%
    mutate(`pred.: Obese`=Obese) %>%
    mutate("fraction pred. obese (lCI)"=calculateFractionNumeric(Normal, Obese)[1],
           "fraction pred. obese (point)"=calculateFractionNumeric(Normal, Obese)[2],
           "fraction pred. obese (uCI)"=calculateFractionNumeric(Normal, Obese)[3])
  classCounts$Normal <- NULL
  classCounts$Obese <- NULL
  classCounts
}



if (FALSE) {
  
  # is AUC of test set predictions inferior to AUC of training set?
  
  aucTrain <- vector(mode="numeric", length=100)
  aucTest <- vector(mode="numeric", length=100)
  fracObeseTrain <- vector(mode="numeric", length=100)
  fracObeseTest <- vector(mode="numeric", length=100)
  meanAgeTrain <- vector(mode="numeric", length=100)
  meanAgeTest <- vector(mode="numeric", length=100)
  
  for (i in 1:100) {
    set.seed(i)
    valSelectA <- sample(x=1:nrow(data_white), size=round(0.2*nrow(data_white)))
    data_white_test <- data_white[valSelectA,]
    data_white_train <- data_white[-valSelectA,]
    modelWhiteTrain <- 
      trainRidgeLASSO(effect=formulationWhite$effect, type=formulationWhite$type, 
                      transformation=formulationWhite$transformation, 
                      new.data=oversample(data_white_train))
    data_white_train$predicted <- 
      exp(predict(object=modelWhiteTrain, 
                  newx=makeMatrix(data_white_train, 
                                  includeInteraction=formulationWhite$effect=="interaction")$mat,
                  type="response")[,"s0"])
    data_white_val$predicted <- 
      exp(predict(object=modelWhiteTrain, 
                  newx=makeMatrix(data_white_val, includeInteraction=formulationWhite$effect=="interaction")$mat,
                  type="response")[,"s0"])
    rocTrainWhite <- roc(data_white_train$ObesityClass, data_white_train$predicted, 
                         levels=c("Normal weight", "Obese"))
    cutoffWhite <- coords(rocTrainWhite, x="best")[1,"threshold"]
    data_white_val <- data_white_val %>% 
      mutate(selectROC=c("omitted", "control", "case")[1+(ObesityClass=="Normal weight")+2*(ObesityClass=="Obese")]) %>%
      mutate(Reclassified=c("Normal", "Obese")[1+(predicted>cutoffWhite)])
    rocWhiteBMIclass <- roc(formula=ObesityClass~predicted, data=data_white_val, 
                            levels=c("Normal weight", "Obese"))
    aucTest[i] <- auc(rocWhiteBMIclass)
    statsWhite <- ci.coords(rocWhiteBMIclass, x=cutoffWhite, input="threshold", 
                            ret=c("sensitivity", "specificity"))
    aucTrain[i] <- auc(rocTrainWhite)
    
    classCountsTrain <- table(data_white_train$ObesityClass)
    fracObeseTrain[i] <- classCountsTrain["Obese"]/(classCountsTrain["Normal weight"]+classCountsTrain["Obese"])
    
    classCountsTest <- table(data_white_test$ObesityClass)
    fracObeseTest[i] <- classCountsTest["Obese"]/(classCountsTest["Normal weight"]+classCountsTest["Obese"])
    
  }
  ggplot(data=data.frame("AUC train"=aucTrain, 
                         "AUC validation"=aucTest,
                         "fraction"=fracObeseVal,
                         check.names=FALSE),
         mapping=aes(x=`AUC train`, y=`AUC validation`, col=fraction)) +
    scale_fill_continuous(type="gradient", palette="PuBuGn") +
    geom_point() +
    geom_abline(slope=1, intercept=0, lty="dashed")
  quantile(aucTrain, probs=c(0.025, 0.5, 0.975))
  quantile(aucTest, probs=c(0.025, 0.5, 0.975))
}


