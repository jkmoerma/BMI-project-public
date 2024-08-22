#' Illustrates the interpretation of the figures on metabolic differences of the NO outliers
#' 
#' @return A ggplot object
illustrNO <- function() {
  ggplot(data=data.frame("OO-NO"=c(0,1, -1, 2), "NO-NN"=c(1,0, 2, -1), check.names=FALSE),
         mapping=aes(x=`OO-NO`, y=`NO-NN`)) + 
    geom_point() +
    geom_abline(intercept=0, slope=-1, lty="dashed") + 
    geom_abline(intercept=1, slope=-1, lty="dashed") + 
    geom_vline(xintercept=0) + 
    geom_hline(yintercept=0) +
    geom_segment(data=data.frame("OO-NO"=0, "NO-NN"=0, 
                                 "OO-NO end"=1, "NO-NN end"=0, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="both", type="open"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="purple")) + 
    geom_text(data=data.frame("OO-NO"=0.5, "NO-NN"=0, 
                              "OO-NO end"=1, "NO-NN end"=0, 
                              check.names=FALSE),
              aes(col="purple"),
              label="OO-NN",
              nudge_y= 0.05) +
    geom_segment(data=data.frame("OO-NO"=c(0,1), "NO-NN"=c(1,0), 
                                 "OO-NO end"=c(0, 1)+0.4, "NO-NN end"=c(1,0)+0.4, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="blue")) + 
    geom_text(data=data.frame("OO-NO"=c(0, 1)+0.4, "NO-NN"=c(1, 0)+0.4, 
                              check.names=FALSE),
              aes(col="blue"),
              label=c("ideal model", "random guessing BMI"),
              nudge_x= 0.15, nudge_y= 0.15) +
    geom_segment(data=data.frame("OO-NO"=1, "NO-NN"=0, 
                                 "OO-NO end"=0, "NO-NN end"=1, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="green")) + 
    geom_text(data=data.frame("OO-NO"=0.5, "NO-NN"=0.5, 
                              check.names=FALSE),
              aes(col="green"),
              label="suboptimal model",
              nudge_x= 0.1, nudge_y= 0.1, angle=-45) +
    geom_segment(data=data.frame("OO-NO"=c(1,0), "NO-NN"=c(0,1), 
                                 "OO-NO end"=c(2, -1), "NO-NN end"=c(-1,2), 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="red")) + 
    geom_text(data=data.frame("OO-NO"=c(-0.5, 1.5), "NO-NN"=c(1.5, -0.5), 
                              check.names=FALSE),
              aes(col="red"),
              label=c("abberant in NO", "abberant in NO"),
              nudge_x= 0.1, nudge_y= 0.1, angle=-45) +
    scale_color_manual(values=c("blue", "green", "purple", "red")) +
    theme(legend.position = "none") +
    coord_fixed(ratio=1) +
    labs(title="NO outliers: observations and expectations")
}

#' Illustrates the interpretation of the figures on metabolic differences of the ON outliers
#' 
#' @return A ggplot object
illustrON <- function() {
  ggplot(data=data.frame("OO-NO"=c(0,-1, 1, -2), "NO-NN"=c(1,0, 2, -1), check.names=FALSE),
         mapping=aes(x=`OO-NO`, y=`NO-NN`)) + 
    geom_point() +
    geom_abline(intercept=0, slope=1, lty="dashed") + 
    geom_abline(intercept=1, slope=1, lty="dashed") + 
    geom_vline(xintercept=0) + 
    geom_hline(yintercept=0) +
    geom_segment(data=data.frame("OO-NO"=0, "NO-NN"=0, 
                                 "OO-NO end"=0, "NO-NN end"=1, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="both", type="open"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="purple")) + 
    geom_text(data=data.frame("OO-NO"=0, "NO-NN"=0.5, 
                              check.names=FALSE),
              aes(col="purple"),
              label="OO-NN",
              nudge_x= -0.05, angle=90) +
    geom_segment(data=data.frame("OO-NO"=c(0,-1), "NO-NN"=c(1,0), 
                                 "OO-NO end"=c(0, -1)-0.4, "NO-NN end"=c(1,0)+0.4, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0.1, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="blue")) + 
    geom_text(data=data.frame("OO-NO"=c(0, -1)-0.4, "NO-NN"=c(1, 0)+0.4, 
                              check.names=FALSE),
              aes(col="blue"),
              label=c("ideal model", "random guessing BMI"),
              nudge_x= -0.15, nudge_y= 0.15) +
    geom_segment(data=data.frame("OO-NO"=-1, "NO-NN"=0, 
                                 "OO-NO end"=0, "NO-NN end"=1, 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="green")) + 
    geom_text(data=data.frame("OO-NO"=-0.5, "NO-NN"=0.5, 
                              check.names=FALSE),
              aes(col="green"),
              label="suboptimal model",
              nudge_x= -0.1, nudge_y= 0.1, angle=45) +
    geom_segment(data=data.frame("OO-NO"=c(-1,0), "NO-NN"=c(0,1), 
                                 "OO-NO end"=c(-2, 1), "NO-NN end"=c(-1,2), 
                                 check.names=FALSE),
                 arrow=arrow(angle=30, length=unit(0, "inches"),
                             ends="last", type="closed"),
                 aes(x=`OO-NO`, y=`NO-NN`, xend=`OO-NO end`, yend=`NO-NN end`, col="red")) + 
    geom_text(data=data.frame("OO-NO"=c(0.5, -1.5), "NO-NN"=c(1.5, -0.5), 
                              check.names=FALSE),
              aes(col="red"),
              label=c("abberant in ON", "abberant in ON"),
              nudge_x= -0.1, nudge_y= 0.1, angle=45) +
    scale_color_manual(values=c("blue", "green", "purple", "red")) +
    theme(legend.position = "none") +
    coord_fixed(ratio=1) +
    labs(title="ON outliers: observations and expectations", x="NN-ON", y="OO-ON")
}


#' Calculates the difference in metabolic levels between ON, NO, OO and NN.
#' A plot is made of the scaled metabolic differences for comparison.
#' 
#' @param data Data with metabolite levels and assigned prediction group "predictionGroup"
#' @param ethnicity Ethnicity of the patients in the data set
#' @return A list of 2 ggplot objects for NO ("normal") and ON ("obese") outliers
#' @examples
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(glmnet)
#' library(pROC)
#' library(doParallel)
#' library(ggplot2)
#' 
#' metabolites <- c("met1", "met2", "met_coll", "met_110")
#' 
#' data_normal <- data.frame(ObesityClass=rep("Normal weight", times=200),
#'                           Age=rnorm(n=200, mean=30, sd=3),
#'                           BMI=runif(n=200, min=18, max=25),
#'                           met1=rnorm(n=200, mean=0.5, sd=0.5),
#'                           met2=rnorm(n=200, mean=5, sd=1),
#'                           met_110=rnorm(n=200, mean=100, sd=100))
#' data_overweight <- data.frame(ObesityClass=rep("Overweight", times=100),
#'                           Age=rnorm(n=100, mean=30, sd=3),
#'                           BMI=runif(n=100, min=25, max=30),
#'                           met1=rnorm(n=10, mean=0.7, sd=0.5),
#'                           met2=rnorm(n=100, mean=4, sd=1),
#'                           met_110=rnorm(n=100, mean=100, sd=100))
#' data_obese <- data.frame(ObesityClass=rep("Obese", times=50),
#'                           Age=rnorm(n=50, mean=30, sd=3),
#'                           BMI=runif(n=50, min=30, max=40),
#'                           met1=rnorm(n=50, mean=1, sd=0.5),
#'                           met2=rnorm(n=50, mean=3, sd=1),
#'                           met_110=rnorm(n=50, mean=100, sd=100))
#'                         
#' data_train <- bind_rows(data_normal, data_overweight, data_obese)
#' data_train$BMI <- log(data_train$BMI)
#' data_train$met_coll <- data_train$met1 + 7.2*data_train$met2 + rnorm(n=nrow(data_train), mean=0, sd=0.1)
#' 
#' model <- trainRidgeLASSO(effect="main", type="Ridge", transformation="Log", new.data=data_train)
#' data_train$predicted <- 
#'   exp(predict(model, 
#'               newx=makeMatrix(data_train, 
#'                               includeInteraction=FALSE)$mat)[,"s0"])
#' 
#' rocTrain <- roc(data_train$ObesityClass, data_train$predicted, 
#'                 levels=c("Normal weight", "Obese"))
#' cutoff <- coords(rocTrain, x="best")[1,"threshold"]
#' data_train <- data_train %>% 
#'   mutate(predictionGroup = ifelse(ObesityClass=="Normal weight"&predicted<cutoff, "NN", 
#'                              ifelse(ObesityClass=="Normal weight"&predicted>cutoff, "NO", 
#'                                ifelse(ObesityClass=="Obese"&predicted>cutoff, "OO", 
#'                                  ifelse(ObesityClass=="Obese"&predicted<cutoff, "ON", NA)))))
#' 
#' ggplot(data_train, aes(x=predicted, y=exp(BMI), col=predictionGroup)) + 
#'   geom_point() + 
#'   labs(x="predicted BMI", y="observed BMI")
#' 
#' plotANOVA(data_train, "White")
#' 
plotANOVA <- function(data, ethnicity) {
  for (met in c("Age", metabolites)) {
    data[[met]] <- scale(data[[met]])
  }
  groupLevels <- subset(data, subset=!is.na(predictionGroup), 
                        select=c("Age", metabolites, "predictionGroup"))
  
  levelDiffsNormal <- matrix(nrow=length(metabolites)+1, ncol=6)
  rownames(levelDiffsNormal) <- c("Age", metabolites)
  colnames(levelDiffsNormal) <- c("NO-NN", "NO-NN (lCI)", "NO-NN (uCI)",
                                  "OO-NO", "OO-NO (lCI)", "OO-NO (uCI)")
  
  levelDiffsObese <- matrix(nrow=length(metabolites)+1, ncol=6)
  rownames(levelDiffsObese) <- c("Age", metabolites)
  colnames(levelDiffsObese) <- c("OO-ON", "OO-ON (lCI)", "OO-ON (uCI)",
                                 "NN-ON", "NN-ON (lCI)", "NN-ON (uCI)")
  
  legend <- c()
  for (met in c("Age", metabolites)) {
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
    geom_point() + 
    #geom_pointrange(aes(ymin=`NO-NN (lCI)`, ymax=`NO-NN (uCI)`), fatten=1, lwd=0.1) +
    #geom_pointrange(aes(xmin=`OO-NO (lCI)`, xmax=`OO-NO (uCI)`), fatten=1, lwd=0.1) +
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
    geom_point() + 
    #geom_pointrange(aes(ymin=`OO-ON (lCI)`, ymax=`OO-ON (uCI)`), fatten=1, lwd=0.1) +
    #geom_pointrange(aes(xmin=`NN-ON (lCI)`, xmax=`NN-ON (uCI)`), fatten=1, lwd=0.1) +
    geom_text(vjust=0, nudge_y=0.05, size=2) +
    labs(title=paste0("ON (", ethnicity, " ethnicity)"), 
         x="Scaled log conc. diff. NN-ON", 
         y="Scaled log conc. diff. OO-ON")
  
  list(normal=pN, obese=pO)
}