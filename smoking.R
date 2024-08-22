
#' Takes a wilcox.test object with option "conf.int=TRUE" and returns the pseudo-median 
#' difference between two population in character format.
#' 
#' @param wilcoxTest a wilcox.test object with option "conf.int=TRUE"
#' @param negative Logical adressing whether the roles of the two groups must be reversed 
#' @return An estimate of the pseudo median difference + 95pc CI
#' @examples
#' 
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' 
#' dataA <- data.frame(group=rep("A", times=20),
#'                     val=rnorm(n=20, mean=0, sd=1))
#' dataB <- data.frame(group=rep("B", times=20),
#'                     val=rnorm(n=20, mean=1, sd=1))
#' data <- bind_rows(dataA, dataB)
#' wilcoxTest <- wilcox.test(formula=val~group, data=data, conf.int=TRUE)
#' confIntMedian(wilcoxTest, negative=FALSE)
#' confIntMedian(wilcoxTest, negative=TRUE)
#' 
confIntMedian <- function(wilcoxTest, negative=FALSE) {
  if (negative) {
    return(sprintf("%.2f[%.2f %.2f]", 
                   -wilcoxTest$estimate, 
                   -wilcoxTest$conf.int[2], 
                   -wilcoxTest$conf.int[1]))
  } else {
    return(sprintf("%.2f[%.2f %.2f]", 
                   wilcoxTest$estimate, 
                   wilcoxTest$conf.int[1], 
                   wilcoxTest$conf.int[2]))
  }
  
}


#' Visualize the main metabolite contributions in the discrepancies in prediction 
#' between smokers and non-smokers.
#' 
#' @param data_smoking A data frame with metabolite measurements and predicted BMI of smoking patients
#' @param data_test A data frame with metabolite measurements and predicted BMI of non-smoking test set patients
#' @return A ggplot object visualizing the main metabolic contributions of smoking to the discrepancy in predicted BMI
#' @examples
#' 
plotSmokingContributions <- function(data_smoking, data_test, model) {
  metDiff <- vector(mode="numeric", length=length(metabolites)+1)
  names(metDiff) <- c("Age", metabolites)
  metDiffSd <- vector(mode="numeric", length=length(metabolites)+1)
  names(metDiffSd) <- c("Age", metabolites)
  for (met in c("Age", metabolites)) {
    nonSmokers <- subset(data_test, subset = ObesityClass=="Normal weight")[[met]]
    smokers <- subset(data_smoking, subset = ObesityClass=="Normal weight")[[met]]
    mean1 <- mean(nonSmokers)
    mean2 <- mean(smokers)
    met_beta <- model$beta[met, "s0"]
    metDiff[met] <- (mean(smokers)-mean(nonSmokers))*met_beta
    metDiffSd[met] <- sqrt(var(smokers)/length(smokers)+var(nonSmokers)/length(nonSmokers))*met_beta
    
  }
  
  fracs <- data.frame(met=names(metDiff), 
                      diffs=metDiff, 
                      lCI=metDiff-1.96*metDiffSd, 
                      uCI=metDiff+1.96*metDiffSd)
  w <- order(abs(fracs$diffs), decreasing=TRUE)   # from large effects to small effects
  fracs <- fracs[w,]
  fracs <- fracs[1:20,]
  
  ggplot(data=fracs, 
         aes(x=reorder(met,abs(diffs)), y=diffs, ymin=lCI, ymax=uCI)) + 
    geom_pointrange(stat="identity") +
    geom_hline(yintercept=0, lty="dashed") +
    #theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.title = element_blank()) +
    xlab("") + ylab("contribution") +
    coord_flip()
}

