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
