calculateFraction <- function(n1, n2, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:2), size=n1+n2, prob=c(n1, n2)/(n1+n2), replace=TRUE))
    repCounts[2]/(repCounts[1]+repCounts[2])
  })
  sprintf("%.2f[%.2f-%.2f]", quantile(fracs, prob=0.5), 
          quantile(fracs, prob=0.025),quantile(fracs, prob=0.975))
}

countClasses <- function(data, race) {
  classCounts <- data %>% 
    group_by(ObesityClass, Reclassified) %>% summarise(n=n())
  classCounts <- pivot_wider(classCounts, names_from="Reclassified", 
                             values_from="n")
  classCounts <- classCounts %>% 
    mutate(Race=race) %>%
    mutate(ObesityClass=as.character(ObesityClass)) %>%
    mutate(`pred.: Normal`=Normal) %>%
    mutate(`pred.: Obese`=Obese) %>%
    mutate("fraction pred. obese"=calculateFraction(Normal, Obese))
  classCounts$Normal <- NULL
  classCounts$Obese <- NULL
  classCounts
}
