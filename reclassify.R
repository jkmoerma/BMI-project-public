calculateFraction <- function(n1, n2, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:2), size=n1+n2, prob=c(n1, n2)/(n1+n2), replace=TRUE))
    repCounts[2]/(repCounts[1]+repCounts[2])
  })
  sprintf("%.2f[%.2f-%.2f]", quantile(fracs, prob=0.5), 
          quantile(fracs, prob=0.025),quantile(fracs, prob=0.975))
}
