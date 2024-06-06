calculateFraction <- function(n1, n2, n3, reps=1000) {
  fracs <- replicate(n=reps, {
    repCounts <- table(sample(factor(1:3), size=n1+n2+n3, prob=c(n1, n2, n3)/(n1+n2+n3), replace=TRUE))
    repCounts[3]/(repCounts[1]+repCounts[3])
  })
  sprintf("%.2f[%.2f-%.2f]", quantile(fracs, prob=0.5), 
          quantile(fracs, prob=0.025),quantile(fracs, prob=0.975))
}
