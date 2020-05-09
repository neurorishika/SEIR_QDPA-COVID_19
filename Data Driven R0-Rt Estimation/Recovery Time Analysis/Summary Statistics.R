setwd(".")
reporting = read.csv("reporting_delay.csv")
fitvalues = as.numeric(unlist(c(reporting['lag'])))

bsmean <- function(x,B){
  bstrap <- c()
  for (i in 1:B){
    bstrap <- c(bstrap,mean(sample(x,length(x),replace=T)))
  }
  return(bstrap)
}
output<-bsmean(x = fitvalues,B=1000)
mean(output)
quantile(output,.025)
quantile(output,.975)

bssd <- function(x,B){
  bstrap <- c()
  for (i in 1:B){
    bstrap <- c(bstrap,sd(sample(x,length(x),replace=T)))
  }
  return(bstrap)
}
output<-bssd(x = fitvalues,B=1000)
mean(output)
quantile(output,.025)
quantile(output,.975)

bsmedian <- function(x,B){
  bstrap <- c()
  for (i in 1:B){
    bstrap <- c(bstrap,median(sample(x,length(x),replace=T)))
  }
  return(bstrap)
}
output<-bsmedian(x = fitvalues,B=10000)
mean(output)
quantile(output,.025)
quantile(output,.975)

bsiqr <- function(x,B){
  bstrap <- c()
  for (i in 1:B){
    bstrap <- c(bstrap,IQR(sample(x,length(x),replace=T)))
  }
  return(bstrap)
}
output<-bsiqr(x = fitvalues,B=10000)
mean(output)
quantile(output,.025)
quantile(output,.975)