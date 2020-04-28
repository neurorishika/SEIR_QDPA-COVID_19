library(R0)

setwd("..")
df = read.csv("Data/R0_India.csv")
epid.count = as.numeric(unlist(c(df['active'])))
import = as.numeric(unlist(c(df['imported'])))
GT.flu = generation.time("gamma", c(3.96,4.75))

sense = sensitivity.analysis(epid.count, import = import, GT.flu, begin=1:(length(epid.count)/2), end=((length(epid.count)/2)+1):length(epid.count), est.method="ML", sa.type="time")
png(file="Output/R0/India/MLTS.png",width = 4, height = 4, res = 300, units = 'in',pointsize = 6)
plot(sense)
dev.off()

best1<-plot(sense, what=c("criterion"))

sense2 = sensitivity.analysis(epid.count, import = import, begin =best1$best.fit$Begin.dates, end = best1$best.fit$End.dates, GT.type="gamma", GT.mean=seq(2,7,1), GT.sd.range=3, est.method="ML", sa.type="GT")
print(sense2)
png(file="Output/R0/India/ML_Gamma.png",width = 4, height = 4, res = 300, units = 'in',pointsize = 6)
plot(as.numeric(sense2[,2]), as.numeric(sense2[,4]),
     ylim=range(c(as.numeric(sense2[,5]), 
                  as.numeric(sense2[,6]))),
     xlim=range(c(as.numeric(sense2[,2])-as.numeric(sense2[,3]), 
                  as.numeric(sense2[,2])+as.numeric(sense2[,3]))),
     pch=19, xlab="Gamma (Mean +/- SD)", ylab="R0 (with 95% CI)",
     main="Sensitivity Analysis for Gamma Distribution")
arrows(as.numeric(sense2[,2]), as.numeric(sense2[,5]), as.numeric(sense2[,2]), as.numeric(sense2[,6]), length=0.05, angle=90, code=3)
arrows(as.numeric(sense2[,2])-as.numeric(sense2[,3]), as.numeric(sense2[,4]), as.numeric(sense2[,2])+as.numeric(sense2[,3]), as.numeric(sense2[,4]), length=0.05, angle=90, code=3)
dev.off()

sense = sensitivity.analysis(epid.count, import = import, GT.flu, begin=1:(length(epid.count)/2), end=((length(epid.count)/2)+1):length(epid.count), est.method="EG", sa.type="time")
png(file="Output/R0/India/EGTS.png",width = 4, height = 4, res = 300, units = 'in',pointsize = 6)
plot(sense)
dev.off()


sense2 = sensitivity.analysis(epid.count, import = import, begin =best2$best.fit$Begin.dates, end = best2$best.fit$End.dates, GT.type="gamma", GT.mean=seq(2,7,1), GT.sd.range=3, est.method="EG", sa.type="GT")
print(sense2)
png(file="Output/R0/India/EG_Gamma.png",width = 4, height = 4, res = 300, units = 'in',pointsize = 6)
plot(as.numeric(sense2[,2]), as.numeric(sense2[,4]),
     ylim=range(c(as.numeric(sense2[,5]), 
                  as.numeric(sense2[,6]))),
     xlim=range(c(as.numeric(sense2[,2])-as.numeric(sense2[,3]), 
                  as.numeric(sense2[,2])+as.numeric(sense2[,3]))),
     pch=19, xlab="Gamma (Mean +/- SD)", ylab="R0 (with 95% CI)",
     main="Sensitivity Analysis for Gamma Distribution")
arrows(as.numeric(sense2[,2]), as.numeric(sense2[,5]), as.numeric(sense2[,2]), as.numeric(sense2[,6]), length=0.05, angle=90, code=3)
arrows(as.numeric(sense2[,2])-as.numeric(sense2[,3]), as.numeric(sense2[,4]), as.numeric(sense2[,2])+as.numeric(sense2[,3]), as.numeric(sense2[,4]), length=0.05, angle=90, code=3)
dev.off()


fileConn<-file("Output/R0/India/r0output.txt")
writeLines(c(as.character(best1$best.fit$R),
             as.character(best1$best.fit$CI.lower),
             as.character(best1$best.fit$CI.upper),
             as.character(best1$best.fit$Rsquared),
             as.character(best2$best.fit$R),
             as.character(best2$best.fit$CI.lower),
             as.character(best2$best.fit$CI.upper),
             as.character(best2$best.fit$Rsquared)), fileConn)
close(fileConn)

  #plot(res.R)
plotfit(res.R)
