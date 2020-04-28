library(R0)

setwd("..")
states <- c("KL", "MH", "DL", "RJ", "TN", "UP", "TG", "AP", "KA", "JK", "HR", "PB")#, "GJ", "WB", "MP","BR")
df = read.csv("Data/States.csv")
df[is.na(df)] <- 0
sink("Output/R0/State/state_output.txt")
for (val in states){
  cat(val[1])
  dir.create(paste("E:/Coronavirus/Fitting Model/final model/Rt Analysis/R0/",val,sep = ''), showWarnings = FALSE)
  epid.count = as.numeric(unlist(c(df[val])))
  names(epid.count) = as.Date(as.character(unlist(c(df['ï..date']))))
  #epid.count = epid.count[epid.count<1000]
  GT.flu = generation.time("gamma", c(3.96,4.75))
  res.R = estimate.R(epid.count, GT=GT.flu, methods=c("ML"))
  
  
  fileConn<-file(paste("Output/R0/State/",val,"/output.txt",sep = ''))
  writeLines(c(as.character(res.R$estimates$ML$R),
               as.character(res.R$estimates$ML$conf.int)), fileConn)
  close(fileConn)
  
  cat(",")
  cat(res.R$estimates$ML$R)
  cat(",")
  cat(res.R$estimates$ML$conf.int[1])
  cat(",")
  cat(res.R$estimates$ML$conf.int[2])
  cat(",")
  cat(res.R$estimates$ML$Rsquared)
  cat("\n")
}
sink()
warnings()

