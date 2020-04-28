library(EpiEstim)
setwd("..")
df = read.csv("Data/India.csv")
epid.count = setNames(data.frame(matrix(ncol = 2, nrow = 48)),c("local", "imported"))
epid.count$local = as.numeric(unlist(c(df['new_cases'])))
epid.count$imported = as.numeric(unlist(c(df['imported'])))
window = 5
nrows = 48
res.R = estimate_R(epid.count,method = "uncertain_si",config = make_config(list(t_start=2:(nrows-window),t_end=(2+window):nrows,mean_si = 3.96, std_mean_si = 0.215,
                                                                                min_mean_si = 3.53, max_mean_si = 4.39,
                                                                                std_si = 4.75, std_std_si = 0.145,
                                                                                min_std_si = 4.46, max_std_si = 5.07,
                                                                                n1 = 468, n2 = 468,mean_prior=2.6,
                                                                                std_prior=2)))
write.csv(res.R$R,'Output/Rt/India/rtoutput.csv',row.names = TRUE)
plot(res.R)

