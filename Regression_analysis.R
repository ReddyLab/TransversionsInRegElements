library(reshape2)
library(car)
data <- read.table("jaspar_Ts_Tv.tab", header=T)
long_data <- melt(data, id=c("name","pos"))
long_data$Tv <- 1
long_data[which(long_data$variable == "AG"),]$Tv = 0
long_data[which(long_data$variable == "CT"),]$Tv = 0
summary(data)
full_model = lm(value ~ variable + name + pos, data=long_data)
summary(full_model)

pdf("Ts_Tv_in_jaspar.pdf")
crPlot(full_model, "variable")
dev.off()

ts_tv_model = lm(value ~ Tv + pos + name, data=long_data)
summary(ts_tv_model)
crPlots(ts_tv_model)
