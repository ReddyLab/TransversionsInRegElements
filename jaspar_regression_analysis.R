library(reshape2)
library(ggplot2)
library(plyr)
library(car)

sink("jaspar_regression_analysis.log")
data <- read.table("jaspar_Ts_Tv.tab", header=T)
long_data <- melt(data, id=c("name","pos","IC", "DegCons"))
long_data$mut_type <- factor("Tv", levels=c("Ts","Tv"))
long_data[which(long_data$variable == "AG"),]$mut_type = "Ts"
long_data[which(long_data$variable == "CT"),]$mut_type = "Ts"
summary(long_data)

#
# Diagnostic plots
#

# Diagnostic plot 1: effect size relative to motif center.
# Observation: effects are strong in motif center, as expected.
pdf("positional_effect_of_mutations_on_pssm_score.pdf")
boxplot(value ~ as.integer(5*pos), data=long_data, xlab = "Normalized distance from motif center (quintiles)", ylab="|difference in log-odds score|")
dev.off();

# Diagnostic plot 2: effect size of transversion
# Observation: transversions appear to have a stronger effect
pdf("boxplot_of_transversion_effects_on_pssm_score.pdf")
boxplot(value ~ mut_type, data=long_data)
by(long_data$value, long_data$mut_type, summary)
dev.off()

#
# Main figure #1: The effect is particularly pronounced near the center of the motif
#

# binning:
bin_size = 0.2
long_data$bin=as.factor(round_any((long_data$pos), bin_size))

# regression within each bin
regression_stats <- ddply(long_data, .(bin), summarize, 
      estimate  = summary(lm(value ~ mut_type))$coefficients[2,1], 
      std_error = summary(lm(value ~ mut_type))$coefficients[2,2],
      tstat     = summary(lm(value ~ mut_type))$coefficients[2,3], 
      pval      = summary(lm(value ~ mut_type))$coefficients[2,4])

# plot theme
plot_theme <- theme_bw() +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
theme(axis.line = element_line(colour = "black"))

# box plots
pdf("Figure_1_A_Transversions_alter_TF_binding_more.pdf", width=10, height=3)
ggplot(long_data, aes(x=bin, y=value, fill=mut_type)) + 
geom_boxplot(notch=T) + 
plot_theme + 
scale_fill_manual(values = c("grey","white")) + 
labs(x="Normalized distance from motif center", y="|difference in log-odds score|")+
annotate("text", x = regression_stats[,1], y = 22, label=sprintf("%1.1g", regression_stats$pval))
dev.off()

#
# Main figure #2: The effect is also pronounced in most informative positions
#
IC_bin_size = 0.1
long_data$IC_bin=as.factor(round_any(long_data$IC/min(long_data$IC), IC_bin_size))

# regression within each bin
regression_stats <- ddply(long_data, .(IC_bin), summarize, 
      estimate  = summary(lm(value ~ mut_type))$coefficients[2,1], 
      std_error = summary(lm(value ~ mut_type))$coefficients[2,2],
      tstat     = summary(lm(value ~ mut_type))$coefficients[2,3], 
      pval      = summary(lm(value ~ mut_type))$coefficients[2,4])

pdf("Figure_1_B_Transversions_alter_degenerate_sites_more.pdf", width=10, height=3)
ggplot(long_data, aes(x=IC_bin, y=value, fill=mut_type)) + 
geom_boxplot(notch=T) + 
plot_theme + 
scale_fill_manual(values = c("grey","white")) + 
labs(x="Normalized information content", y="|difference in log-odds score|") +
annotate("text", x = regression_stats[,1], y = 22, label=sprintf("%1.1g", regression_stats$pval))
dev.off()

save.image("Jaspar_ts_tv.rdata")


