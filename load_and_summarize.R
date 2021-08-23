pdf(file="bias_plots.pdf")
par(mfrow=c(2,2))
mean.bias.g1 = dget("mean.bias.g1.txt")
make_mean_bias_plot = function(mean.bias.g1, main="Group 1 Mean Bias", alt='cc', daylength=11){
mean.bias.g1 = mean.bias.g1[,,1:daylength]
summarized.mean.bias.g1 = apply(mean.bias.g1[1:500,,], c(2,3), mean, na.rm=T)
days = seq(0,40, length.out = 11)[1:length(summarized.mean.bias.g1['mi',])]
mcerror = summarized.mean.bias.g1 - apply(mean.bias.g1[501:1000,,], c(2,3), mean, na.rm=T)
print(mcerror)
plot(days, rm.na(summarized.mean.bias.g1['mi',]), type="l", xlab="Days", ylab="Bias", col="purple", main=main, ylim=range(min(summarized.mean.bias.g1[alt,], na.rm=T),0))
lines(days, summarized.mean.bias.g1[alt,1:length(days)], col="red")
legend("bottomleft", col=c("red", "purple","black"), lty=c(1,1,2),legend=c(alt,"MI","MC Error"))
lines(days, summarized.mean.bias.g1['mi',1:length(days)] + mcerror['mi',1:length(days)], lty=2, col="purple")
lines(days, summarized.mean.bias.g1['mi',] - mcerror['mi',], lty=2, col="purple")
lines(days, summarized.mean.bias.g1[alt,] + mcerror[alt,], lty=2, col="red")
lines(days, summarized.mean.bias.g1[alt,] - mcerror[alt,], lty=2, col="red")
}
make_mean_bias_plot(mean.bias.g1)

mean.bias.g2 = dget("mean.bias.g2.txt")
apply(mean.bias.g2[1:500,,1:8], c(2,3), mean, na.rm=T)
make_mean_bias_plot(mean.bias.g2[,,1:8], main="Group 2 Mean Bias", daylength = 8)

linear.estimates.g1= dget("linear.estimates.g1.txt")
summarized.linear.estimates.g1 = apply(linear.estimates.g1[1:500,,],c(2,3),mean, na.rm=T)
summarized.linear.estimates.g1
lest.mcerror.g1 = apply(linear.estimates.g1[1:500,,],c(2,3),mean, na.rm=T) - apply(linear.estimates.g1[501:1000,,],c(2,3),mean, na.rm=T)
abs(lest.mcerror.g1)


linear.estimates.g2 = dget("linear.estimates.g2.txt")
apply(linear.estimates.g2[1:500,,],c(2,3),mean, na.rm=T)
summarized.linear.estimates.g2 = apply(linear.estimates.g2[1:500,,],c(2,3),mean, na.rm=T)
summarized.linear.estimates.g2
lest.mcerror.g2 = apply(linear.estimates.g2[1:500,,],c(2,3),mean, na.rm=T) - apply(linear.estimates.g2[501:1000,,],c(2,3),mean, na.rm=T)
abs(lest.mcerror.g2)



source("plot_pointwise_bias.R")
dev.off()

