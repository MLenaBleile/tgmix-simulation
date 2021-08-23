

pointwise.bias.g1 = dget("pointwise.bias.g1.txt")
plot_pointwise_bias = function(pointwise.bias.g1, main="Group 1 Pointwise Bias", xlim=c(0,40)){
summarized.pointwise.bias.g1 = apply(pointwise.bias.g1[1:500,,], c(2,3), mean, na.rm=T)
mcerr.pointwise.bias.g1 = abs(summarized.pointwise.bias.g1 - apply(pointwise.bias.g1[501:1000,,], c(2,3), mean, na.rm=T))
days = seq(0,40, by=4)
plot(days, summarized.pointwise.bias.g1['mi',], col="purple", type="l", main=main, xlab="Days", ylab="Bias",xlim=xlim, ylim=c(-5,0))
lines(days, mcerr.pointwise.bias.g1['mi',]+summarized.pointwise.bias.g1['mi',], col="purple", lty=2)
lines(days, -mcerr.pointwise.bias.g1['mi',]+summarized.pointwise.bias.g1['mi',], col="purple", lty=2)

lines(days, summarized.pointwise.bias.g1['lvcf',], col="red")
lines(days, mcerr.pointwise.bias.g1['lvcf',]+summarized.pointwise.bias.g1['lvcf',], col="red", lty=2)
lines(days, -mcerr.pointwise.bias.g1['lvcf',]+summarized.pointwise.bias.g1['lvcf',], col="red", lty=2)
legend("bottomleft", legend=c("LVCF","MI","MC Error"), col=c("red","purple","black"), lty=c(1,1,2))

}

plot_pointwise_bias(pointwise.bias.g1)

pointwise.bias.g2 = dget("pointwise.bias.g2.txt")
plot_pointwise_bias(pointwise.bias.g2, main="Group 2 Pointwise Bias", xlim=c(0,25))
