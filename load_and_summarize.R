mean.bias.g1 = dget("mean.bias.g1.txt")
aggregate(mean.bias.g1[,,])
mean.bias.g2 = dget("mean.bias.g2.txt")
aggregate(mean.bias.g2[,,])

linear.estimatest.g1= dget("linear.estimates.g1.txt")
apply(linear.estimates.g1[,,],c(2,3),mean, na.rm=T)

linear.estimates.g2 = dget("linear.estimates.g2.txt")
apply(linear.estimates.g2[,,],c(2,3),mean, na.rm=T)


aggregate(dget("pointwise.bias.g1.txt"))
pointwise.bias.g1 = dget("pointwise.bias.g1.txt")
apply(pointwise.bias.g1, c(2,3), mean, na.rm=T)
pointwise.bias.g2 = dget("pointwise.bias.g2.txt")
apply(pointwise.bias.g2, c(2,3), mean, na.rm=T)


