M=1000
n=5
np1 = n+1
twon = 2*n
time=seq(0,40,by=4)
theta = c(5,0.25,1,1)
theta_a = c(5,0.35, 1,1)
simulation_iteration = as.integer(read.table("simulation_iteration.txt"))
set.seed(1998+simulation_iteration)
cat(1998 + simulation_iteration," ",simulation_iteration,"\n", file="random_seeds.txt", append=TRUE)

mouse.sim = array(data=NA, dim=c(M,n*2, length(time)+1))

for(m in 1:M){
  for(i in 1:n){
    mouse.sim[m,i,] = c("group 1", generate.one(time, theta))
  }
  for(i in np1:twon){
    mouse.sim[m, i,] = c("group 2",generate.one(time, theta_a))
  }
}

if (simulation_iteration == 1){
TotalTime=length(time)
pointwise.bias.g1 = array(dim=c(M, 2, TotalTime), dimnames=list(NULL, c("mi","lvcf"), paste("Day",time)))
pointwise.bias.g2 = array(dim=c(M, 2, TotalTime), dimnames=list(NULL, c("mi","lvcf"), paste("Day",time)))
mean.bias.g1 =array(dim=c(M, 2, TotalTime), dimnames=list(NULL, c("mi","cc"), paste("Day",time)))
mean.bias.g2 = array(dim=c(M, 2, TotalTime), dimnames=list(NULL, c("mi","cc"), paste("Day",time)))
linear.estimates.g1 = array(dim=c(M, 2, 3), dimnames=list(NULL, c("mi","cv"), c("intercept","estimate", "mse")))
linear.estimates.g2 = array(dim=c(M, 2, 3), dimnames=list(NULL, c("mi","cv"), c("intercept","estimate", "mse")))

main.analysis.gp.diff = data.frame(lvcf=NA, cc=NA, cv=NA, mi=NA)}else if(simulation_iteration >1){
    pointwise.bias.g1 = dget("pointwise.bias.g1.txt")
    pointwise.bias.g2 = dget("pointwise.bias.g2.txt")
    mean.bias.g1 = dget("mean.bias.g1.txt")
    mean.bias.g2 = dget("mean.bias.g2.txt")
    linear.estimates.g1 = dget("linear.estimates.g1.txt")
    linear.estimates.g2 = dget("linear.estimates.g2.txt")
    main.analysis.gp.diff = dget("main.analysis.gp.diff.txt")
  
}