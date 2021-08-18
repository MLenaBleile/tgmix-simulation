

source("C:/Users/s198663/Documents/TGM/simulation/preamble.R")

source("C:/Users/s198663/Documents/TGM/simulation/generate_simulation.R")






M=1000
for(m in simulation_iteration:M){
one.dataset=mouse.sim[m,,]
beta.v = c(0.1,1,2,-30)
valid_dataset = F
while(!valid_dataset){
data_with_dropout = t(apply(one.dataset[,-c(1)], 1,dropout.process, beta=beta.v))
valid_dataset = is.matrix(make_cv_dataset(data_with_dropout[1:n,]))*sum(rowSums(!is.na(make_cc_dataset(data_with_dropout[1:n,])))>=1)

valid_dataset = valid_dataset*is.matrix(make_cv_dataset(data_with_dropout[-seq(1,n),]))*sum(rowSums(!is.na(make_cc_dataset(data_with_dropout[-seq(1,n),])))>=1)

valid_dataset = valid_dataset*(sum(is.na(make_cc_dataset(data_with_dropout[1:n,])))>0)*(sum(is.na(make_cc_dataset(data_with_dropout[-seq(1,n),])))>0)
}
cat("\n \n Valid dataset 1: ")
data_with_dropout[1:n,]
cat("\n \n Valid dataset 2: ")
data_with_dropout[-seq(1,n),]
ccdataset.g1 = make_cc_dataset(data_with_dropout[1:n,colSums(is.na(data_with_dropout[1:n,]))<n])
data_with_dropout[-seq(1,n),]
ccdataset.g2 = make_cc_dataset(data_with_dropout[-seq(1,n),colSums(is.na(data_with_dropout[-seq(1,n),]))<n])

lvcfdataset.g1 = lvcf(data_with_dropout[1:n,])
lvcfdataset.g2 = lvcf(data_with_dropout[-seq(1,n),])

cvdataset.g1 = make_cv_dataset(data_with_dropout[1:n,])
cvdataset.g2= make_cv_dataset(data_with_dropout[-seq(1,n),])

data_for_mi=filter_all_na(as.data.frame(data_with_dropout))
data_for_mi$Animal=1:nrow(data_with_dropout)
colnames(data_for_mi)[-c(ncol(data_for_mi))] = paste("Day", time, sep="")[1:(ncol(data_for_mi)-1)]
data_for_mi=data_for_mi[,c(ncol(data_for_mi),1:(ncol(data_for_mi)-1))]

midataset.g1 = tgmix_imputation(data_for_mi[1:n,colSums(is.na(data_for_mi[1:n,]))<n], plots=F)
midataset.g2 = tgmix_imputation(data_for_mi[-seq(1,n),colSums(is.na(data_for_mi[-seq(1,n),]))<n], plots=F)
    
trueMeans.g1 = colMeans(as.data.frame(t(apply(one.dataset[1:n, -c(1)],1,as.numeric))))
trueMeans.g2 = colMeans(as.data.frame(t(apply(one.dataset[-seq(1,n), -c(1)],1,as.numeric))))
ccMeans.g1 = colMeans(ccdataset.g1[rowSums(is.na(ccdataset.g1))<1,])
ccMeans.g2 = colMeans(ccdataset.g2[rowSums(is.na(ccdataset.g2))<1,])
tgmixMeans.g1 = colMeans(midataset.g1[,-c(1:2)])
tgmixMeans.g2 = colMeans(midataset.g2[,-c(1:2)])
tgmixMeanBias.g1 = tgmixMeans.g1 - trueMeans.g1[1:length(tgmixMeans.g1)]
tgmixMeanBias.g2 = tgmixMeans.g2 - trueMeans.g2[1:length(tgmixMeans.g2)]
ccMeanBias.g1 = ccMeans.g1 - trueMeans.g1[1:length(ccMeans.g1)]
ccMeanBias.g2 = ccMeans.g2 - trueMeans.g2[1:length(ccMeans.g2)]
    
mean.bias.g1[m,'cc',1:length(ccMeanBias.g1)] = ccMeanBias.g1
mean.bias.g1[m, 'mi',1:length(tgmixMeanBias.g1)] = tgmixMeanBias.g1
mean.bias.g2[m,'cc',1:length(ccMeanBias.g2)] = ccMeanBias.g2

mean.bias.g2[m, 'mi',1:length(tgmixMeanBias.g2)] = tgmixMeanBias.g2
    
    
    
    
truedataset.g1= as.data.frame(t(apply(one.dataset[1:n, -c(1)],1,as.numeric)))
truedataset.g2 = t(apply(one.dataset[-seq(1,n), -c(1)],1,as.numeric))
lvcfBias.g1 = colMeans(as.data.frame(lvcfdataset.g1[[1]]) - truedataset.g1)
lvcfBias.g2 = colMeans(as.data.frame(lvcfdataset.g2[[1]]) - truedataset.g2)

    
lastcol.g1 =ncol(midataset.g1) - 2
lastcol.g2 =ncol(midataset.g2)-2
miBias.g1 = get_mi_bias(midataset.g1, truedataset.g1[,1:lastcol.g1])
if(sum(miBias.g1[1:2])>0){
    cat("\n \n MI Bias: ")
    miBias.g1
    dput(miBias.g1, "mibiasg1.txt")
    cat("\n\n MI complete dataset:")
    midataset.g1
    dput(midataset.g1, "midatasetg1.txt")
    cat("\n \n True Dataset:")
    truedataset.g1
    dput(truedataset.g1, "truedatasetg1.txt")
    break
}
miBias.g2 = get_mi_bias(midataset.g2, truedataset.g2[,1:lastcol.g2])

pointwise.bias.g1[m,'mi',1:length(miBias.g1)] = miBias.g1
pointwise.bias.g2[m, 'mi',1:length(miBias.g2)] = miBias.g2

pointwise.bias.g1[m,'lvcf', 1:length(lvcfBias.g1)] = lvcfBias.g1
pointwise.bias.g2[m, 'lvcf',1:length(lvcfBias.g2)] = lvcfBias.g2
    
cv.fit.g2 = lme4::fixef(fit_lmer(cvdataset.g2)[[1]])
cv.mse.g2 = fit_lmer(cvdataset.g2)[[2]]

mi.est.and.mse = fit_lmer_multiple(midataset.g2)
mi.fit.est = mi.est.and.mse[[1]]
mi.fit.mse = mi.est.and.mse[[2]]

flmg1 = fit_lmer_multiple(midataset.g1)
linear.estimates.g1[m,'mi',c('intercept','estimate')] = flmg1[[1]]
linear.estimates.g1[m, 'mi',c('mse')] = flmg1[[2]]

fleg1 = fit_lmer(cvdataset.g1)

linear.estimates.g1[m,'cv',c('intercept','estimate')] = lme4::fixef(fleg1[[1]])
linear.estimates.g1[m,'cv',c('mse')] = fleg1[[2]]

flmg2 = fit_lmer_multiple(midataset.g2)
linear.estimates.g2[m,'mi',c('intercept','estimate')] = flmg2[[1]]
linear.estimates.g2[m, 'mi',c('mse')] = flmg2[[2]]

fleg2 = fit_lmer(cvdataset.g2)

linear.estimates.g2[m,'cv',c('intercept','estimate')] = lme4::fixef(fleg2[[1]])
linear.estimates.g2[m,'cv',c('mse')] = fleg2[[2]]


lvcfdataset.g1.new=as.data.frame(lvcfdataset.g1[1])

lvcfdataset.g2.new=as.data.frame(lvcfdataset.g2[1])
lvcf.all = rbind(lvcfdataset.g1.new, lvcfdataset.g2.new)
    

maxcol.cc =  min(c(ncol(ccdataset.g1), ncol(ccdataset.g2)))
cc.all = rbind(as.data.frame(ccdataset.g1[,1:maxcol.cc]), as.data.frame(ccdataset.g2[,1:maxcol.cc]))
maxcol = min(c(ncol(cvdataset.g1), ncol(cvdataset.g2)))
cv.all = rbind(as.data.frame(cvdataset.g1[,1:maxcol]), as.data.frame(cvdataset.g2[,1:maxcol]))
idx=1
for (df in list(lvcf.all, cc.all,cv.all)){
    main.analysis.gp.diff[m,idx] = get_main_analysis(df)[[2]][1]$coefficients[1]
    idx=idx+1
    #print(df)
}

for(imputation in 1:max(midataset.g1$imputation)){
    colcount = min(c(ncol(midataset.g1), ncol(midataset.g2)))
    mi.all = rbind(midataset.g1[midataset.g1$imputation==imputation,1:colcount], midataset.g2[midataset.g2$imputation==imputation,1:colcount])
    main.analysis.gp.diff[m, idx] = get_main_analysis(mi.all[,-c(1,2)])[[2]][1]$coefficients[1]
}

dput(mean.bias.g1, "mean.bias.g1.txt")
dput(mean.bias.g2, "mean.bias.g2.txt")
dput(pointwise.bias.g1, "pointwise.bias.g1.txt")
dput(pointwise.bias.g2, "pointwise.bias.g2.txt")
dput(linear.estimates.g1, "linear.estimates.g1.txt")
dput(linear.estimates.g2,"linear.estimates.g2.txt")
dput(main.analysis.gp.diff, "main.analysis.gp.diff.txt")

m_next = m+1
cat(m_next,"\n", file="simulation_iteration.txt")

}

aggregate= function(df, M=1000){
    out = apply(df,c(2,3), sum, na.rm=T)/M
out
}

#main.analysis.gp.diff
cat("Pointwise bias: group 1")
aggregate(pointwise.bias.g1)
cat("Pointwise bias: group 2")
aggregate(pointwise.bias.g2)

cat("Mean bias: group 1")
aggregate(mean.bias.g1)
cat("Mean bias: group 2")
aggregate(mean.bias.g2)

cat("linear estimates: group 1")
aggregate(linear.estimates.g1)
cat("linear estimates: group 2")
aggregate(linear.estimates.g2)

cat("\n group results")
colMeans(main.analysis.gp.diff, na.rm=T)
