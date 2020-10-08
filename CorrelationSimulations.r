
#simulating a nx10 matrix using different PDFs
matrixcreator <- function(n = 500){
  x1 = rnorm(n, mean= 4, sd = 4)
  x2 = rweibull(n, shape = 3.7, scale= 2.3)
  x3 = rexp(n, rate = 2)
  x4 = rgamma(n, shape=0.8)
  x5 = rexp(n, rate=7)
  x6 = rlogis(n, location = 2, scale= 10)
  x7 = rnorm(n, mean=4, sd = .5)
  x8 = rweibull(n, shape = 1, scale = 0.2)
  x9 = rgamma(n, shape= 4, scale = 1.7)
  x10 = rnorm(n)
  X = matrix(c(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10), ncol= 10)
  return(X)
}

matrixcreator(10)

#now creating a y column that is linearly correlated to the matrix
ycreator = function(X){
  y = 0.5*X[,1] + 0.2*X[,2] - 0.8*X[,3] + rnorm(nrow(X))
  return(y)
}

ycreator(matrixcreator(1000))

#creating n.samps matrices and y values, then checking coverage probabilities as well as false positives and false negatives
sampler = function(n.samps, n = 1000, seed = 1453){
  set.seed(seed)
  ints = matrix(ncol = 10, nrow = n.samps)
  falneg = matrix(ncol=3, nrow = n.samps)
  for (i in 1:n.samps){ 
    X = matrixcreator(n)
    y = ycreator(X)
    mod = lm(y ~ X)
    j = confint(mod)
    in_interval = c(j[2,1]<0.5 & 0.5<j[2,2], j[3,1]<0.2 & 0.2<j[3,2], j[4,1] < -0.8 & -0.8 < j[4,2], apply(j[5:11,], 1, function(x) prod(x) < 0))
    ints[i,] = in_interval
    falneg[i,] = apply(j[2:4,], 1, function(x) prod(x) < 0)
  }
  confint_perc = apply(ints, 2, function(x) sum(x)/n.samps)
  false_positives = 1- confint_perc[4:10]
  false_negatives = apply(falneg, 2, function(x) sum(x)/n.samps)
  
  return(list('confint_perc' = confint_perc, 'false_positives' = false_positives, 'false_negatives' = false_negatives))
}

sampler(500, 10000, 400)

#changing the number of samples and considering how coverage percentage, false negatives, and false positives change as we increase the number of samples
varying_samps <- function(runs, n=500, seed = 500){
  set.seed(seed)
  CI = matrix(nrow= length(runs), ncol = 10)
  FP = matrix(nrow= length(runs), ncol = 7)
  FN = matrix(nrow = length(runs), ncol = 3)
  for (i in 1:length(runs)){
    iteration = sampler(runs[i], n)
    CI[i, ] = iteration$confint_perc
    FP[i, ] = iteration$false_positives
    FN[i, ] = iteration$false_negatives
  }
  return(list('CI' = CI, 'FP' = FP, 'FN' = FN))
}

runs = seq(50, 2050, 200)
varying_samps(runs)

#creating a correlation between x1 and x4 then checking how coverage changes
correlationcreator <- function(rho, n.samps=500, n= 2000, seed=600){
  set.seed(seed)
  ints = matrix(ncol = 10, nrow = n.samps)
  falneg = matrix(ncol=3, nrow = n.samps)
  for (i in 1:n.samps){ 
    X = matrixcreator(n)  
    X[,4] = rho*X[,1] + sqrt(1 - rho^2)*rnorm(n)
    y = ycreator(X)
    mod = lm(y ~ X)
    j = confint(mod)
    in_interval = c(j[2,1]<0.5 & 0.5<j[2,2], j[3,1]<0.2 & 0.2<j[3,2], j[4,1] < -0.8 & -0.8 < j[4,2], apply(j[5:11,], 1, function(x) prod(x) < 0))
    ints[i,] = in_interval
    falneg[i,] = apply(j[2:4,], 1, function(x) prod(x) < 0)
  }
  confint_perc = apply(ints, 2, function(x) sum(x)/n.samps)
  false_positives = 1- confint_perc[4:10]
  false_negatives = apply(falneg, 2, function(x) sum(x)/n.samps)
  
  return(list('confint_perc' = confint_perc, 'false_positives' = false_positives, 'false_negatives' = false_negatives))
}

correlationcreator(0.3)

#seeing how coverage changes at different levels of correlation
varying_rho <- function(rho_sequence, n.samps=500, n =1000, seed = 100){
  set.seed(seed)
  CI = matrix(nrow= length(rho_sequence), ncol = 10)
  FP = matrix(nrow= length(rho_sequence), ncol = 7)
  FN = matrix(nrow = length(rho_sequence), ncol = 3)
  for (i in 1:length(rho_sequence)){
    iteration = correlationcreator(rho_sequence[i], n.samps, n)
    CI[i, ] = iteration$confint_perc
    FP[i, ] = iteration$false_positives
    FN[i, ] = iteration$false_negatives
  }
  return(list('CI' = CI, 'FP' = FP, 'FN' = FN))
}

rho_sequence = seq(0, 0.9, 0.1)
rho_sequence
varying_rho(rho_sequence)

#plotting both coverage over changing samps and rho
runs = seq(50, 4050, 200)
CI_plot = varying_samps(runs)
library(viridis)
col = viridis(3)
par(mfrow=c(1,1), mgp=c(1.5,0,0), tcl=0, mar=c(3,3,1,1), cex.lab=1, cex.axis=1)
plot(runs, CI_plot$CI[,1], type= 'l', xlab = "Number of Samples", ylab = "Confidence Estimate", col = col[1], ylim = c(min(CI_plot$CI[,1:3]) - 0.003, max(CI_plot$CI[,1:3]) + 0.003))
for (i in 2:3) points(runs, CI_plot$CI[,i], col = col[i], type= 'l', lty = i)
abline(h=0.950, lty = 2)
legend('topright', legend = c('x1', 'x2', 'x3'), col = col, lty = 1:3 , cex = 0.8, bty = 'n')
pdf('varying_samps_plot.pdf')
dev.off()

rho_sequence = seq(0, 0.9, 0.05)
rho_plot = varying_rho(rho_sequence, 1000, 2000)
plot(rho_sequence, rho_plot$CI[,1], type = 'l', xlab = 'Rho Values', ylab = 'X1 Confidence Estimates')
abline(h = 0.95, lty = 2)

par(mfrow=c(1,1), mgp=c(1.5,0,0), tcl=1, mar=c(3,3,1,1), cex.lab=1, cex.axis=1)
plot(rho_sequence, rho_plot$FP[,1], type = 'l', xlab = 'Rho Values', ylab = 'X4 False Positive Rate')
abline(h = 0.95, lty = 2)

rho_plot


