library(LaplacesDemon)
setwd("~/Documents/BOSS/R")
source("power_quantiles.R")

mu_alt = 1
n_samples = 10000
alt_pvs = 1- pnorm(alt_test_stat)

alt_test_stat = rnorm(n_samples,mu_alt,1)
neg_log_alt_pvs = -log(alt_pvs)

nuexp_mean = mean(neg_log_alt_pvs); print(nuexp_mean)

alt_pvs = 1- pnorm(alt_test_stat)
neg_log_alt_pvs = -log(alt_pvs)


find_beta_kld = function(beta){
px <- dexp(seq(2,5,.1),1/beta)
py <- pdf_y(seq(2,5,.1),1)
KLD(px,py)$intrinsic.discrepancy 
}

nuexp_mean = optim(4,find_beta_kld, method = "Brent", lower = .1, upper = 10)$par ; print(nuexp_mean)
