#Alternate hypothesis simulation


setwd("/Users/egoossen/Documents/Research/SEOS/")

source("seos_covariance.R")
source("seos_p_act.R")
source("seos_cdf_sim.R")
source("seos_gumbel_pval.R")
#rm(list = ls())
library(selectiveInference)
library(doParallel)




tot_num = 15
alt_num = 4
nul_num = tot_num - alt_num

iters = 1000
pvals = matrix(NA, nrow = iters, ncol = alt_num)
alt_exp = matrix(NA, nrow = iters, ncol = alt_num)
nul_exp = matrix(NA, nrow = iters, ncol = nul_num)


for(j in 1:iters){
pvals[j,] = 1 - pnorm(rnorm(alt_num , mean = 1.5, sd = 1) )
alt_exp[j,] = -log(pvals[j,])
nul_exp[j,] = rexp(n = nul_num, rate = 1)
tot_exp = c(alt_exp[j,],nul_exp[j,] )
all_ord_exp = tot_exp[order(tot_exp, decreasing = T)]
alt_ord_exp = alt_exp[j,order(alt_exp[j,], decreasing = T)]
nul_ord_exp = nul_exp[j,order(nul_exp[j,], decreasing = T)]
SEOS_real = rep(NA,tot_num)
SEOS_aprx = rep(NA,tot_num)

for(k in 1:tot_num){
  SEOS_real[k] = sum(all_ord_exp[1:k])
  summ =  0
  l = 1
  while(l <= k ){
    if(l <= alt_num){
      summ = summ +   alt_ord_exp[l]
    }else{
      summ = summ +   nul_ord_exp[l - alt_num]
    }
    l = l + 1
  }
  SEOS_aprx[k] = summ
}

sum_ord_stat_exp =  rep(NA,tot_num)
all_SEOS_pval = rep(NA,tot_num)
aprox_SEOS_pval = rep(NA,tot_num)
for(h in 1:length(sum_ord_stat_exp)){
  SEOS_real
  #p_val_ord_stat_exp[h] = p_sum_ord_exp_sim(p-h,sum_ord_stat_exp[h],p, iters_ = 100000)
  all_SEOS_pval[h] = p_sum_ord_exp_sim(i_ = tot_num-h, t_ = SEOS_real[h], n_  = tot_num, iters_ = 1000, sim_type = "importance")
  aprox_SEOS_pval[h] = p_sum_ord_exp_sim(i_ = tot_num-h, t_ = SEOS_aprx[h], n_  = tot_num, iters_ = 1000, sim_type = "importance")
  }


print(c(min(all_SEOS_pval),min(aprox_SEOS_pval)))
}
length(which(pvals > .05))/(2*iters)
mean(alt_exp)


-log(runif(2,0,.1))
alt_rate = .5

nul_exp = rexp(n = nul_num, rate = 1)
alt_exp = rexp(n = alt_num, rate = alt_rate)
tot_exp = c(alt_exp,nul_exp) ; tot_exp

all_ord_exp = tot_exp[order(tot_exp, decreasing = T)]
nul_ord_exp = nul_exp[order(nul_exp, decreasing = T)]

top_5 = sum(all_ord_exp[1:5]) ; top_5
alt_plus_top_one = sum(c(alt_exp, nul_ord_exp[1])) ; alt_plus_top_one

n = 100
sd = 2







