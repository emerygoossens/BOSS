



#BOSS Monte Carlo Analysis
setwd("~/Documents/BOSS/R/")
source("exp_covariance.R")
source("boss_cdf_sim.R")
library(Hmisc)
library(doParallel)


registerDoParallel(7)
simulate_sum_ordered_exp = function(r, n_iters){
  
  output = matrix(NA,nrow = n_iters, ncol = r) #store observed values
  for( n in 1:n_iters){
    ordered_exps = sort(rexp(r,1), decreasing = T)
    sum_ordered_exps = sum_ordered_subset(ordered_exps)
    output[n,] = sum_ordered_exps 
  }
  output  
}

mc_sum_ord_exp = function(test_stat,sim_exp){
  if( missing(sim_exp)){
    sim_exp = simulate_sum_ordered_exp(r = length(test_stat), n_iters = 1000)
  }
  r = length(test_stat)
  mc_estimate= rep(NA, r)
  mc_var= rep(NA, r)
  lower_quantile = rep(NA, r)
  upper_quantile = rep(NA, r)
  n_iters = nrow(sim_exp)
  indicators = matrix(0,nrow = n_iters, ncol = r)
  for(k in 1:length(test_stat)){
    one_inds = which(sim_exp[,k]>test_stat[k])
    indicators[one_inds,k] = 1
    mc_estimate[k] = length(one_inds)/n_iters
    mc_var[k] =  var(indicators[,k])
    lower_quantile[k] = quantile(indicators[,k], probs = c(.05))
    upper_quantile[k] = quantile(indicators[,k], probs = c(.95))
  }
  list(mc_estimate = mc_estimate, mc_var = mc_var,lower_quantile = lower_quantile, upper_quantile = upper_quantile, indicators = indicators)
}

mc_estimate_parallel = function(test_stat, n_iters){
  r = length(test_stat)
  cores_ = detectCores() - 1
  registerDoParallel(cores = cores_)
  summation = foreach(icount(n_iters), .combine = "+") %dopar% {
    ordered_exps = sort(rexp(r,1), decreasing = T)
    sum_ordered_exps = sum_ordered_subset(ordered_exps)
    sum_ordered_exps > test_stat
  }
summation/n_iters
}


rbmc_sum_ord_exp = function(test_stat, iters_ = 1000){
  
  r = length(test_stat)
  cores_ = detectCores() - 1
  registerDoParallel(cores = cores_)
  estimates = foreach(icount(iters_), .combine = rbind) %dopar% {
    #for(k in 1:iters_){
    scaled_exp_rvs = c(rexp((r -1),rate = 1)*c(1/((2:r))),0)
    sum_scaled_exp_rvs = rep(NA, r)
    sum_scaled_exp_rvs[r] = scaled_exp_rvs[r]
    for(j in (r-1):1){
      sum_scaled_exp_rvs[j] = sum_scaled_exp_rvs[j+1] + scaled_exp_rvs[j]
    }
    
    #summation = summation  + pgamma(t_ - (r - (r - 1):0)*sum_scaled_exp_rvs, shape = 1:r, rate = 1, lower.tail = F )
    output = pgamma(test_stat - (r - (r - 1):0)*sum_scaled_exp_rvs, shape = 1:r, rate = 1, lower.tail = F )
    output
  }
  lower_quantile = rep(NA, r)
  upper_quantile = rep(NA, r)
  rbmc_estimate= rep(NA, r)
  rbmc_var= rep(NA, r)
  for(k in 1:length(test_stat)){
    rbmc_estimate[k] = mean(estimates[,k])
    rbmc_var[k] =  var(estimates[,k])
    lower_quantile[k] = quantile(estimates[,k], probs = c(.05))
    upper_quantile[k] = quantile(estimates[,k], probs = c(.95))
  }
  list(rbmc_estimate = rbmc_estimate, rbmc_var = rbmc_var, lower_quantile = lower_quantile, upper_quantile = upper_quantile, estimates = estimates)
}



#First make function to simulate each set of order combination via monte carlo separately.


r = 10 # Total number of p-values
k = 4  # Number non-null parameters
n_iters = 100 # Number of iterations
t = 2 # We find the probability for which Ordered subset is greate than t

outer_iters = 10

test_statistic =  sum_ordered_subset(sort(c(rexp(k,1/3), rexp(r-k,1)), decreasing = T)) ; test_statistic


boss_output = matrix(NA,nrow = outer_iters, ncol = 6)

for(b in 1:outer_iters){
print(b)
print("MC")
mc_estimates = mc_estimate_parallel(test_stat = test_statistic, 100000)
boss_output[b,1] = min(mc_estimates)
boss_output[b,4] = which(mc_estimates == min(mc_estimates))[1]
print(mc_estimates)
print("RBMC Some Iters")

some_iters = p_sum_ord_exp_sim_simul(t_ = test_statistic, iters_ = 1000, paralellize = T)$ord_comb_pvs
boss_output[b,2] = min(some_iters)
boss_output[b,5] = which(some_iters == min(some_iters))[1]
print(some_iters)

print("RBMC More Iters")

more_iters = p_sum_ord_exp_sim_simul(t_ = test_statistic, iters_ = 10000, paralellize = T)$ord_comb_pvs
boss_output[b,3] = min(more_iters)
boss_output[b,6] = which(more_iters == min(more_iters))[1]
print(more_iters)
}

boss_output




ind_mc = 1:10-1/3 
ind_rbmc =1:10

ind_true = 1:10+1/3

mean_mc = mean(boss_output[,1])
var_mc = var(boss_output[,1])

mean_some = mean(boss_output[,2])
var_some = var(boss_output[,2])

mean_more = mean(boss_output[,3])
var_more = var(boss_output[,3])


df = cbind(boss_output[,1],rep(1,outer_iters))
for(j in 2:3){
 df = rbind(df,  cbind(boss_output[,j],rep(j,outer_iters)))
}
df = data.frame(df)

# 
# errbar(1, mean_mc , mean_mc - 2*var_mc, mean_mc + 2*var_mc, col = "blue", lty = 2, xlim = c(1,3))
# errbar(2, mean_some, mean_some - 2*var_some, mean_some + 2*var_some, col = "red", lty = 2, add = T)
# errbar(3, mean_more, mean_more - 2*var_more, mean_more + 2*var_more, col = "green", lty = 2, add = T)
# 


boxplot(X1~X2, data= df)

