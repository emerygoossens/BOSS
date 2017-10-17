#Monte Carlo Analysis
setwd("~/Documents/BOSS/R/")
source("exp_covariance.R")
source("boss_cdf_sim.R")
library(Hmisc)

#First make function to simulate each set of order combination via monte carlo separately.


r = 10 # Total number of p-values
k = 4  # Number to include in Ordered Subset
n_iters = 100 # Number of iterations
t = 2 # We find the probability for which Ordered subset is greate than t

test_statistic =  sum_ordered_subset(sort(rexp(r,1), decreasing = T))

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

sim_out = simulate_sum_ordered_exp(10,100000)

mc_output = mc_sum_ord_exp(test_stat = test_statistic, sim_out)

mc_indicators = mc_output$indicators
mc_estimates = mc_output$mc_estimate ; mc_estimates
mc_var = mc_output$mc_var ; mc_var
mc_lq =  mc_output$lower_quantile ; mc_lq
mc_uq =  mc_output$upper_quantile ;mc_uq


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

rbmc_output = rbmc_sum_ord_exp(test_statistic)
rbmc_indicators = rbmc_output$estimates
rbmc_estimates = rbmc_output$rbmc_estimate ;rbmc_estimates
rbmc_var = rbmc_output$rbmc_var ; rbmc_var
rbmc_lq =  rbmc_output$lower_quantile ; rbmc_lq
rbmc_uq =  rbmc_output$upper_quantile ;rbmc_uq


true_output = rbmc_sum_ord_exp(test_statistic, 100000)
true_indicators = true_output$estimates
true_estimates = true_output$rbmc_estimate ;true_estimates
true_var = true_output$rbmc_var ; true_var
true_lq =  true_output$lower_quantile ; true_lq
true_uq =  true_output$upper_quantile ;true_uq



ind_mc = 1:10-1/3 
ind_rbmc =1:10

ind_true = 1:10+1/3
  
errbar(ind_mc, mc_estimates, mc_estimates - 2*mc_var, mc_estimates + 2*mc_var, col = "blue", lty = 2, xlim = c(0.5, 10.5))
errbar(ind_rbmc, rbmc_estimates, rbmc_estimates - 2*rbmc_var, rbmc_estimates + 2*rbmc_var, add = T, col = "red", lty = 2)
errbar(ind_true, true_estimates, true_estimates - 2*true_var, true_estimates + 2*true_var, add = T, col = "green", lty = 2) 
       



