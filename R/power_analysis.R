#Approximate Power Analysis
setwd("/Users/egoossen/Documents/BOSS/R")
source("core_stat_funcs.R")
source("exp_covariance.R")
source("boss_cdf_sim.R")
library(doParallel)
registerDoParallel(cores = detectCores()-1)
#Suppose we have a single test statistics with mean 1, std 1

n_samples = 1000
n_iters = 5000

r = 2
num_nonnull = c(1,2)

alpha = .05

mu_alt_list = c(.5,1,1.5,2,3,4)
alt_output = array(NA, dim= c(n_iters, length(mu_alt_list),  length(num_nonnull) ))
apx_output = array(NA, dim= c(n_iters, length(mu_alt_list) ,  length(num_nonnull)))

alt_power_out = matrix(NA, nrow = length(mu_alt_list), ncol = length(num_nonnull))
apx_power_out = matrix(NA, nrow = length(mu_alt_list), ncol = length(num_nonnull))
ent = boss.effective.num.tests(r)$boss_ent

for(m in 1:length(mu_alt_list)){
  mu_alt = mu_alt_list[m]
  #Simulated mean match, can be done theoretically. 
  
  alt_test_stat = rnorm(n_samples,mu_alt,1)
  
  #We now compute the 2-sided p-value
  
  
  # alt_pvs = compute_marginal_pvalues(alt_test_stat)
  # alt_pvs = 2*(1 - pnorm(abs(alt_test_stat)))
  alt_pvs = 1- pnorm(alt_test_stat)
  neg_log_alt_pvs = -log(alt_pvs)
  
  
  #The mean and variance of these the -log(pvs) is given by
  print("Mean of Neg Log Alt PVS")
  mean_neg_log_alt = mean(neg_log_alt_pvs)
  print(mean_neg_log_alt)
  print("Sqrt Var of Neg Log Alt PVS")
  sqrt_var_neg_log_alt = sqrt(var(neg_log_alt_pvs))
  print(sqrt_var_neg_log_alt)
  
  # #If we were to moment match the mean, the variance of the non-unit exponential would be
  # print("First Moment")
  # nuexp_mean = mean(neg_log_alt_pvs); print(nuexp_mean)
  # print("Sqrt Variance")
  # nuexp_sqrt_var = sqrt(var(neg_log_alt_pvs)) ;print(nuexp_sqrt_var)
  # print("Mean of First Moment and Sqrt Variance")
  # nuexp_mean = mean(c(sqrt(var(neg_log_alt_pvs)), nuexp_mean)) ; print(nuexp_mean)
  
  
  find_beta_kld = function(beta){
    px <- dexp(seq(-log(alpha),10,.1),1/beta)
    py <- pdf_y(seq(-log(alpha),10,.1),mu_alt)
    KLD(px,py)$intrinsic.discrepancy 
  }
  
  
  pdf_y = function(y, sqrtn_delta=mu_alt){
    zexpnegy = qnorm(1-exp(-y))
    (dnorm(zexpnegy - sqrtn_delta)/dnorm(zexpnegy))*exp(-y)
  }

  print("Mean of Exp")
  nuexp_mean = optim(sqrt_var_neg_log_alt,find_beta_kld, method = "Brent", lower =sqrt_var_neg_log_alt/2 , upper = 4*mean_neg_log_alt)$par  ; print(nuexp_mean)

  beta_dist = function(x,beta = nuexp_mean){
    exp(-x/beta)/beta
  }

  curve(pdf_y,0,5)
  curve(beta_dist,0,5, add=T)
  
for(k in num_nonnull){

#which is much higher.
# nuexp_mean = nuexp_sqrt_var; print(nuexp_mean)
#Let's compare the power of simulated null and alternative p-values vs. exponentials when using BOSS

alt_apx_power = foreach(icount(n_iters), .combine = rbind) %dopar% {
# for(j in 1:n_iters){

null_pvs = runif(r-k, 0,1)
alt_test_stats = rnorm(k,mu_alt,1)
alt_pvs = compute_marginal_pvalues(alt_test_stats)

all_pvs = c(null_pvs, alt_pvs)

neg_log_alt_pvs = -log(all_pvs)
exp_approx_neg_log_pvs = c(-log(null_pvs),rexp(k, rate = 1/nuexp_mean))
sort_neg_log_alt_pvs = sort(neg_log_alt_pvs,decreasing = T)
sort_neg_log_apx_pvs = sort(exp_approx_neg_log_pvs, decreasing = T)

min_alt_boss = min(ord_exp_to_boss(sort_neg_log_alt_pvs))
min_apx_boss = min(ord_exp_to_boss(sort_neg_log_apx_pvs))
c(min_alt_boss, min_apx_boss)
# alt_output[j,m,k] = min_alt_boss
# apx_output[j,m,k] = min_apx_boss
}
alt_output[,m,k] = alt_apx_power[,1]
apx_output[,m,k] = alt_apx_power[,2]

print("Mean of Alternative:")
print(mu_alt)
print("Non-Null Parameters:")
print(num_nonnull[k])
print("True Power")
alt_power = length(which(alt_apx_power[,1]<(.05/ent)))/n_iters; print(alt_power)
print("Approximate Power")
apx_power = length(which(alt_apx_power[,2]<(.05/ent)))/n_iters; print(apx_power)

alt_power_out[m,k] = alt_power
apx_power_out[m,k] = apx_power

}
}

save(alt_output, file = paste("../Data/Power/alt_power_pvs_",n_iters,"_iters",length(num_nonnull),"_nparams_",Sys.Date(),".Rda",sep = ""))

save(apx_output, file = paste("../Data/Power/apx_power_pvs_",n_iters,"_iters",length(num_nonnull),"_nparams_",Sys.Date(),".Rda",sep = ""))


save(alt_power_out, file = paste("../Data/Power/alt_power_mat_",n_iters,"_iters",length(num_nonnull),"_nparams_",Sys.Date(),".Rda",sep = ""))
save(apx_power_out, file = paste("../Data/Power/apx_power_mat_",n_iters,"_iters",length(num_nonnull),"_nparams_",Sys.Date(),".Rda",sep = ""))

























