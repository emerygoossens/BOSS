rm(list = ls())


setwd("//noback/egoossen/phd_research/BOSS")
load("output/strug_output/strug_analysis_precompute.rda")

source('exp_covariance.R')
source('seos_cdf_sim.R')
source('importance_sampler.R')

library(gap)
library(profvis)
library(mvnfast)
library(doParallel)


library(doParallel)
library(mvtnorm)
params = c(10,20,30,50,100)

p = 100
corr_mat = corrmat_apical_PD[1:p,1:p]

corr_mat = diag(rep(1,p))


max_iter = 10000
inner_iter = 1000
# null_iters = 10000
#null_iter_pvs = rep(NA,null_iters_)

test_stat = rmvnorm(n = 1, sigma = corr_mat)
marginal_pvs = compute_marginal_pvalues(test_stat)


seos_naive_pv = compute_min_SEOS_pv(marginal_pvs, sim_iters = 10000, parallelize = F, compute_ent = F);# null_naive
naive_pv =seos_naive_pv$min_pv
corr_pv = min(naive_pv*log(p),.999)


min_pv_test_stat = naive_pv







maha_is = function(min_pv_test_stat, test_stat,  simul_type = "RB", is_scale = 10, parallelize_is = F, inner_iter = 100, max_iter = 1000, A_mat, cores_ = 12){
  #p = 10; iters = 10; unused_cores = 1 ; emp_cov_mvn = Posdef(10)
  #test_stat = test_stat_apical ; emp_cov_mvn = null_test_cov ; unused_cores = 1; cores = 23 ;
  # min_pv = seos_naive_apical_pv$min_pv; test_stat = test_stat_apical;
  # emp_cov_mvn = diag(diag(corrmat_apical)); simul_type = "RB";is_scale = 6
  # parallelize_is = F; inner_iter = 2; max_iter = 3; A_mat = chol_apical
  # is_scale = is_scale_
  
  p = length(test_stat)
  sign_test_stat = sign(test_stat)
  
  
  zero_mean = rep(0,p)
  is_mean = rep(0,p)
  
  registerDoParallel(cores = 8)
  
  
  
  simulations = foreach(icount(max_iter), .combine = rbind) %dopar% {
    test_stat_sim = rmvn(n = 1, mu = is_mean, sigma = A_mat, isChol = T)
    
    pvs_sim = compute_marginal_pvalues(test_stat_sim)
    min_pv_sim = compute_min_SEOS_pv_c(pvs_sim, sim_iters = inner_iter, compute_ent = F, parallelize = F)$min_pv
    
    if(min_pv_sim <min_pv_test_stat){
      is_weight = 1 
      
    }else{
      is_weight = 0 
    }
    
    c(min_pv_sim,is_weight)
    
  }
  sum(simulations[,2])/max_iter
  min_is_mc_pv = min(sum(simulations[which(simulations[,1]<min_pv_test_stat),2])/qq,.999)
  list(min_is_mc_pv = min_is_mc_pv, simulations = simulations)
}

corr_pv_is_rb_mc_weight = function(min_pv_test_stat, test_stat,  simul_type = "RB", is_scale = 10, parallelize_is = F, inner_iter = 100, max_iter = 1000, A_mat, cores_ = 12){
  #p = 10; iters = 10; unused_cores = 1 ; emp_cov_mvn = Posdef(10)
  #test_stat = test_stat_apical ; emp_cov_mvn = null_test_cov ; unused_cores = 1; cores = 23 ;
  # min_pv = seos_naive_apical_pv$min_pv; test_stat = test_stat_apical;
  # emp_cov_mvn = diag(diag(corrmat_apical)); simul_type = "RB";is_scale = 6
  # parallelize_is = F; inner_iter = 2; max_iter = 3; A_mat = chol_apical
  # is_scale = is_scale_

  p = length(test_stat)
  sign_test_stat = sign(test_stat)


  zero_mean = rep(0,p)
  is_mean = rep(0,p)

  registerDoParallel(cores = 8)



    simulations = foreach(icount(max_iter), .combine = rbind) %dopar% {
      test_stat_sim = rmvn(n = 1, mu = is_mean, sigma = A_mat, isChol = T)

      pvs_sim = compute_marginal_pvalues(test_stat_sim)
      min_pv_sim = compute_min_SEOS_pv_c(pvs_sim, sim_iters = inner_iter, compute_ent = F, parallelize = F)$min_pv

      if(min_pv_sim <min_pv_test_stat){
       is_weight = 1 
 
      }else{
        is_weight = 0 
      }

      c(min_pv_sim,is_weight)

    }
    sum(simulations[,2])/max_iter
  min_is_mc_pv = min(sum(simulations[which(simulations[,1]<min_pv_test_stat),2])/qq,.999)
  list(min_is_mc_pv = min_is_mc_pv, simulations = simulations)
}
