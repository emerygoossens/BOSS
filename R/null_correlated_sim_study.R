rm(list = ls())



# setwd("//noback/egoossen/phd_research/BOSS")
setwd("/Users/egoossen/Documents/BOSS/R")
load("strug_analysis_precompute.rda")

source('exp_covariance.R')
source('boss_cdf_sim.R')
source('importance_sampler.R')
source("core_stat_funcs.R")


library(gap)
library(profvis)
library(mvnfast)
library(doParallel)
library(reliaR)

library(doParallel)
library(mvtnorm)
params = c(10,20,30,50,100)

for(q in params){
p = q
corr_mat = corrmat_apical_PD[1:p,1:p]



iterss = 1000
inner_iterss = 100
is_scale_ = round(sqrt(p))
null_iters_ = 100
#null_iter_pvs = rep(NA,null_iters_)
registerDoParallel(cores = detectCores()-2)
null_iter_pvs = foreach(icount(null_iters_), .combine = rbind) %dopar% {

  
  
# corr_mat = diag(diag(corr_mat))
  
chol_mat = chol(corr_mat)
num_sig = 10
test_stat = rmvnorm(n = 1, c(rep(3,num_sig),rep(0,p -num_sig)),sigma = corr_mat)  
marginal_pvs = compute_marginal_pvalues(test_stat,diag(corr_mat))

#marginal_pvs = runif(10,0,.1)

boss_naive_pv = compute_min_boss_pv(marginal_pvs, sim_iters = 10000, parallelize = T, compute_ent = T);# null_naive
boss_naive_pv
contrib_ind = as.numeric(names(boss_naive_pv$var_contributing_pvs))

min_pv_mc_rbmc = corr_pv_mc_rbmc(boss_naive_pv$min_pv, A_mat = chol_mat, max_iter = 10000, cores_ = 6)
min_pv_mc_rbmc$mc_min_boss_pv
min_pv_is_null_test = corr_pv_is_rb_mc(min_pv = boss_naive_pv$min_pv, contrib_ind = contrib_ind,
                                   test_stat = test_stat, 
                                   emp_cov_mvn = corr_mat,  
                                   simul_type = "RB",
                                   is_scale = p, parallelize_is = T,
                                   inner_iter = inner_iterss, 
                                   max_iter = 100000, A_mat = chol_mat, cores_ = 6)
min_pv_is_null_test$min_is_rbmc_pv
}

save(null_iter_pvs,file = paste("output/null_simulation_study_correlated/null_corr_sims_",null_iters_,"_p_",p,"_iters",iterss,"_is_scale_",is_scale_,"_inner_iter_",inner_iterss,"_",Sys.time(),".Rda", sep = ""))


}

