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
# is_scale_ = round(sqrt(p))
null_iters_ = 100
#null_iter_pvs = rep(NA,null_iters_)
registerDoParallel(cores = detectCores()-2)
null_iter_pvs = foreach(icount(null_iters_), .combine = rbind) %dopar% {

  
  
# corr_mat = diag(diag(corr_mat))
  
chol_mat = chol(corr_mat)
num_sig = 3
test_stat = rmvnorm(n = 1, c(rep(3,num_sig),rep(0,p -num_sig)),sigma = corr_mat)  


is_weight = exp(dmvn(X = test_stat, mu = rep(0,p), sigma = chol_mat, isChol = T, 
                     log = T) - dmvn(X = test_stat, mu = rep(0,p), sigma = 1.3*chol_mat, isChol = T, log =T))

marginal_pvs = compute_marginal_pvalues(test_stat,diag(corr_mat))

#marginal_pvs = runif(10,0,.1)

boss_naive_pv = compute_min_boss_pv(marginal_pvs, sim_iters = 10000, parallelize = T, compute_ent = T);# null_naive
boss_naive_pv

start1 = Sys.time()
ar_pv = ar_is(boss_naive_pv$min_pv, A_mat = chol_mat,inner_iter = inner_iterss, max_iter = 5000, cores_ = 6, corr_mat = corr_mat)
end1 = Sys.time()
simulations = ar_pv$simulations


contrib_ind = as.numeric(names(boss_naive_pv$var_contributing_pvs))

start2 = Sys.time()
min_pv_mc_rbmc = corr_pv_mc_rbmc(boss_naive_pv$min_pv,inner_iter = inner_iterss, A_mat = chol_mat, max_iter = 10000, cores_ = 6)
end2 = Sys.time()


start3 = Sys.time()
is_var_pv= is_var_scale(is_scale =2,inner_iter = inner_iterss, boss_naive_pv$min_pv, A_mat = chol_mat, max_iter = 10000, cores_ = 6)
end3 = Sys.time()


start4 = Sys.time()
ar_is_var_pv= ar_is_scale_var(var_scale = 2,inner_iter = inner_iterss,corr_mat = corr_mat, boss_naive_pv$min_pv, A_mat = chol_mat, max_iter = 10000, cores_ = 6)
end4 = Sys.time()


print(end1 - start1)
print(end2 - start2)
print(end3 - start3)
print(end4 - start4)


ar_pv$mc_min_boss_pv
min_pv_mc_rbmc$mc_min_boss_pv
is_var_pv$min_is_rbmc_pv
ar_is_var_pv$min_is_rbmc_pv


length(which(ar_pv$simulations[,1]<boss_naive_pv$min_pv))
length(which(min_pv_mc_rbmc$simulations[,1]<boss_naive_pv$min_pv))
length(which(is_var_pv$simulations[,1]<boss_naive_pv$min_pv))
length(which(ar_is_var_pv$simulations[,1]<boss_naive_pv$min_pv))



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

