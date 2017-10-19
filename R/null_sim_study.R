rm(list = ls())
library(gap)
library(doParallel)

#setwd("/noback/egoossen/phd_research/SEOS/R/functions/BOSS")
setwd("/Users/egoossen/Documents/BOSS/R")
source("null_sim_study_functions.R")
source("exp_covariance.R")
source("boss_cdf_sim.R")
source("core_stat_funcs.R")
start = Sys.time()
simulation_study_null(iterations = 10000,num_params = c(30), unused_cores = 2)#,100,200))
end = Sys.time()
end - start

