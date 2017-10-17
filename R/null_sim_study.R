rm(list = ls())
library(gap)
library(doParallel)

#setwd("/noback/egoossen/phd_research/SEOS/R/functions/BOSS")
setwd("/Users/egoossen/Documents/BOSS/R")
source("null_sim_study_functions.R")
source("exp_covariance.R")
source("boss_cdf_sim.R")

start = Sys.time()
simulation_study_null(iterations = 100000,num_params = c(10,20,30,50), unused_cores = 2)#,100,200))
end = Sys.time()
end - start

