
#Installation of the following packages may be required. 
install.packages("inline")
install.packages("compiler")
install.packages("gap")
install.packages("reliaR")
install.packages("doParallel")

library(inline)
library(compiler)
library(gap)
library(reliaR)
library(doParallel)

setwd("~/Downloads/Sciome/")
source("exp_covariance.R")
source("boss_cdf_sim.R")

p_values = c(exp(-rexp(4,1/15)),runif(16))
names(p_values) = 1:20

compute_min_boss_pv(p_values)

