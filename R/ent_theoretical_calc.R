library(reliaR)
library(doParallel)
#Suppose we have the sum (X) and max (Y) of two exponential random variables
#Want to verify the probability of a type I error under the null
#
source('~/Documents/BOSS/R/exp_covariance.R')



theoretical_ent = function(ent2,alpha = alpha){
  
  alpha_corrected = alpha/ent2
  #alpha_ = 0.1/2
  #The threshold for the sum of two exponentials would be
  
  sum_thresh = qgamma(p = 1-alpha_corrected,shape = 2, rate = 1, lower.tail = T)
  
  #The threshold for the max would be
  
  max_thresh = qgen.exp(p = 1-alpha_corrected, alpha = 2, lambda = 1, lower.tail = T)
  
  #We are then interested in the 
  #P(min(BOSS(X), BOSS(Y)) < alpha_corrected) = 1 - P(X > sum_thresh, Y > max_thresh)
  
  # 1- 3.7*exp(-3.7) - exp(-3.7)
  
  geg1 = function(z){
    1 - exp(-z)*(z + 1)
  }
  
  geg2 = function(z,w){
    exp(-z)*(z-1) + exp(-w)*(w - 2*z - 1)
  }
  
  (1 - geg1(max_thresh) - geg2(max_thresh,sum_thresh)) - alpha
  
}

find_ent2 = function(alpha = 0.05){
uniroot(theoretical_ent,c(1,2), alpha = alpha)$root
}

# ent2 - boss.effective.num.tests(2)$li
