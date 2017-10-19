#core statistcs functions


# compute_marginal_pvalues = function(test_statistic, cov_mat){
#   2*pnorm(abs(test_statistic), sd = diag(cov_mat),  lower.tail = F, log.p = F)
# }


correct_zero_combo_pvalues = function(combo_pvals){
  ifelse(combo_pvals == 0,min(ifelse(combo_pvals>0,combo_pvals,1))*2, combo_pvals)
}

compute_marginal_pvalues = function(test_statistic, sigma = rep(1,length(test_statistic))){
  2*pnorm(abs(test_statistic), sd = sigma,  lower.tail = F, log.p = F)
}


fishers_method = function(ps){
  pchisq(sum(-2*log(ps)), df = 2*length(ps), lower.tail = F)
}

fishers_for_all_ordered_subsets = function(pvalues){
  pvalues = pvalues[order(pvalues,decreasing = F)]
  num_ps = length(pvalues)
  combo_pvs = rep(NA,num_ps)
  for(j in 1:num_ps){
    combo_pvs[j] = fishers_method(pvalues[1:j])
  }
  combo_pvs
}

simulate_lm_pvalues = function(p = 10, mean = rep(0,p), sigma = diag(rep(1,p)), beta = rep(0,p),n = 1000){
  x=matrix(rmvnorm(n = n, mean = mean, sigma =  sigma),n,p, byrow = T)
  x=scale(x,T,F)
  y=x%*%beta+rnorm(n)
  y=y-mean(y)
  
  univariate_pvals = matrix(NA, nrow = 1, ncol = p)
  
  for(j in 1:p){
    univariate_pvals[,j] = coef(summary(lm(y ~ 0 + x[,j])))[1,4]
    
  }
  list(univariate_pvals = univariate_pvals, x = x, y = y)
}

simulate_mvn_pvalues = function(p = 10, mean = rep(0,p), sigma = diag(rep(1,p)),n = 1000){
  test_stat = rmvnorm(n,mean = mean, sigma = sigma)
  univariate_pvals = compute_marginal_pvalues(test_stat)
  list(univariate_pvals = univariate_pvals, test_stat = test_stat)
}


