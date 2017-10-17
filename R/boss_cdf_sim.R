

cpp_if_src <- '
  Rcpp::NumericVector xa(a);
  int n_xa = xa.size();
  for(int i=0; i < n_xa; i++) {
    if(xa[i]<0) xa[i] = 0;
  }
  return xa;
'
cpp_if <- cxxfunction(signature(a="numeric"), cpp_if_src, plugin="Rcpp")


relative_error = function(estimate, value){
    abs(1 - estimate/value)
}


unif_ord_corrected_pvs = function(p_values){
	number_pvs = length(p_values)
	for(j in 1:number_pvs){

	}
}

correct_zero_combo_pvalues = function(combo_pvals){
  ifelse(combo_pvals == 0,min(ifelse(combo_pvals>0,combo_pvals,1))*2, combo_pvals)
}


pvalues_to_boss = function(univariate_pvs){
  
  ord_stat_pvals = sort(univariate_pvs)
  ord_stat_exp = -log(ord_stat_pvals)
  sum_ord_stat_exp =  rep(NA,length(univariate_pvs))
  sum_ord_stat_exp[1] = ord_stat_exp[1]
  for(j in 2:length(univariate_pvs)){
    sum_ord_stat_exp[j] = ord_stat_exp[j] + sum_ord_stat_exp[j-1]
  }

  p_sum_ord_exp_sim_simul(sum_ord_stat_exp, iters = 10000)$ord_comb_pvs

  
}

sum_ordered_subset = function(ordered_subset){
  sum_ordered_subsets =  rep(NA,length(ordered_subset))
  sum_ordered_subsets[1] = ordered_subset[1]
  for(j in 2:length(ordered_subset)){
    sum_ordered_subsets[j] = ordered_subset[j] + sum_ordered_subsets[j-1]
  }
  sum_ordered_subsets
}

ord_exp_to_boss = function(ord_stat_exp){
  sum_ord_stat_exp =  rep(NA,length(ord_stat_exp))
  sum_ord_stat_exp[1] = ord_stat_exp[1]
  for(j in 2:length(ord_stat_exp)){
    sum_ord_stat_exp[j] = ord_stat_exp[j] + sum_ord_stat_exp[j-1]
  }
  p_sum_ord_exp_sim_simul(sum_ord_stat_exp, iters = 1000)$ord_comb_pvs
  

}



p_sum_ord_exp_sim_simul = function( t_ ,lambda_ = 1, iters_ = 1000, sim_type= "RB", paralellize = F){ #, sim_type = "importance"){

summation = 0
n_ = length(t_)
if(paralellize){
  cores_ = detectCores() - 1
  registerDoParallel(cores = cores_)
  summation = foreach(icount(iters_), .combine = "+") %dopar% {
  #for(k in 1:iters_){
    scaled_exp_rvs = c(rexp((n_ -1),rate = 1/lambda_)*c(1/((2:n_))),0)
    sum_scaled_exp_rvs = rep(NA, n_)
    sum_scaled_exp_rvs[n_] = scaled_exp_rvs[n_]
    for(j in (n_-1):1){
      sum_scaled_exp_rvs[j] = sum_scaled_exp_rvs[j+1] + scaled_exp_rvs[j]
    }
    
    #summation = summation  + pgamma(t_ - (n_ - (n_ - 1):0)*sum_scaled_exp_rvs, shape = 1:n_, rate = 1, lower.tail = F )
    output = pgamma(t_ - (n_ - (n_ - 1):0)*sum_scaled_exp_rvs, shape = 1:n_, rate = 1, lower.tail = F )
    output
  }
  #summation = colSums(simulations)
  
}else{
for(k in 1:iters_){
scaled_exp_rvs = c(rexp((n_ -1),rate = 1/lambda_)*c(1/((2:n_))),0)
sum_scaled_exp_rvs = rep(NA, n_)
sum_scaled_exp_rvs[n_] = scaled_exp_rvs[n_]
for(j in (n_-1):1){
 sum_scaled_exp_rvs[j] = sum_scaled_exp_rvs[j+1] + scaled_exp_rvs[j]
}

summation = summation  + pgamma(t_ - (n_ - (n_ - 1):0)*sum_scaled_exp_rvs, shape = 1:n_, rate = 1, lower.tail = F )
}
}
ord_comb_pvs = summation/iters_
true_pv_min = (1-(1-exp(-t_[1]))^n_)
true_pv_fishers = pgamma(t_[n_], shape = n_ , rate = 1, lower.tail = F)
pv_min_rel_err = relative_error(ord_comb_pvs[1], true_pv_min)
pv_fishers_rel_err = relative_error(ord_comb_pvs[n_], true_pv_fishers)

list(ord_comb_pvs = ord_comb_pvs, true_pv_min = true_pv_min, true_pv_fishers = true_pv_fishers, pv_min_rel_err = pv_min_rel_err, pv_fishers_rel_err = pv_fishers_rel_err)

}


p_sum_ord_exp_sim_simul_c <- cmpfun(p_sum_ord_exp_sim_simul)


compute_min_boss_pv = function(pvalues, sim_iters = 20000, compute_ent = T, parallelize = F, sim_type = "RB"){
  #pvalues = pv_apical; sim_iters = 2
  length_pv = length(pvalues)
  #names(pvalues) = 1:length_pv
  neg_log_ord_pvs = -log(sort(pvalues, decreasing = F))


  ordered_combo_output = compute_sum_combo_pvs_c(neg_log_ord_pvs, type = "boss", iterations = sim_iters, sim_type = sim_type, parallelize = parallelize)
  ordered_combo = correct_zero_combo_pvalues(ordered_combo_output$ord_comb_pvs)
  rel_error_fisher  = ordered_combo_output$pv_fishers_rel_err
  rel_error_min = ordered_combo_output$pv_min_rel_err
  true_fisher = ordered_combo_output$true_pv_fishers
  true_pv_min = ordered_combo_output$true_pv_min
  min_pv = min(ordered_combo)[1]
  min_pv_ind = min(which(ordered_combo == min_pv))
  if(compute_ent){
  effective_num_tests = boss.effective.num.tests(length(pvalues))$boss_ent
  }else{
    effective_num_tests = 1
  }
  corr_min_pv = effective_num_tests*min_pv
  vars_contributing = names(neg_log_ord_pvs)[1:min_pv_ind]
  contributing_pvs = pvalues[vars_contributing]
  list(corr_min_pv = corr_min_pv, contributing_pvs = contributing_pvs,
    effective_num_tests = effective_num_tests, min_pv = min_pv,
    ordered_combo = ordered_combo, rel_error_fisher = rel_error_fisher,
    rel_error_min = rel_error_min, true_fisher = true_fisher, true_pv_min = true_pv_min)
}

compute_sum_combo_pvs = function(ordered_test_stats, type_test = "boss", sim_type = "RB",iterations = 10000, parallelize = T){
  num_tests  = length(ordered_test_stats)
  pvs_of_sum = rep(NA,num_tests)
  sum_test_stats = rep(NA,num_tests)
  sum_test_stats[1]  = ordered_test_stats[1]
  for(h in 1:(num_tests-1)){
    sum_test_stats[h+1] = sum_test_stats[h] + ordered_test_stats[h+1]
  }
  # if(parallelize){
  # 
  # for(h in 1:num_tests){
  #   if(type_test == "boss"){
  #   #pvs_of_sum[h] = p_sum_ord_exp_sim(i_ = num_tests-h, t_ = sum_test_stats[h], n_  = num_tests, iters = iterations, sim_type = sim_type)
  #   }
  #   if(type_test == "SES"){
  #     #pvs_of_sum[h] = 1 - pchisq(sum_test_stats[h], df = 2*h) # Chisquare Fisher's method
  #     pvs_of_sum[h] = 1- pgamma(sum_test_stats[h], shape = h, rate = 1 )
  #   }
  # }
  # }else{
    pvs_of_sum = p_sum_ord_exp_sim_simul_c(t_ = sum_test_stats, iters_ = iterations, sim_type= "RB", paralellize = parallelize)

  # }
  pvs_of_sum
}


compute_sum_combo_pvs_c = cmpfun(compute_sum_combo_pvs)


compute_min_boss_pv_c = cmpfun(compute_min_boss_pv)



#Simulate order of exponentials

#t_ = seq(1,2,.5) - top_boss_niid

max_exp_niid_dist_fun = function(lambdas_, thresh){
  1- prod((1-exp(-lambdas_*thresh)))
}

max_exp_spacing_niid_dist_fun = function(lambdas_, thresh){
  1- prod((1-exp(-lambdas_*thresh)))
}

simulate_order = function(lambdas = seq(1,6,2)){
	  num_lambda = length(lambdas)
  	D_vec = rep(NA,length(lambdas))
  	lambdas_copy = lambdas
  	for(j in 1:num_lambda){
    	probs = lambdas_copy/sum(lambdas_copy); #print(probs)
    	D_vec[j] = which(rmultinom(1,1,probs) == 1)
    	lambdas_copy[D_vec[j]] = 0
  	}
  	D_vec
  }

compute_weighted_boss_pv = function(pvs, weights,iterations = 100){
#pvs = pv_trial; weights = weights_trial; iterations = 100

num_lambda  = length(weights)
weights = weights/sum(weights)
exp_lambda = 1/weights

pv_estimates = rep(0,num_lambda)

sorted_weighted_neg_log_pvs = sort(-log(pvs)*weights, decreasing = T)
boss_weighted_neg_log_pvs = rep(NA,num_lambda)
boss_weighted_neg_log_pvs[1] = sorted_weighted_neg_log_pvs[1]
for(q in 2:num_lambda){
  boss_weighted_neg_log_pvs[q] = boss_weighted_neg_log_pvs[q-1] + 
    sorted_weighted_neg_log_pvs[q]
}


for(i in 1:iterations){

exp_order = simulate_order(lambdas = exp_lambda)

A_vec = rep(NA,num_lambda)
for(q in 1:num_lambda){
  A_vec[q] = 1/sum(exp_lambda[exp_order[q:num_lambda]])
}

exp_rvs = rexp(num_lambda-1,rate = 1)

eos_niid = rep(0,num_lambda)
top_boss_niid = rep(NA,num_lambda)

eos_niid[1] = exp_rvs[1]*A_vec[1]
if(num_lambda>2){
for(q in 2:(num_lambda-1)){
 eos_niid[q] = eos_niid[q-1] + A_vec[q]*exp_rvs[q]
}
}

for(q in 1:num_lambda){
  top_boss_niid[q] = sum(eos_niid[num_lambda:(1+ num_lambda-q)])
}

t_ = sapply(boss_weighted_neg_log_pvs - rep(eos_niid[num_lambda-1],num_lambda) -top_boss_niid, function(x){max(x,0)})

top_boss_niid_pv_est = rep(NA,2)
for(q in 1:num_lambda){
  top_boss_niid_pv_est[q] = pexp(t_[q], rate = 1/A_vec[num_lambda], lower.tail = F)
}
pv_estimates = pv_estimates + top_boss_niid_pv_est
}

pv_estimates/iterations
}


# 
# 
# pv_trial = runif(2);pv_trial
# pv_trial = c(.001,.002)
# weights_trial = c(10,10)
# compute_weighted_boss_pv(pv_trial,weights_trial, iterations = 20000)
# weights_trial = c(2,1,9)
# compute_weighted_boss_pv(pv_trial,weights_trial, iterations = 20000)
# weights_trial = c(1,1,1)
# compute_weighted_boss_pv(pv_trial,weights_trial, iterations = 20000)
# weights_trial = c(1,1,1)/3
# compute_weighted_boss_pv(pv_trial,weights_trial, iterations = 20000)
# 
# 
# compute_min_boss_pv(pv_trial, parallelize = F, compute_ent = F, sim_iters = 20000)
# 
# compute_weighted_boss_pv(pv_trial,weights = c(1,1,1), iterations = 20000)


