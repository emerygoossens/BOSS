#library(gpuR)

compute_is_proposal_dist = function(cutoff, covmat, tst_stat){
  mu_is = solve(covmat)%*%cutoff
  #mgf_is  = exp(.5*(t(mu_is)%*%covmat%*%mu_is))
  weight_is = exp(.5*(t(mu_is)%*%covmat%*%mu_is) - t(mu_is)%*%tst_stat )
  list(mu_is = mu_is, weight_is = weight_is)
}

compute_is_weight = function(test_stat_sim_, mu_zero, mu_prop, chol_sigma){
  exp(dmvn(X = test_stat_sim_, mu = mu_zero, sigma = chol_sigma, isChol = T, log = T) - dmvn(X = test_stat_sim_, mu = mu_prop, sigma = chol_sigma, isChol = T, log = T))
}
compute_is_weight2 = function(test_stat_sim_, mu_zero, mu_prop, chol_sigma){
  dmvn(X = test_stat_sim_, mu = mu_zero, sigma = chol_sigma, isChol = T, log = F)#/dmvn(X = test_stat_sim_, mu = mu_prop, sigma = chol_sigma, isChol = T, log = F)
}

compute_is_weight3 = function(test_stat_sim_, mu_prop, sigma_inv){
  exp(-.5*(mu_prop%*%sigma_inv%*%(t(2*test_stat_sim_ - mu_prop))))
}

compute_is_weight4 = function(test_stat_sim_, mu_prop, mu_times_sigma_inv){
  exp(mu_times_sigma_inv%*%(t(2*test_stat_sim_ - mu_prop)))
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

corr_pv_mc_rbmc = function(min_pv_test_stat, inner_iter = 100, max_iter = 10, A_mat, cores_ = 12){
  p = nrow(A_mat)
  registerDoParallel(cores = cores_)
  
    simulations = foreach(icount(max_iter), .combine = rbind) %dopar% {
    
    # simmed_stat = simulate_mvn_pvalues(p = p, sigma = emp_cov_mvn, n = 1)
    test_stat_sim = rmvn(n = 1, mu = rep(0,p),sigma = A_mat, isChol = T)
    pvs_sim = compute_marginal_pvalues(test_stat_sim)
      boss_output = compute_min_boss_pv_c(pvs_sim, sim_iters = inner_iter, compute_ent = F, parallelize = F)
      
  boss_output$min_pv

    }
    mc_min_boss_pv = length(which(simulations < min_pv_test_stat))/max_iter
  list(simulations = simulations, mc_min_boss_pv = mc_min_boss_pv)
}

corr_pv_is_rb_mc = function(min_pv_test_stat, test_stat, contrib_ind, emp_cov_mvn, simul_type = "RB", is_scale = 10, parallelize_is = F, inner_iter = 100, max_iter = 1000, A_mat, cores_ = 12){
  #p = 10; iters = 10; unused_cores = 1 ; emp_cov_mvn = Posdef(10)
  #test_stat = test_stat_apical ; emp_cov_mvn = null_test_cov ; unused_cores = 1; cores = 23 ; 
  # min_pv = seos_naive_apical_pv$min_pv; test_stat = test_stat_apical;
  # emp_cov_mvn = diag(diag(corrmat_apical)); simul_type = "RB";is_scale = 6
  # parallelize_is = F; inner_iter = 100; max_iter = 3; 
  p = nrow(emp_cov_mvn)
  sign_test_stat = sign(test_stat)

  is_mean = test_stat/is_scale
  zero_mean = rep(0,p)
  registerDoParallel(cores = cores_)

  if(parallelize_is == T){

    qq = max_iter
    simulations = foreach(icount(max_iter), .combine = rbind) %dopar% {

      test_stat_sim = rmvn(n = 1, mu = is_mean, sigma = A_mat, isChol = T)
      pvs_sim = compute_marginal_pvalues(test_stat_sim)
      
      boss_output = compute_min_boss_pv_c(pvs_sim, sim_iters = inner_iter, compute_ent = F, parallelize = F)
      min_pv_sim = boss_output$min_pv
      if(min_pv_sim <min_pv_test_stat){
      is_weight = exp(dmvn(X = test_stat_sim, mu = zero_mean, sigma = A_mat, isChol = T, 
       log = T) - dmvn(X = test_stat_sim, 
        mu = is_mean, sigma = A_mat, isChol = T, log = T))}
      else{
        is_weight = 0
      }
      # sign_test_stat_sim = sign(test_stat_sim)
      
      # pvs_sim = ifelse(sign_test_stat == 1, 1 - pnorm(test_stat_sim, 
        # sd = diag(emp_cov_mvn)), pnorm(test_stat_sim, sd = diag(emp_cov_mvn)))
      

      

      c(min_pv_sim,is_weight)
    }
  }else{
    
    simulations  = matrix(NA, nrow = max_iter, ncol = 2)
    rel_error = Inf
    current_estimate = 10
    small_change = 0
    qq = 1
    accurate = F

    while(!accurate && qq < max_iter){

    
      test_stat_sim = rmvn(n = 1, mu = is_mean, sigma = A_mat, isChol = T)

      is_weight = exp(dmvn(X = test_stat_sim, mu = rep(0,p), 
        sigma = emp_cov_mvn, log = T) - dmvn(X = test_stat_sim, 
        mu = is_mean, sigma = emp_cov_mvn, log = T))
      
      sign_test_stat_sim = sign(test_stat_sim)
      pvs_sim = ifelse(sign_test_stat == 1, 
        1 - pnorm(test_stat_sim, sd = diag(emp_cov_mvn)), 
        pnorm(test_stat_sim, sd = diag(emp_cov_mvn)))
      seos_output = compute_min_SEOS_pv(pvs_sim, sim_iters = inner_iter, compute_ent = F, parallelize = F)
      min_pv_sim = seos_output$min_pv
      #min_pv_ind = min(which(ordered_combo_sim == min_pv_sim))
      # fishers_all_combo_pvs = fishers_for_all_ordered_subsets(pvs_sim)
      # min_pv_sim = min(fishers_all_combo_pvs)
      # min_pv_ind = which(fishers_all_combo_pvs == min_pv_sim)

      # fishers_all_combo_pvs = fishers_for_all_ordered_subsets(pvs_sim)
      # min_pv_sim = min(fishers_all_combo_pvs)
      # min_pv_ind = which(fishers_all_combo_pvs == min_pv_sim)

      
      # if(pool_neighbors){
      #   neighbor_ind = c(min_pv_ind -1:3, min_pv_ind, min_pv_ind + 1:3)
      #   neighbor_ind = neighbor_ind[ neighbor_ind <= p]
      #   neighbor_ind = neighbor_ind[ neighbor_ind >= 1]
      #   min_pv_sim = mean(ordered_combo_sim[neighbor_ind])
      # }

      simulations[qq,] = c(min_pv_sim,is_weight)
      # if(qq > min_iter && qq %% step == 0){
      #   past_estimate = current_estimate
      #   weights_to_sum = sum(simulations[which(simulations[1:qq,1]<min_pv_test_stat),2])
      #   current_estimate  = weights_to_sum/qq
      #   delta_rel_error = abs(rel_error - relative_error(current_estimate, past_estimate))
      #   rel_error = relative_error(current_estimate, past_estimate)
      #   #print(paste("current estimate:",current_estimate))
      #   #print(paste("relative error:",rel_error))
      #   #print(paste("ordered_combo sample:"))
      #   #print(ordered_combo_sim)
      #   if(rel_error < rel_err_cutoff && current_estimate >0){
      #     small_change = small_change + 1
      #   }
      #   if(small_change>1){
      #     accurate = T
      #   }
      #   if(current_estimate > .2){
      #     accurate = T
      #   }
        
      # }
      qq  = qq + 1
    }
    if(abs(qq - max_iter) <3){
      qq = max_iter
    }else{
      qq = qq - 1
    }
    }
  rownames(simulations) = NULL
  colnames(simulations) = NULL
  min_is_rbmc_pv = min(sum(simulations[which(simulations[,1]<min_pv_test_stat),2])/qq,.999)
  list(min_is_rbmc_pv = min_is_rbmc_pv, simulations = simulations)
}

relative_error = function(estimate, value){
    abs(1 - estimate/value)
}


relative_error(.001, .0001)




