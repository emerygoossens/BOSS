#mvn tail importance sampling

test_stat = 
  
  
  
  maha_is = function(min_pv_test_stat, test_stat,  simul_type = "RB", is_scale = 10, parallelize_is = F, inner_iter = 100, max_iter = 1000, A_mat, cores_ = 12){
    #p = 10; iters = 10; unused_cores = 1 ; emp_cov_mvn = Posdef(10)
    #test_stat = test_stat_apical ; emp_cov_mvn = null_test_cov ; unused_cores = 1; cores = 23 ;
    # min_pv = seos_naive_apical_pv$min_pv; test_stat = test_stat_apical;
    # emp_cov_mvn = diag(diag(corrmat_apical)); simul_type = "RB";is_scale = 6
    # parallelize_is = F; inner_iter = 2; max_iter = 3; A_mat = chol_apical
    # is_scale = is_scale_
    
    p = length(test_stat)

    
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