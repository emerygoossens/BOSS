#boss Functions
#Expectation, Variance, and Covariance 
#of the Sum of Exponential Order Statistics

eos.iid.var = function(index = 5,n_total = 6, lambda = 1){

  if( index <= n_total){
    a_j = 1/(lambda*(n_total - seq(1,index,1) + 1))
    var = sum((a_j)^2)
    var
  }else{
    print("Error: index must be less than or equal to n_total")
  }
  
}
eos.iid.cov = function(index_1 = 3, index_2 = 4, n_total_cov){
  if(index_1 < index_2){
    covar = eos.iid.var(index = index_1, n_total_cov) + eos.iid.var(index_2, n_total_cov)
    covar
  }else{
    print("Error: index_1 must be less than index_2")
  }
}

eos.iid.cov.mat = function(total_variables){
  
  covar.mat = matrix(NA, nrow = total_variables, ncol = total_variables)
  cor.mat = covar.mat
  
  for(j in 1:total_variables){
    temp = rep(eos.iid.var(j, n_total = total_variables),length(j:total_variables))
    covar.mat[j,j:total_variables] = temp
    covar.mat[j:total_variables,j] = temp  
  }
  
  for(j in 1:total_variables){
    for(k in 1:total_variables){
      temp = covar.mat[j,k]/sqrt(covar.mat[j,j]*covar.mat[k,k])
      cor.mat[j,k] = temp
      cor.mat[k,j] = temp
    }
  }
  
  list("covariance.matrix" = covar.mat, "correlation.matrix" = cor.mat)
}


boss.iid.cov  = function(index_1 = 2,index_2 = 5, eos_cov_mat){  
  
  n_total_boss = ncol(eos_cov_mat)
  
  comp_1 = 0
  for(u in index_1:(index_2-1)){
    comp_1 = comp_1 +  sum(eos_cov_mat[u,index_2:n_total_boss])
  }
  comp_2 = sum(diag(eos_cov_mat)[index_2:n_total_boss])
  boss_cov = comp_1 + comp_2
  boss_cov
}

boss.iid.cov.mat = function(total_num_cov = 5, parallelize = F){

  l = total_num_cov
  eos_iid_cov_mat = eos.iid.cov.mat(total_num_cov)$covariance.matrix
  eos_diag = diag(eos_iid_cov_mat)
  boss_iid_cov_mat = matrix(NA, ncol =  total_num_cov, nrow = total_num_cov)
  boss_iid_cor_mat = matrix(NA, ncol =  total_num_cov, nrow = total_num_cov)

  for(j in 1:total_num_cov){
    boss_iid_cov_mat[j,j] = sum(eos_diag[j:total_num_cov]) + 2*sum(eos_iid_cov_mat[j:l,j:l][upper.tri(eos_iid_cov_mat[j:l,j:l], diag = F)])
  }
  upper = eos_iid_cov_mat
  upper[lower.tri(upper)]  = 0
  diag(upper) = 0 
  if(parallelize == F){
  for(j in 1:total_num_cov){
    for(k in min((j+1),total_num_cov):total_num_cov){
        temp = boss_iid_cov_mat[k,k] + sum(upper[j:( k - 1),k:l])
        #boss.iid.cov(j,k,eos_iid_cov_mat)
        boss_iid_cov_mat[j,k] = temp
        boss_iid_cov_mat[k,j] = temp
        temp2 = boss_iid_cov_mat[j,k]/sqrt(boss_iid_cov_mat[j,j]*boss_iid_cov_mat[k,k]) 
        boss_iid_cor_mat[j,k] = temp2
        boss_iid_cor_mat[k,j] = temp2
      
    }
  }

  diag(boss_iid_cor_mat) = 1
  }else{

    registerDoParallel(cores = detectCores() -1)
    temp_cov_mat = foreach(j=1:total_num_cov, .combine = rbind) %dopar% {
    par_row = rep(NA,total_num_cov)
    for(k in min((j+1),total_num_cov):total_num_cov){
      
        par_row[k] = boss_iid_cov_mat[k,k] + sum(upper[j:( k - 1),k:l])

    }
    par_row
  }
  for(j in 1:total_num_cov){
    for(k in min((j+1),total_num_cov):total_num_cov){
        boss_iid_cov_mat[j,k] = temp_cov_mat[j,k]
        boss_iid_cov_mat[k,j] = temp_cov_mat[j,k]
        temp = boss_iid_cov_mat[j,k]/sqrt(boss_iid_cov_mat[j,j]*boss_iid_cov_mat[k,k]) 
        boss_iid_cor_mat[j,k] = temp
        boss_iid_cor_mat[k,j] = temp
      }}

  diag(boss_iid_cor_mat) = 1
  }
  list("covariance.matrix" = boss_iid_cov_mat, "correlation.matrix" = boss_iid_cor_mat)
}


geg1 = function(z){
  1 - exp(-z)*(z + 1)
}

geg2 = function(z,w){
  exp(-z)*(z-1) + exp(-w)*(w - 2*z - 1)
}

theoretical_ent = function(ent2,alpha = alpha){

  alpha_corrected = alpha/ent2
  #alpha_ = 0.1/2
  #The threshold for the sum of two exponentials would be

  sum_thresh = qgamma(p = 1-alpha_corrected,shape = 2, rate = 1, lower.tail = T)

  #The threshold for the max would be

  max_thresh = qgen.exp(p = 1-alpha_corrected, alpha = 2, lambda = 1, lower.tail = T)

  #We are then interested in the
  #P(min(X, Y) < alpha_corrected) = 1 - P(X > sum_thresh, Y > max_thresh)

  # 1- 3.7*exp(-3.7) - exp(-3.7)



  (1 - geg1(max_thresh) - geg2(max_thresh,sum_thresh)) - alpha

}

find_ent2 = function(alpha = 0.05){
  uniroot(theoretical_ent,c(1,2), alpha = alpha)$root
}



boss.effective.num.tests= function(num_tests = 5, alpha = 0.05){
  cor_mat = boss.iid.cov.mat(num_tests, parallelize = F)$correlation.matrix
  eigen_vals = eigen(cor_mat, symmetric = T, only.values = T)$values
  li = num_tests - sum((eigen_vals-1)*ifelse(eigen_vals>1,1,0))
  ji_li = sum(ifelse(abs(eigen_vals)>1,1,0) + abs(eigen_vals) - floor(abs(eigen_vals)))
  V_lambda = sum(((eigen_vals - 1)^2)/(num_tests-1)) 
  chev = 1 + (num_tests-1)*( 1- V_lambda/num_tests)
  log_test = log(num_tests + 1)
  
  find_ent2 = function(alpha = 0.05){
    uniroot(theoretical_ent,c(1,2), alpha = alpha)$root
  }
  cor_mat2 = boss.iid.cov.mat(2, parallelize = F)$correlation.matrix
  eigen_vals2 = eigen(cor_mat2, symmetric = T, only.values = T)$values
  li2 = 2 - sum((eigen_vals2-1)*ifelse(eigen_vals2>1,1,0))
  
  delta = max(find_ent2(alpha = alpha) - li2 ,0)
  boss_ent = delta  + li
  list(chev = chev, ji_li = ji_li , li = li,  log_test = log_test, boss_ent = boss_ent)
}

boss.effective.num.tests.corr.mat = function(corr_mat){
  M = nrow(corr_mat)
  eigen_vals = eigen(corr_mat, symmetric = T, only.values = T)$values
  li = length(diag(corr_mat)) - sum((eigen_vals-1)*ifelse(eigen_vals>1,1,0))
  ji_li = ifelse(abs(eigen_vals)>1,1,0) + abs(eigen_vals) - floor(abs(eigen_vals))
  V_lambda = sum(((eigen_vals - 1)^2)/(M-1)) 
  chev = 1 + (M-1)*( 1- V_lambda/M)
}



empirical.cov.mat  = function(total_exp = 5, n_obs = 100000){
  output = matrix(NA,nrow = n_obs, ncol = total_exp)
  for(k in 1:n_obs){
    unifs = sort(runif(total_exp), decreasing = T)
    eos = -log(unifs)
    boss = eos
    for(j in 1:length(eos)){
      boss[j] = sum(eos[j:length(eos)])
    }
    output[k,] = boss
  }
  list("covariance" = cov(output), "correlation"  = cor(output))
}


ses.iid.cov.mat = function(num_exp){
  num_exp = 10
  variance = seq(1,num_exp,1) 
  var_mat = diag(variance)
  for(i in 1:num_exp){
    for(j in (i):(num_exp)){
      var_mat[i,j] = i
    }
  }
  var_mat[lower.tri(var_mat)] = t(var_mat)[lower.tri(var_mat)]
  cor_mat = var_mat
  for(i  in 1:num_exp){
    for(j in 1:num_exp){
      cor_mat[i,j] = var_mat[i,j]/sqrt(var_mat[i,i]*var_mat[j,j])
    }
  }
  
  list(var_mat = var_mat, cor_mat = cor_mat)
}




