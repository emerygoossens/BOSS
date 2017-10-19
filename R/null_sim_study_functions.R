

sim_null = function(n = 1000, p = 10, sigma= diag(rep(1,p)),iters = 10, boss_iters = 10000, unused_cores = 2){
  num_cores =detectCores()- unused_cores
  print(paste("Starting Null Simulation for ",p," Parameters, ", iters, "Iterations Using ", num_cores, " Cores"), sep = "")
  effective_num_tests = boss.effective.num.tests(p)$boss_ent
  registerDoParallel(cores=num_cores)
  
  simulations = foreach(icount(iters), .combine = rbind) %dopar% {
    output= rep(NA,5)
    sim_setup = simulate_lm_pvalues(p = p, sigma = sigma)
    univariate_pvals = sim_setup$univariate_pvals
    boss_output = pvalues_to_boss(univariate_pvals, boss_iters = boss_iters)
    output[1] = boss_output$true_pv_fishers
    options(scipen = 99)
    fsum = summary(lm(sim_setup$y~0 +sim_setup$x))
    fstat = fsum$fstatistic
    f_pval = 1 -pf(fstat[1],fstat[2],fstat[3])
    output[2] = f_pval
    output[3] = boss_output$ord_comb_pvs[1]
    output[4] = min(min(boss_output$ord_comb_pvs)*effective_num_tests,.9999)
    output[5] = boss_output$num_selected
    output
  }
  
  rownames(simulations) = NULL
  colnames(simulations) = c("Fishers","F-test","Min-PV","BOSS", "#BOSS")
  simulations
}
  
sim_null_par = function(n = 1000, p = 10, sigma= diag(rep(1,p)),iters = 10, boss_iters = 10000, unused_cores = 1){
  
  #n = 100; p = 30; sigma = 1;iters = 100; iter = 1; output = matrix(NA, nrow = iters, ncol = 4) ; true_cor_mat = boss.iid.cov.mat(p)$correlation.matrix
  #output = matrix(NA, nrow = iters, ncol = 4)
  #true_cor_mat = boss.iid.cov.mat(p)$correlation.matrix
  num_cores =detectCores()- unused_cores
  print(paste("Starting Null Simulation for ",p," Parameters, ", iters, "Iterations Using ", num_cores, " Cores"), sep = "")
  
  effective_num_tests = boss.effective.num.tests(p)$boss_ent
  
  
  registerDoParallel(cores=num_cores)
  simulations = foreach(icount(iters), .combine = rbind) %dopar% {
    output= rep(NA,5)
    
    # x=matrix(rnorm(n*p),n,p)
    # x=scale(x,T,F)
    # beta= rep(0,p)
    # y=x%*%beta+sigma*rnorm(n)
    # y=y-mean(y)
    # 
    # univariate_pvals = matrix(NA, nrow = 1, ncol = p)
    # 
    # for(j in 1:p){
    #   univariate_pvals[,j] = coef(summary(lm(y ~ 0 + x[,j])))[1,4]
    #   
    # }
    sim_setup = simulate_lm_pvalues(p = p, sigma = sigma)
    univariate_pvals = sim_setup$univariate_pvals
    boss_output = pvalues_to_boss(univariate_pvals, boss_iters = boss_iters)
    # combo_pval_all = fishers_method(univariate_pvals)
    output[1] = boss_output$true_pv_fishers
    # ord_stat_pvals = sort(univariate_pvals)
    # ord_stat_exp = -log(ord_stat_pvals)
    # sum_ord_stat_exp =  rep(NA,p)
    # sum_ord_stat_exp[1] = ord_stat_exp[1]
    # for(j in 2:p){
    #   sum_ord_stat_exp[j] = ord_stat_exp[j] + sum_ord_stat_exp[j-1]
    # }
    # 
    # 
    # p_val_ord_stat_exp_sim = p_sum_ord_exp_sim_simul(sum_ord_stat_exp, iters = boss_iters)$ord_comb_pvs
    options(scipen = 99)
    fsum = summary(lm(sim_setup$y~0 +sim_setup$x))
    fstat = fsum$fstatistic
    f_pval = 1 -pf(fstat[1],fstat[2],fstat[3])
    output[2]= f_pval
    output[3] = boss_output$ord_comb_pvs[1]
    output[4] = min(min(boss_output$ord_comb_pvs)*effective_num_tests,.9999)#ifelse(min(p_val_ord_stat_exp_sim)<.1,ent_corrected,output[3])##correct_pv(p_val_ord_stat_exp_sim, fwer = 0.05)
    
    num_selected = which(min(boss_output$ord_comb_pvs) == boss_output$ord_comb_pvs)

    output[5] = num_selected[1]
    output
  }
  rownames(simulations) = NULL
  colnames(simulations) = c("Fishers","F-test","Min-PV","BOSS", "#BOSS")
  
  list(p_val_diagnostics = simulations)
}



simulation_study_null = function(iterations = 100, num_params = c(10), unused_cores = 1){
  start = Sys.time()
  
  for(l in 1:length(num_params)){
    
    sim_out = sim_null_par(p = num_params[l], iters = iterations, unused_cores = unused_cores)
    
    save(sim_out, file = paste("../Output/null_sim_study/null_sim_out_",num_params[l],"_iters",iterations,"_",Sys.Date(),".Rda",sep = ""))
  }
}


make_null_sim_plots = function(output_file){

    sim_out = load(output_file)
    output = sim_out$p_val_diagnostics
    num_params = max(output[,5])
    names = colnames(output)
    
    #c("Fishers","F-test","Min-PV","BOSS", "#BOSS")
    print(paste("Making Plots for Simulation with ",num_params,"Parameters"))
    #Plot 1
    png(paste("plots/null_simulation_study/pv_comp_p_",num_params,"_iters",iterations,"_",Sys.Date(),".png",sep = ""),width = 1400, height = 1050)
    par(mfrow = c(1,3),oma = c(0, 0, 2, 0), mar=rep(7,4))
    ramp <- colorRamp(c("red", "black"))
    plot_colors = rgb( ramp(seq(0, 1, length = num_params)), max = 255)
    
    for(j in 1:2){
      indices = which(output[,4]==1)
      plot(output[indices,3], output[indices,j], xlab = paste(names[3]," PV",sep = ""),ylab =paste(names[j]," PV",sep = ""), col = plot_colors[1],pch = 19, xlim = c(0,1), ylim = c(0,1), cex = 3, cex.axis = 3, cex.lab = 3)
      if(j ==1){
        legend(x = .2,y =.1, c("1 Selected BOSS",paste(num_params,"Selected BOSS")), cex = 3, pch = c(19,19),col=c('red','black'), box.col   = "white")
      }
      abline(c(0,0),c(1,1))
      for(k in 2:max_params){
        indices = which(output[,4]==k)
        points(output[indices,3],output[indices,j], col = plot_colors[k],pch = 19,cex = 3)
      }
    }
    
    indices = which(output[,4]==1)
    
    plot(output[indices,3], output[indices,5],cex = 3, cex.axis = 3, xlab = paste(names[3]," PV",sep = ""),ylab =paste(names[5]," PV",sep = ""), col = plot_colors[1], pch = 19, xlim = c(0,1), ylim = c(0,1), cex.lab = 3)
    abline(c(0,0),c(1,1))
    
    for(k in 2:max_params){
      indices = which(output[,4]==k)
      points(output[indices,3],output[indices,5], cex = 3, cex.axis = 3, col = plot_colors[k],pch = 19)
    }
    mtext(paste("Comparison of P-Values for BOSS, Fisher's, and Min. p = ",num_params, ", Iterations = ",iterations, sep = ""), outer = TRUE, cex = 3, line = -2)
    dev.off()
    
    
    #Plot 2
    png(paste("plots/null_simulation_study/qq_hist_",num_params,"_iters",iterations,"_",,Sys.Date(),".png",sep = ""),width = 1400, height = 800)
    par(mfrow = c(2,3),oma = c(0, 0, 2, 0), mar=rep(8,4))
    
    max_log_p = ceiling(max(-log(output[,1:3], base = 10)))
    
    qqunif(output[,1], log = T, main="Fisher", xlim = c(0,max_log_p), ylim = c(0,max_log_p),cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    qqunif(output[,2], log = T, main="Min", xlim = c(0,max_log_p), ylim = c(0,max_log_p),cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    qqunif(output[,3], log = T, main="BOSS", xlim = c(0,max_log_p), ylim = c(0,max_log_p),cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    
    hist(output[,1], breaks = 25, xlim = c(0,1), main = "Fisher", xlab="P-Value", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    hist(output[,2], breaks = 25, xlim = c(0,1), main ="Min" , xlab="P-Value", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    hist(output[,3], breaks = 25, xlim = c(0,1), main = "BOSS", xlab="P-Value", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3)
    
    mtext(paste("Comparison of P-Values for BOSS, Fisher's, and F-test p = ",num_params, sep = ""), outer = TRUE, cex = 3, line = -2)
    dev.off()
    
    #Plot 3
    png(paste("plots/null_simulation_study/low_quantile_p_",num_params,"_iters",iterations,"_",,Sys.Date(),".png",sep = ""),width = 1400, height = 1050)
    par(mfrow = c(1,1),oma = c(0, 0, 2, 0), mar=rep(8,4))
    
    pquantile = seq(0,.1,.00001)
    pobserved = pquantile
    for(j in 1:length(pquantile)){
      pobserved[j] = length(which(output[,3]<pquantile[j]))/iterations
    }
    
    plot(pquantile,pobserved, xlab = "Expected", ylab = "Observed", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3, xlim = c(0,.1), ylim = c(0,.1) )
    abline(c(0,0),c(1,1))
    mtext(paste("Observed vs. Expected Quantiles for Lower Tail p = ",num_params, sep = ""), outer = TRUE, cex = 3, line = -2)
    dev.off()
  }
  