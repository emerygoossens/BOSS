library(gap)
setwd("~/Documents/BOSS/R")


source("exp_covariance.R")
#load("/Users/egoossen/Downloads/null_sim_innter_30k_30_iters1e+05_2017-01-07.Rda")
load("/Users/egoossen/Documents/BOSS/Data/ENT/null_sim_out_30_iters1e+05_2017-01-06.Rda")
sim_out30 = sim_out$p_val_diagnostics
load("/Users/egoossen/Documents/BOSS/Data/ENT/null_sim_out_20_iters1e+05_2017-01-06.Rda")
sim_out20 = sim_out$p_val_diagnostics
load("/Users/egoossen/Documents/BOSS/Data/ENT/null_sim_out_10_iters1e+05_2017-01-06.Rda")

sim_out10 = sim_out$p_val_diagnostics
load("/Users/egoossen/Documents/BOSS/Data/ENT/null_sim_inner_10k_50_iters1e+05_2017-01-08.Rda")
sim_out50 = sim_out$p_val_diagnostics
library(doParallel)
  


load("/Users/egoossen/Documents/BOSS/Data/ENT/boss_matrix.rda")

# qqunif(boss_matrix[,1])
indices = c(2,5,10,15,20,25,30,35,40,45,50)
max_ind = 50
multiplier = rep(NA, length(indices))
for(j in 1:ncol(boss_matrix)){
  sequence = seq(1,6,.01)
  pvs = boss_matrix[,j]
  for(s in sequence){
    
    type1_err = length(which(pvs*s<.05))/nrow(boss_matrix)
    if(type1_err <= 0.05){
      multiplier[j] = s
      break
    }
  }
}


ent_vs_ind = matrix(NA,nrow = length(indices), ncol = 6)
for(j in 1:length(indices)){
  ent_vs_ind[j,] = c(unlist(boss.effective.num.tests(indices[j])),multiplier[j])
}

ent_vs_ind[,3] + ent_vs_ind[1,5] - ent_vs_ind[1,3]

png(paste("../Plots/ent_plot",Sys.Date(),".png",sep = ""),width = 900, height = 900)
par(mfrow = c(1,1),oma = c(2, 2, 2, 2), mar = c(5,6,4,3))
   

plot(indices,ent_vs_ind[,1], type = "l", lwd = 3,xlim = c(1,50),ylim = c(0,9),
     col = 2, main = "Log Relationship and Effective Number Tests", 
     xlab = "Number of Tests", ylab = "Estimated Effective Number of Tests", 
     cex.lab = 2,cex = 2, cex.axis = 2, cex.main = 3)

lines(indices,ent_vs_ind[,2], col = 3, lwd = 3)
lines(indices,ent_vs_ind[,3], col = 4, lwd = 3)
lines(indices,ent_vs_ind[,4], col = 5, lwd = 3)
# curve( log(x), from = 2, to = max_ind, n = max_ind, add = T, col = 5, lwd = 3)
lines(indices,multiplier, col = 8, lwd = 3, lty = 3)

lines(indices,ent_vs_ind[,5] , col = 1,lty = 1, lwd = 3)

#points(indices,ent_vs_ind[,4])
#title(main = "Log Relationship and Effective Number Tests", cex = 4)
legend(x = 5, y = 8.5, legend = c("Cheverud (2001)", "Li & Ji (2005)", "Li (2012)", "Log(Tests +1)", "Empirical", "Proposed"), 
       col = c(2,3,4,5,8,1), lty = c(1,1,1,1,2,1), cex = 2, lwd = c(3,3,3,3,3,3))
dev.off()
# 
# sim_out$p_val_diagnostics[,4]
# 
# # length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(30)*log(30) < .1))/100000
# # length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(30)*log(30)<.01))/100000
# # length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(30)*log(30)<.001))/100000
# # length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(30)*log(30)<.0001))/100000
# # 
# 
# output = sim_out$p_val_diagnostics
# iterations = nrow(output)
# num_params = 50
# png(paste("low_quantile_ent_p_",num_params,"_iters",iterations,"_",Sys.Date(),".png",sep = ""),width = 1400, height = 1050)
# par(mfrow = c(1,1),oma = c(0, 0, 2, 0), mar=rep(8,4))
#     
# pquantile = seq(0,.1,.00001)
# pobserved = pquantile
# for(j in 1:length(pquantile)){
# pobserved[j] = length(which(output[,4]<pquantile[j]))/iterations
# }
# 
# plot(pquantile,pobserved, xlab = "Expected", ylab = "Observed", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3, xlim = c(0,.1), ylim = c(0,.1) )
# abline(c(0,0),c(1,1))
# mtext(paste("Effective # Tests: Observed vs. Expected Quantiles for Lower Tail p = ",num_params, sep = ""), outer = TRUE, cex = 3, line = -2)
# dev.off()
# 
# 
# new_output = output
# new_output[,4] = new_output[,4]*log(num_params)/boss.effective.num.tests(num_params)
# 
# png(paste("low_quantile_log_p_",num_params,"_iters",iterations,"_",Sys.Date(),".png",sep = ""),width = 1400, height = 1050)
# par(mfrow = c(1,1),oma = c(0, 0, 2, 0), mar=rep(8,4))
#     
# pquantile = seq(0,.1,.00001)
# pobserved = pquantile
# for(j in 1:length(pquantile)){
# pobserved[j] = length(which(new_output[,4]<pquantile[j]))/iterations
# }
# 
# plot(pquantile,pobserved, xlab = "Expected", ylab = "Observed", cex.lab = 3,cex = 3, cex.axis = 3, cex.main = 3, xlim = c(0,.1), ylim = c(0,.1) )
# abline(c(0,0),c(1,1))
# mtext(paste("Log Correction: Observed vs. Expected Quantiles for Lower Tail p = ",num_params, sep = ""), outer = TRUE, cex = 3, line = -2)
# dev.off()
# 
# 
# 
# length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(num_params)*log(num_params)<.1))/100000
# length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(num_params)*log(num_params)<.01))/100000
# length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(num_params)*log(num_params)<.001))/100000
# length(which(sim_out$p_val_diagnostics[,4]/boss.effective.num.tests(num_params)*log(num_params)<.0001))/100000
# 
# hist(sim_out$p_val_diagnostics[,5], breaks = 100)
# 
# first_moment = mean(sim_out$p_val_diagnostics[,5])
# sec_moment = var(sim_out$p_val_diagnostics[,5]) + first_moment^2
# 
# bb_mle = function(count_data, N ){
#   m = mean(count_data)
#   s2= var(count_data)
#   pie = m/N
#   theta = (s2 - N*pie*(1-pie))/((N^2)*pie*(1-pie) - s2)
#   list(pie = pie, theta = theta)
# }
# 
# bb_mle(sim_out10[,5], 10)
# bb_mle(sim_out20[,5], 20)
# bb_mle(sim_out30[,5], 30)
# 
# mean(sim_out10[,5])
# 
# x = sim_out30[,5]
# N = 30
# LL = function(alpha, beta){
#   R = dbb(x, u = alpha, v = beta, N = N)
#   -sum(log(R))
# }
# 
# mle_est30 = mle(LL, start = list(alpha = .5, beta = .5))
# est30_coef = coef(mle_est30)
# est30_coef[1]/(sum(est30_coef))
# 
# 
# x = sim_out20[,5]
# N = 20
# LL = function(alpha, beta){
#   R = dbb(x, u = alpha, v = beta, N = N)
#   -sum(log(R))
# }
# mle_est20 = mle(LL, start = list(alpha = .5, beta = .5))
# est20_coef = coef(mle_est20)
# est20_coef[1]/(sum(est20_coef))
# 
# 
# x = sim_out10[,5]
# N = 10
# LL = function(alpha, beta){
#   R = dbb(x, u = alpha, v = beta, N = N)
#   -sum(log(R))
# }
# mle_est10 = mle(LL, start = list(alpha = .4, beta = .6))
# est10_coef = coef(mle_est10)
# est10_coef[1]/(sum(est10_coef))
# 
# 
# 
# 
# 
# 
# 
# x = sim_out50[,5]
# N = 50
# LL = function(alpha, beta){
#   R = dbb(x, u = alpha, v = beta, N = N)
#   -sum(log(R))
# }
# mle_est50 = mle(LL, start = list(alpha = .4, beta = .6))
# est50_coef = coef(mle_est50)
# est50_coef[1]/(sum(est50_coef))*50
# 
# qqunif(sim_out50[,4]/boss.effective.num.tests(50)*log(50), log = T)
# 
# 
# 
# mean(sim_out10[,5])/10
# mean(sim_out20[,5])/20
# mean(sim_out30[,5])/30
# mean(sim_out50[,5])/50
# 
# 
