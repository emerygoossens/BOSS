

# library(gap)
setwd("~/Documents/BOSS/R")
source("effective_num_test_plot.R")



alpha_list = c(.1,.05,.01,.001,.0001)

num_par_list = c(10,20,30,50)

data_list = list(sim_out10, sim_out20, sim_out30, sim_out50)

num_iters = nrow(data_list[[1]])

test_list = c("Fishers", "F-test", "Min-PV","BOSS")
j = 1
k = 1
a = 1
empirical_typeIerror = array(NA, dim = c( length(test_list), length(alpha_list),length(num_par_list)))

for(t in 1:length(test_list)){
for(a in 1:length(alpha_list)){
for(p in 1:length(num_par_list)){


empirical_typeIerror[t,a,p] = length(which(data_list[[p]][,t]<alpha_list[a]))/num_iters

}}}

