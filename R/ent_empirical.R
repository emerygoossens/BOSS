load("/Users/egoossen/Documents/phd_research/SEOS/BOSS/boss_matrix.rda")
library(gap)
qqunif(boss_matrix[,1])
indices = c(2,5,10,15,20,25,30,35,40,45,50)
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


mult_mat = cbind(multiplier,indices)
rownames(mult_mat) = indices
mult_mat[10]
