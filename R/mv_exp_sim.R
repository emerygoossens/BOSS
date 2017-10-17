





nc_exp_density = function(e,nobs = 20, delta = .1){
  num = exp(-e)*dnorm(qnorm(1-exp(-e)) -sqrt(nobs)*delta)
  denom = dnorm(qnorm(1-exp(-e)))
  num/denom
}


exp(-x)*dnorm(qnorm(1-exp(-x)) -sqrt(20)*1)/dnorm(qnorm(1-exp(-x)))


curve(exp(-x)*dnorm(qnorm(1-exp(-x)) -sqrt(100)*.1)/dnorm(qnorm(1-exp(-x))),0,30)
curve(exp(-x)*dnorm(qnorm(1-exp(-x)) -sqrt(25)*1)/dnorm(qnorm(1-exp(-x))),0,30, add = T)
curve(exp(-x/4)/4, add = T)
