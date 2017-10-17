
t = 2
root_n_delta =1

beta_match = -t/log(1-pnorm(qnorm(1-exp(-t))-root_n_delta))

pdf_2sided_alt = function(p, sqrtn_delta=root_n_delta){
  dnorm(qnorm(1-p/2) - sqrtn_delta/4)/dnorm(qnorm(1-p/2))
}

pdf_1sided_alt = function(p, sqrtn_delta=root_n_delta){
  dnorm(qnorm(1-p) - sqrtn_delta)/dnorm(qnorm(1-p))
}



  

pdf_y = function(y, sqrtn_delta=root_n_delta){
  zexpnegy = qnorm(1-exp(-y))
  (dnorm(zexpnegy - sqrtn_delta)/dnorm(zexpnegy))*exp(-y)
}

y_pdf_y = function(y){
  y*pdf_y(y = y)
}
beta_dist = function(x,beta = beta_match){
    exp(-x/beta)/beta
  }

curve(pdf_y,0,5)
curve(beta_dist,0,5, beta = nuexp_mean, add=T)
# beta = t/(-log(pnorm(qnorm(exp(-t))-root_n_delta))); beta
# exp(-(1)*t)
# 
# pnorm(qnorm(0))

# -t/log(pnorm(1-qnorm(1-exp(-t)) - root_n_delta))


neg_log_alt_pv_ev_integrand = function(v){
  -log(pnorm(-v))*dnorm(v - 1)
}

curve(neg_log_alt_pv_ev_integrand,-10,10)

integrate(neg_log_alt_pv_ev_integrand,-10,10)
integrate(y_pdf_y,0,20)


curve(pdf_2sided_alt, 0,1)
curve(pdf_1sided_alt, 0,1, add = T)

integrate(pdf_1sided_alt,0,1)
integrate(pdf_2sided_alt,0,1)

box_cox = function(x, lambda = 1/2){
  (x^lambda -1)/lambda
}

curve(box_cox, from = 0, to = 2)
