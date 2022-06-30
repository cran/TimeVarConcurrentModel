##################################
#### Hypothesis Test Function ####
##################################
HypothesisTestFunctionalData <- function(x,y, B = 500, n.basis,show_iter=FALSE)UseMethod("HypothesisTestFunctionalData")

HypothesisTestFunctionalData.default <- function(x,y, B = 500, n.basis,show_iter=FALSE)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  t = dim(x)[1] ## observed over t days
  n = dim(x)[2] ## number of "countries"
  d = dim(x)[3] ## number of variables
  
  
  fbasis =  create.fourier.basis(c(0,t), n.basis) #t-1)
  ySmooth = smooth.basis(1:t, y, fbasis)
  
  
  fit = estimate_beta_functions(x,y, ySmooth, fbasis)
  CI_bands = estimate_CI_bands(fit, ySmooth, y)
  
  p = length(fit$betaestlist)
  
  p.values = rep(1,p)
  for (j in 1:p) 
  {
    integral_H0 = get_bootstrap_null_distr(x,y, index_to_test = j, B = B, n.basis= n.basis,show_iter=show_iter)
    
    integral = sintegral(CI_bands[[j]]$argvals,pmax(CI_bands[[j]]$lower,0))$int
    integral = integral - sintegral(CI_bands[[j]]$argvals,pmin(CI_bands[[j]]$upper,0))$int  
    
    p.values[j] = length(integral_H0[which(integral_H0 > integral)])/B
  }
  
  
  return(p.values)
}