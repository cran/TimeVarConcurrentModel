############################
#### Bootstrap Function ####
############################

### get the "bootstrapped" distribution of the test statistic by inducing the null hypothesis randomly
## The test statistic is the integral of the confidence bands that do not include 0
get_bootstrap_null_distr <- function(x,y, index_to_test = 2, B = 500, n.basis,show_iter=FALSE)  ##index_to_test = 1 would correspond to intercept
{
  t = dim(x)[1]
  n = dim(x)[2]
  
  if (show_iter==TRUE){
    if (index_to_test==1){
      message("\n Bootstrap sample for testing intercept: ")
    }
    else{
      message("\n Bootstrap sample for testing covariate ",index_to_test-1, ":" )
    }
  }
  
  integral = 0
  for (b in 1:B)
  {
    ### induce null hypothesis	
    y_aux = y
    ####### get beta0 ##################################################
    fbasis =  create.fourier.basis(c(0,t), n.basis) #T-1)
    ySmooth = smooth.basis(1:t, y_aux, fbasis)			
    fit = estimate_beta_functions(x, y_aux, ySmooth, fbasis)
    CI_bands = estimate_CI_bands(fit, ySmooth, y_aux)
    betaestlist = fit$betaestlist
    p = length(fit$betaestlist)
    betavec = list()
    for (j in 1:p)
    {
      if (is.fdPar(betaestlist[[j]]))
        betavec[[j]] = eval.fd(CI_bands[[j]]$argvals, betaestlist[[j]]$fd)
      else 
      {
        if (is.fd(betaestlist[[j]])) 
          betavec[[j]] = eval.fd(CI_bands[[j]]$argvals, betaestlist[[j]])
        else 
          stop("BETAESTLIST does not contain a functional parameter or data object.")
      }
    }
    ####################################################################################################
    for (i in 1:n)
    {
      for (j in 1:p)
        if (j != index_to_test)
          y_aux[,i] = y[,i] - betavec[[j]]  ## removing beta0(t) and all the other X(t)beta(t) except the one being tested
    }
    y_aux = matrix(sample(y_aux), t, n)
    
    fbasis =  create.fourier.basis(c(0,t), n.basis) #T-1)
    ySmooth = smooth.basis(1:t, y_aux, fbasis)
    
    fit = estimate_beta_functions(x, y_aux, ySmooth, fbasis)
    CI_bands = estimate_CI_bands(fit, ySmooth, y_aux)
    
    integral[b] = sintegral(CI_bands[[index_to_test]]$argvals,pmax(CI_bands[[index_to_test]]$lower,0))$int
    integral[b] = integral[b] - sintegral(CI_bands[[index_to_test]]$argvals,pmin(CI_bands[[index_to_test]]$upper,0))$int  
    ## minus here because I have the bands below 0 and I want the positive integral of it. Could have been abs
    
    if (show_iter==TRUE){
      replaceMessage(b)
    }
    
  }
  
  return(integral)
}