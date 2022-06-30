############################################
#### Estimate Confidence Bands Function ####
############################################

estimate_CI_bands <- function(fit, ySmooth, y)
{
  t = dim(y)[1]
  #n = dim(y)[2]
  
  y2cmap = ySmooth$y2cMap
  
  # This is the matrix that goes from the observations to the coefficients. 
  
  # We now need a covariance matrix for the errors of the original observed 
  # precipitation from the functional linear model
  
  Errmat = y - eval.fd(1:t, fit$yhatfdobj) 
  
  SigmaE2 = cov(t(Errmat))
  
  # We can now run fRegress.stderr
  
  fit.std = fRegress.stderr(fit, y2cmap, SigmaE2)
  
  
  betaestlist = fit$betaestlist
  betastderrlist = fit.std$betastderrlist
  rangeval = betaestlist[[1]]$fd$basis$rangeval
  argvals = seq(rangeval[1],rangeval[2],len=t)
  p = length(betaestlist)
  
  result = list()
  for (j in 1:p) 
  {
    if (is.fdPar(betaestlist[[j]]))
      betavec = eval.fd(argvals, betaestlist[[j]]$fd)
    else 
    {
      if (is.fd(betaestlist[[j]])) 
        betavec = eval.fd(argvals, betaestlist[[j]])
      else 
        stop("BETAESTLIST does not contain a functional parameter or data object.")
    }
    
    betastderr = eval.fd(argvals, betastderrlist[[j]])
    
    result[[j]] = list()
    result[[j]]$betavec = betavec
    result[[j]]$upper = betavec + 1.96*betastderr
    result[[j]]$lower = betavec - 1.96*betastderr
    result[[j]]$argvals = argvals
  }
  return(result)
}