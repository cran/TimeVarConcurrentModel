FitBandsPlot <- function(x,y,t,n.basis,covnames)UseMethod("FitBandsPlot")

FitBandsPlot.default <- function(x,y,t,n.basis,covnames)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  CovariateNames = covnames
  fbasis =  create.fourier.basis(c(0,t), n.basis)
  ySmooth = smooth.basis(1:t, y, fbasis)
  
  fit = estimate_beta_functions(x, y, ySmooth, fbasis)
  CI_bands = estimate_CI_bands(fit, ySmooth, y)
  
  y2cmap = ySmooth$y2cMap
  Errmat = y - eval.fd(1:t, fit$yhatfdobj) 		
  SigmaE2 = cov(t(Errmat))		
  fit.std = fRegress.stderr(fit, y2cmap, SigmaE2)
  
  myfitbetaplot(fit$betaestlist,fit.std$betastderrlist,CovariateNames)	
}
