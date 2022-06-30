##################################
#### Estimating Beta Function ####
##################################

estimate_beta_functions <- function(x, y, ySmooth, fbasis)
{
  
  t = dim(x)[1]
  n = dim(x)[2]
  d = dim(x)[3]
  
  # Basis and parameter objects
  
  y_fd = ySmooth$fd
  
  # We can retain xlist from the scalar response model. 
  
  # Now we need to set up a list of covariates.
  
  xlist = list(len=d+1)
  
  # First co-variate is just the intercept: a vector of ones
  
  xlist[[1]] = rep(1,n)
  
  
  # Other covariates
  for (j in 1:d)
  {
    xSmooth = smooth.basis(1:t, x[,,j], fbasis)
    x_fd =xSmooth$fd
    
    xlist[[j+1]] = x_fd		
  }
  
  
  
  #### 1. fdPar objects and estimation
  
  bwtlist2 = list(len=d+1)
  
  # The intercept is now a functional parameter as well as beta 1.   Since this
  # is an identifiable model without smoothing, we'll set the smoothing parameter 
  # very low. 
  
  harmLfd = vec2Lfd(c(0,(2*pi/(t))^2,0),rangeval=c(0,t))
  beta.fdPar2 = fdPar(fbasis,harmLfd,1e-5)
  
  for (j in 1:(d+1))
    bwtlist2[[j]] = beta.fdPar2
  
  
  
  # Regression fit
  
  fit = fRegress(y_fd,xlist,bwtlist2)
  
  return(fit)
}
