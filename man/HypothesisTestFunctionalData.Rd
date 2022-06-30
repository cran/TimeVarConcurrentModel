\name{HypothesisTestFunctionalData}
\alias{HypothesisTestFunctionalData}

\title{
Hypothesis Test for Concurrent Multivariate Regression Models With Time-Varying
Coefficients
}
\description{
 This function provides a hypothesis test which is used for concurrent multivariate
regression models that have time varying coefficients. The hypothesis test can be
used to determine if covariates were significant for any time "t" within a given
range of time.
}
\usage{
HypothesisTestFunctionalData(x,y, B = 500, n.basis, show_iter=FALSE)
}

\arguments{
  \item{x}{
  Array containing the data from the different covariates
  }
  \item{y}{
  Matrix containing the data from the response variable
  }
  \item{B}{
  Number of bootstraps used to determine the distribution for the test statistic
  }
  \item{n.basis}{
  Number of bases used for the bootstrapping process
  }
  \item{show_iter}{
  Let's the user choose whether or not they would like to see the progress of the bootstrap
  }
  
}
\details{
With the given parameters, this function will conduct the specified number of bootstraps and then use the results to estimate the beta functions and their associated confidence bands. With this information, the function will then use an integral to calculate the area where the upper confidence bands fall below 0 and where the lower confidence bands rise above 0. The total area calculated will be compared to the expected area based on the bootstraps and a p-value will then be determined for each covariate.

}
\value{
HypothesisTestFunctionalData will return a one-dimensional array with the p-values for each covariate
}

\author{
Adriano Zambom
}

\examples{
library(fda)
library(Bolstad2)

country = c(rep("Brazil",365), rep("US",365), rep("Canada",365), rep("UK",365))
cases = c(1:365, 2*(1:365), 9*(1:365), 7*(1:365))
stringency  = c(seq(10,20,length.out = 365), 10 + sin(1:365), 10 + cos(1:365), abs(tan(1:365)))
stringency2  = c(seq(10,20,length.out = 365)+rnorm(365), 
10 + sin(1:365)+abs(rnorm(365)), 10 + cos(1:365)+abs(rnorm(365)), abs(tan(1:365))+ abs(rnorm(365)))
data = data.frame(country, cases, stringency, stringency2)

n = dim(data)[2]
t = dim(data)[1]/dim(data)[2]

y = matrix(0,t,n)  
for (i in 1:n)
	y[,i] = data$cases[((i-1)*t + 1):(i*t)]


d = 2  ## set this manually for the covariates you want
n.basis = 8 ## sets number of bases to be used

x = array(0,c(t,n,d))  
for (i in 1:n)
{
	x[,i,1] = data$stringency[((i-1)*t + 1):(i*t)]
	x[,i,2] = data$stringency2[((i-1)*t + 1):(i*t)]
}

p.values = HypothesisTestFunctionalData(x,y, B = 5,n.basis)
p.values
}

