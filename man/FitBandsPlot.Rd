\name{FitBandsPlot}
\alias{FitBandsPlot}

\title{
Fitting Beta Functions with Confidence Bands
}
\description{
This function plots the beta function for each covariate along with its corresponding confidence bands. The area being summed for the calculation of the p-values are also highlighted within the plots so that the viewer can see the points in time where 0 is not within the confidence bands.
}
\usage{
FitBandsPlot(x,y,t,n.basis,covnames)
}

\arguments{
   \item{x}{
  Array containing the data from the different covariates.
  }
  \item{y}{
  Matrix containing the data from the response variable
  }
  \item{t}{
  Variable used for time. For example, if the data is observed daily, then T will denote the total oberved number of days.
  }
  \item{n.basis}{
  Number of bases used for the bootstrapping process
  }
  \item{covnames}{
  Column vector containing the names of all the covariates in x array. These names should be in the same order as they appear in the x array.
  }
}
\details{
This function will plot the beta functions with their corresponding confidence bands. Each plot will be labeled accordingly based on the argument "covnames". The beta values for each covariate will be represented on the y-axis, while the variable of time will be represented on the x-axis. Each plot will shade the area that is being computed by the integrals in the the hypothesis test function. The intercept function will not include this shading.
}
\value{
FitBandsPlot will output n+1 plots, where n is the number of covariates in the x array.
}
\references{
fitbetaplot function originally from fda package. Written by J. O. Ramsay, , Hadley Wickham, Spencer Graves, Giles Hooker
Edited by Chance Vandergeugten and renamed myfitbetaplot
}
\author{
Adriano Zambom
Chance Vandergeugten
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
n.basis = 10 ## sets number of bases to be used
covnames = c("Stringency","Stringency 2") ## names of covariates in order

x = array(0,c(t,n,d))  
for (i in 1:n)
{
	x[,i,1] = data$stringency[((i-1)*t + 1):(i*t)]
	x[,i,2] = data$stringency2[((i-1)*t + 1):(i*t)]
}

FitBandsPlot(x,y,t,n.basis,covnames)

}

