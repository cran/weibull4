\name{likelihood}
\alias{likelihood}
\title{Calculates the Log Likelihood
}
\description{Likelihood is an Internal Function of the Weibull4 Package
}
\usage{
likelihood(x, y, param, modes)
}
\arguments{
  \item{x}{
Vector with the x values
}
  \item{y}{
Vector with the y values
  }
  \item{param}{
Vector with shape, scale, location, area and SD parameters for calculating the log of Likelihood for the weibull4 package
}
  \item{modes}{
Sets unimodal (modes=1) or bimodal (modes=2) Weibull's distribution
}
}
\value{
Unitary vector with the sum of the likelihood}
\references{
https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
}
\author{
Florian Hartig - Theoretical Ecology
}
\examples{
function (param)
{
    shape <- param[1]
    scale <- param[2]
    loc <- param[3]
    area <- param[4]
    sd <- param[5]
    pred <- weibull4(xi, shape, scale, loc, area)
    singlelikelihoods <- dnorm(yi, mean = pred, sd = sd, log = T)
    sumll <- sum(singlelikelihoods, na.rm = T)
    return(sumll)
  }
}
\keyword{ likelihood }
\keyword{ weibull4 }
