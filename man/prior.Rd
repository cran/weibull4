\name{prior}
\alias{prior}
\title{
Calculates the Prior Distribution for Metropolis-MCMC
}
\description{
This is an internal function of the Weibull4 package
}
\usage{
prior(param)
}
\arguments{
  \item{param}{
A vector with shape, scale, location, area and SD parameters
}
}
\value{
A vector with prior distribution for Metropolis-MCMC
}
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
    shapeprior <- dunif(shape, min = 1, max = 5, log = T)
    scaleprior <- dnorm(scale, sd = scale/2, log = T)
    locprior <- dunif(loc, min = 1, max = loc * 2, log = T)
    areaprior <- dunif(area, min = area/2, max = area * 2, log = T)
    sdprior <- dunif(sd, min = 1, max = sd * 2, log = T)
    return(shapeprior + scaleprior + locprior + areaprior + sdprior)
  }
}
\keyword{ weibull4 }
\keyword{ prior }
