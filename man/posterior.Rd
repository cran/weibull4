\name{posterior}
\alias{posterior}
\title{
Calculates the posterior distribution for Metropolis-MCMC
}
\description{
This is an internal function of the weibull4 package
}
\usage{
posterior(x, y, param)
}
\arguments{
  \item{x}{
Vector with the x values
}
  \item{y}{
Vector with the y values
  }
  \item{param}{
Vector containning shape, scale, location, area and SD parameters
}
}
\value{
Vector containning the posterior distribution for Metropolis-MCMC
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
    return(likelihood(param) + prior(param))
  }
}
\keyword{ posterior }
\keyword{ metropolis }
