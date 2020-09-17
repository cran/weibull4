\name{proposalfunction}
\alias{proposalfunction}
\title{Proposal Distribution for Metropolis-MCMC
}
\description{
This is an internal function of the weibull4 package
}
\usage{
proposalfunction(param)
}
\arguments{
  \item{param}{Vector containing shape, scale, location, area and SD parameters
}
}
\value{
Vector containing proposal values for shape, scale, location, area and SD
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
    return(rnorm(5, mean = param, sd = param * 0.015))
  }
}
\keyword{ markov }
\keyword{ proposal }
