\name{run_metropolis_MCMC}
\alias{run_metropolis_MCMC}
\title{
Runs the Metropolis-MCMC algorithm for weibull4 package
}
\description{
This is an internal function of the weibull4 package
}
\usage{
run_metropolis_MCMC(x, y, startvalue, iterations)
}
\arguments{
  \item{x}{
Vector with the x values
}
  \item{y}{
Vector with the y values
  }
  \item{startvalue}{
Vector with starting shape, scale, location, area and SD values for Metropolis-MCMC calculations
}
  \item{iterations}{
Number of iterations to be performed in MCMC simulation
}
}
\value{
Matrix with 5 columns and iterations rows with Markov chains for shape, scale, location, area and SD parameters)
}
\references{
https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
}
\author{
Florian Hartig - Theoretical Ecology
}
\examples{
function (startvalue, iterations)
{
    chain <- array(dim = c(iterations + 1, 5))
    chain[1, ] <- startvalue
    for (i in 1:iterations) {
        proposal <- proposalfunction(chain[i, ])
        probab <- exp(posterior(proposal) - posterior(chain[i,
            ]))
        if (runif(1) < probab) {
            chain[i + 1, ] <- proposal
        }
        else {
            chain[i + 1, ] <- chain[i, ]
        }
    }
    return(chain)
  }
}
\keyword{ metropolis }
\keyword{ markov }
