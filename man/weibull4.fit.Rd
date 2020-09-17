\name{weibull4.fit}
\alias{weibull4.fit}
\title{Weibull 4-parameters Metropolis-MCMC non-linear curve fitting function
}
\description{Package to perform non-linear regression in data on the number of daily new
cases and daily new deaths of COVID-19 and other epidemics to the 4-parameters
Weibull distribution using Metropolis-Markov Chain-Monte Carlo Simulations (MCMC),
as described in Moreau, 2020
}
\usage{
weibull4.fit(x, y, shape=NA, scale=NA, loc=NA, area=NA, iter=1000, xmax=0)
}
\arguments{
  \item{x}{Vector: time data (may be date/time or numeric)
}
  \item{y}{Vector: observed/measure event
}
  \item{shape}{Starting value for Weibull's shape parameters. If it is NA,
  weibull4.fit will try to calculate it from x and y data.
}
  \item{scale}{Starting value for Weibull's scale parameters. If it is NA,
  weibull4.fit will try to calculate it from x and y data.
}
  \item{loc}{Starting value for Weibull's location parameters. If it is NA,
  weibull4.fit will try to calculate it from x and y data.
}
  \item{area}{Starting value for Weibull's area parameters or the area under the
  PDF curve. If it is NA, weibull4.fit will try to calculate it from x and y data.
}
  \item{iter}{Number of iterations to perform Metropolis-MCMC.
}
  \item{xmax}{Forecast date to be calculated after x data. It must be in the same format than x.
}
}
\details{This package was specially built to fit COVID-19 data on the number of
daily new cases and deaths in countries. So x must be integer. Alternatively,
Date format is allowed
}
\value{LIST containning:
  \item{Fit data}{Matrix with x and y fitted data}
  \item{Estimates}{Matrix containning shape, scale, location, area and SD of the MEtropolis-MCMC
  in the row 1; and standard deviation for each parameter in the row 2}
  \item{Markov chains}{Matrix containning the Markov chains for shape, scale, location,
  area and SD parameters}
}
\references{MOREAU, V. H. (2020) Forecast projections of COVID-19 pandemic by Weibull distribution from daily new cases and deaths data, submited.
Part of the MCMC code came from the Florian Hartig blog: "https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/"
}
\author{Vitor Hugo Moreau, Ph.D
}
\note{Comments, suggestions and doubts must be sent to vitorhmc@ufba.br
}

\section{Warning }{This package is a secondary product of the reffered science paper.
 Please, note that there is no warrants or professional support on its use}

\examples{
## Perform non-linear curve fitting with World's data for daily new deaths of COVID-19
fit <- weibull4.fit(worldCOVID$date, worldCOVID$new_deaths)
plot(worldCOVID$date, worldCOVID$new_deaths, ylab="World's daily new deaths", xlab="Date")
lines(fit[[1]][,1], fit[[1]][,2], col="red")
}
\keyword{ weibull }% use one of  RShowDoc("KEYWORDS")
