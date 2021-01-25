\name{weibull4.fit}
\alias{weibull4.fit}
\title{Weibull 4-parameters Metropolis-MCMC non-linear curve fitting function
}
\description{Package to perform non-linear regression in data on the number of daily new cases and daily new deaths of COVID-19 and other epidemics to the 4-parameters
Weibull distribution using Metropolis-Markov Chain-Monte Carlo Simulations (MCMC),
as described in Moreau, 2021
}
\usage{
weibull4.fit(x, y, shape=NA, scale=NA, loc=NA, area=NA, shape2=NA, scale2=NA, loc2=NA,
area2=NA, iter=1000, xmax=0, modes=1)
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
  \item{shape2}{Starting value for the shape parameters of the second mode of the
  Weibull's distribution. If it is NA, weibull4.fit will try to calculate it from x and y data. It works only if modes=2.
}
  \item{scale2}{Starting value for the scale parameters of the second mode of the
  Weibull's distribution. If it is NA, weibull4.fit will try to calculate it from x and y data. It works only if modes=2.
}
  \item{loc2}{Starting value for the location parameters of the second mode of the
  Weibull's distribution. If it is NA, weibull4.fit will try to calculate it from x and y data. It works only if modes=2.
}
  \item{area2}{Starting value for the area parameters of the second mode of the
  Weibull's distribution. If it is NA, weibull4.fit will try to calculate it from x and y data. It works only if modes=2.
}
  \item{iter}{Number of iterations to perform Metropolis-MCMC.
}
  \item{modes}{Sets whether data may be fit with unimodal or bimodal Weibull distribution. This option was implemented to fit data feom the second wave of infections and deaths for COVID-19. Use modes=1 for unimodal distribution (single peak) and modes=2 for bimodal distribution (two peaks).
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
\references{MOREAU, V. H. (2021) Forecast projections of COVID-19 pandemic by Weibull distribution from daily new cases and deaths data, Model Assisted Statistics and Applications. in press.
Part of the MCMC code came from the Florian Hartig blog: "https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/"
}
\author{Vitor Hugo Moreau, Ph.D
}
\note{Comments, suggestions and doubts must be sent to vitorhmc@ufba.br
}

\section{Warning }{This package is a secondary product of the reffered science paper.
 Please, note that there is no warrants or professional support on its use}

\examples{
## Perform non-linear curve fitting with World's data for daily new deaths
## of COVID-19
## Examples below are with low number of iterations (1000), because of CRAN rules limitations.
## They will possibily give poor results. For best results set iter=10000.
fit <- weibull4(US_COVID$date, US_COVID$new_deaths, iter=1000)
plot(US_COVID$date, US_COVID$new_deaths, ylab="US' daily new deaths", xlab="Date")
lines(fit[[1]][,1], fit[[1]][,2], col="red")
## Perform non-linear curve fitting with Canada's data for daily new cases of
## COVID-19 using two unimodal Weibull distribution, with split date on Aug,  ## 1st
fit <- weibull4(Canada_COVID$date, Canada_COVID$new_cases, modes=1, iter=1000)
plot(Canada_COVID$date, Canada_COVID$new_cases, ylab="Canada's daily new cases", xlab="Date")
lines(fit[[1]][,1], fit[[1]][,2], col="red")
## Perform non-linear curve fitting with Brazil's data for daily new deaths of
## COVID-19 using a single bimodal Weibull distribution
fit <- weibull4(Brazil_COVID$date, Brazil_COVID$new_deaths, modes=2, iter=1000)
plot(Brazil_COVID$date, Brazil_COVID$new_deaths, ylab="Brazil's daily new deaths", xlab="Date")
lines(fit[[1]][,1], fit[[1]][,2], col="red")
}
\keyword{ weibull.fit }% use one of  RShowDoc("KEYWORDS")
