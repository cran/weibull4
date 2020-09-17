\name{weibull4}
\alias{weibull4}
\title{weibull 4-parameters distribution building function
}
\description{Builds data for a 4-parameters Weibull distribution of a given x data
}
\usage{
weibull4(x=seq(0, 1, length.out = 10), shape=2.5, scale=1, loc=0, area=1)
}
\arguments{
  \item{x}{Vector: data range for calculation of the Weibull distribution. If it is NULL, it will be set to seq(0,1,0.1)
}
  \item{shape}{Weibull's shape parameter
}
  \item{scale}{Weibull's scale parameter
}
  \item{loc}{Weibull's location parameter
}
  \item{area}{Weibull's area parameter: area under the PDF curve
}
}
\details{This package was specially built to fit COVID-19 data on the number of
daily new cases and deaths in countries. So x must be integer. Alternatively,
Date format is allowed
}
\value{Vector: f(x) Weibull distribution's ordinate
}
\references{MOREAU, V. H. (2020) Forecast projections of COVID-19 pandemic by Weibull distribution from daily new cases and deaths data, submited.
}
\author{Vitor Hugo Moreau, Ph.D.
}
\note{Comments, suggestions and doubts must be sent to vitorhmc@ufba.br
}

\section{Warning }{This package is a secondary product of the reffered science paper. Please, note that there is no warrants or professional support on its use}

\examples{
## Build a 4-parameters Weibull distribution with given parameters
weibull4(seq(1,100,1), shape=2.5, scale=30, loc=10, area=1000)
## Build and plot 4-parameters Weibull distribution with given parameters
plot(seq(1,100,1), weibull4(seq(1,100,1), 2, 30, 10, 1), type="l")
## Build and plot 4-parameters Weibull distribution with a time series in the abscissa
Date <- seq(Sys.Date(), as.Date("2022-12-31"),1)
plot(Date, weibull4(Date, 1.6, 100, 100, 100), type="l")
}
\keyword{ weibull }
