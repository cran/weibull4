\name{weibull4.build}
\alias{weibull4.build}
\title{weibull 4-parameters distribution building function
}
\description{Builds data for a 4-parameters Weibull distribution of a given x data
}
\usage{
weibull4.build(x=seq(0,1,length.out=10), shape=2.5, scale=1, loc=0, area=20,
shape2=5, scale2=2, loc2=6, area2=1, modes=1)
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
  \item{shape2}{second mode Weibull's shape parameter. It works only if modes=2.
}
  \item{scale2}{second mode Weibull's scale parameter. It works only if modes=2.
}
  \item{loc2}{second mode Weibull's location parameter. It works only if modes=2.
}
  \item{area2}{second mode Weibull's area parameter. It works only if modes=2.
}
  \item{modes}{Sets whether data may be fit with unimodal or bimodal Weibull distribution. This option was implemented to fit data feom the second wave of infections and deaths for COVID-19. Use modes=1 for unimodal distribution (single peak) and modes=2 for bimodal distribution (two peaks).
}

}
\details{This package was specially built to fit COVID-19 data on the number of
daily new cases and deaths in countries. So x must be integer. Alternatively,
Date format is allowed
}
\value{Vector: f(x) Weibull distribution's ordinate
}
\references{MOREAU, V. H. (2021) Using the Weibull distribution to model COVID-19 epidemic data. Modeling Assited Statistics Applications, 2021, in press.
}
\author{Vitor Hugo Moreau, Ph.D.
}
\note{Comments, suggestions and doubts must be sent to vitorhmc@ufba.br
}

\section{Warning }{This package is a secondary product of the reffered science paper. Please, note that there is no warrants or professional support on its use}

\examples{
## Build a 4-parameters Weibull distribution with given parameters
weibull4.build(seq(1,100,1), shape=2.5, scale=30, loc=10, area=1000, modes=1)
## Build and plot 4-parameters Weibull distribution with given parameters
plot(seq(1,100,1), weibull4.build(seq(1,100,1), 2, 30, 10, 1), type="l")
## Build and plot 4-parameters Weibull distribution with a time series in the abscissa
Date <- seq(Sys.Date(), as.Date("2022-12-31"),1)
plot(Date, weibull4.build(Date, 1.6, 100, 100, 100), type="l")
## Build and plot a bimodal, 4-parameters Weibull distribution with given parameters
weibull4.build(seq(1,100,1), shape=2.5, scale=30, loc=10, area=1000,
shape2=2.5, scale2=60, loc2=40, area2=1000, modes=2)
}
\keyword{ weibull.build }
