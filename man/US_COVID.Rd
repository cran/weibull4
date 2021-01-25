\name{US_COVID}
\alias{US_COVID}
\docType{data}
\title{
US number of new COVID-19 cases
}
\description{
Dataset for weibull4 package. This dataset contain the number of the US' new daily cases of COVID-19 to be fitted to the 4-parameters Weibull distribution in weibull4 package with modes=2, in order to fit to the second wave of COVID-19 infections.
}
\usage{data("US_COVID")}
\format{
  A data frame with 336 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{US_COVID$date
    }
    \item{\code{new_cases}}{US_COVID$new_cases
    }
    \item{\code{new_deaths}}{US_COVID$new_deaths
    }
  }
}
\details{
x values are in Date format. So, xmax must be as.Date()
}
\source{
"https://covid.ourworldindata.org/data/owid-covid-data.csv"
}
\references{
Roser M, Ritchie H, Ortiz-Ospina E, Hasel J. Coronavirus Pandemic (COVID-19) [Internet]. 2020 [cited 2020 May 15]. Available from: https://ourworldindata.org/coronavirus
}
\examples{
data(US_COVID)
## maybe str(US_COVID) ; plot(US_COVID) ...
}
\keyword{datasets}
