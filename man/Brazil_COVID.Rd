\name{Brazil_COVID}
\alias{Brazil_COVID}
\docType{data}
\title{Brazil's COVID-19 number of daily new cases and deaths
}
\description{Dataset for weibull4 package. This dataset contain the number of the Brazil's new daily cases and deaths for COVID-19 to be fitted to the 4-parameters Weibull distribution in weibull4 package.
}
\usage{data("Brazil_COVID")}
\format{
  A data frame with 349 observations on the following 3 variables.
  \describe{
    \item{\code{date}}{Brazil_COVID$date
    }
    \item{\code{new_cases}}{Brazil_COVID$new_cases
    }
    \item{\code{new_deaths}}{Brazil_COVID$new_deaths
    }
  }
}
\details{x values are in Date format. So, xmax must be as.Date() too.
}
\source{"https://covid.ourworldindata.org/data/owid-covid-data.csv"
}
\references{Roser M, Ritchie H, Ortiz-Ospina E, Hasel J. Coronavirus Pandemic (COVID-19) [Internet]. 2020 [cited 2020 May 15]. Available from: https://ourworldindata.org/coronavirus
}
\examples{
data(Brazil_COVID)
## maybe str(Brazil_COVID) ; plot(Brazil_COVID) ...
}
\keyword{datasets}
