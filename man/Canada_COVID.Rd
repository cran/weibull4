\name{Canada_COVID}
\alias{Canada_COVID}
\docType{data}
\title{Canada's COVID-19 number of daily new cases and deaths
}
\description{Dataset for weibull4 package. This dataset contain the number of the Canada's new daily cases and deaths for COVID-19 to be fitted to the 4-parameters Weibull distribution in weibull4 package.
}
\usage{data("Canada_COVID")}
\format{
  A data frame with 349 observations on the following 3 variables.
  \describe{
    \item{\code{date}}{Canada_COVID$date
    }
    \item{\code{new_cases}}{Canada_COVID$new_cases
    }
    \item{\code{new_deaths}}{Canada_COVID$new_deaths
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
data(Canada_COVID)
## maybe str(Canada_COVID) ; plot(Canada_COVID) ...
}
\keyword{datasets}
