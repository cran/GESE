\name{pednew}
\alias{pednew}
\docType{data}
\title{
pednew - an example pedigree structure
}
\description{
A data frame of the complete pedigree strucutre for the families included
}
\usage{data("pednew")}
\format{
  A data frame of 1700 observations on the following 26 variables.
  \describe{
    \item{\code{FID}}{Family ID of class character}
    \item{\code{IID}}{Individual ID of class character}
    \item{\code{faID}}{Father ID, NA if missing}
    \item{\code{moID}}{Mother ID, NA if missing}
    \item{\code{sex}}{Sex, 1 for male, 2 for female and NA if missing.}
  }
}

\examples{
data(pednew)
}
\keyword{datasets}
