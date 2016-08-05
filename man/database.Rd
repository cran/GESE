\name{database}
\alias{database}
\docType{data}
\title{
database file in example
}
\description{
a data frame containing the GENE and MAF information for the variants under consideration in the public reference database.
}
\usage{data("database")}
\format{
  A data frame of 20 observations on the following 3 variables.
  \describe{
    \item{\code{SNP}}{an unique identifier for variant}
    \item{\code{GENE}}{a character vector: Gene name}
    \item{\code{MAF}}{a numeric vector: minor allele frequency of the variants in the referecne database}
  }
}
\details{
A data frame containing the information for all the variants satisfying the same filtering criteria in the chosen reference genome. It should include at least three columns with these names: SNP (unique SNP ID), GENE (gene name), MAF (minor allele frequency for the variant in reference database for the corresponding population).
}
\source{
Randomly simulated data.
}
\examples{
data(database)
}
\keyword{datasets}
