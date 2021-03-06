\name{dataRaw}
\alias{dataRaw}
\docType{data}
\title{
dataRaw - a data frame containing the pedigree, phenotype and genotype information
}
\description{
A data frame that can be created from the \code{.raw} formatted filed generated by PLINK. 
}
\usage{data("dataRaw")}
\format{
  A data frame with 198 observations on the following 26 variables.
  \describe{
    \item{\code{FID}}{Family iD}
    \item{\code{IID}}{Individual ID}
    \item{\code{PAT}}{Father ID}
    \item{\code{MAT}}{Mother ID}
    \item{\code{SEX}}{sex}
    \item{\code{PHENOTYPE}}{Affection status}
    \item{\code{X1}}{Genotype for variant 1}
    \item{\code{X2}}{Genotype for variant 2}
    \item{\code{X3}}{Genotype for variant 3}
    \item{\code{X4}}{Genotype for variant 4}
    \item{\code{X5}}{Genotype for variant 5}
    \item{\code{X6}}{Genotype for variant 6}
    \item{\code{X7}}{Genotype for variant 7}
    \item{\code{X8}}{Genotype for variant 8}
    \item{\code{X9}}{Genotype for variant 9}
    \item{\code{X10}}{Genotype for variant 10}

  }
}
\details{
 The number of rows equal the number of subjects in the data and the number of columns equas the number of markers M + 6. The first six columns with specific column names include: the Family ID (FID), Individual ID (IID), father ID(PAT), mother ID (MAT), sex (SEX) and affection status (PHENOTYPE). The rest of the columns containing the genotypes for the variants listed in the coreesponding \code{mapInfo} file. It is also important to make sure that the recoding is with respect to the minor allele in the population. The affection status of this file will be used as the phenotype.
}

\examples{
data(dataRaw)
}
\keyword{datasets}
