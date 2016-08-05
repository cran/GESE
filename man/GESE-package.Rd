\name{GESE-package}
\alias{GESE-package}
\docType{package}
\title{
\packageTitle{GESE}
}
\description{
\packageDescription{GESE}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{GESE}
\packageIndices{GESE}
computes gene-based segregation tests(GESE and weighted GESE) for family-based sequencing data. The main functions are:
\code{GESE}: computes gene-based segregation information and GESE test p-values (unweighted and weighted version).
\code{trim_oneLineage}: trims the pedigree so that for any subject, either the paternal family or the maternal family is included. Minimal set of sequenced subjects may be removed to ensure one lineage per pedigree only.
\code{trim_unrelated}: trims the pedigree so that only one founder case is kept for each pedigree, and pedigrees with no cases are removed.
\code{condSegProbF}: computes the conditional probability that a variant in the gene is segregating in the family specified, conditional on that the variant is present in the family.

}
\author{
\packageAuthor{GESE}

Maintainer: \packageMaintainer{GESE}
}
\references{
Qiao, D. Lange, C., Laird, N.M., Won, S., Hobbs, B., et al. 2016. {Gene-based segregation method for identifying rare variants for family-based sequencing studies.}

\url{http://scholar.harvard.edu/dqiao/gese}

}
\keyword{package}
\seealso{
\code{\link[GESE:GESE-package]{GESE}}
}
\examples{
data(pednew)
data(mapInfo)
data(dataRaw)
data(database)
results <- GESE(pednew, database, dataRaw, mapInfo, threshold=1e-2)
results

}
