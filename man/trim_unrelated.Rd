\name{trim_unrelated}
\alias{trim_unrelated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trims the pedigree structure to exclude multiple founder cases
}
\description{
Trims the families to include only one founder case in each pedigree. It also exclude families with only one control subject.
}
\usage{
trim_unrelated(seqSub, pednew2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqSub}{ 
A data frame that includes at least three columns: the Family ID (FID), Individual ID (IID), and affection status (PHENOTYPE). This input values should be output from the \code{trim_oneLineage} function.
}
  \item{pednew2}{
A data frame includes the complete pedigree structure information for all sequenced families in the dataset. The required column names of this data frame include: FID (family ID), IID (individual ID, must be of class character), faID (father ID, NA if unavailable), moID (mother ID, NA if unavailable), and sex.

}
}
\details{
For each pedigree,each there are multiple case founders in the pedigree, to satisfy the assumption that the causal rare variant must be introduced by one founder only, we need to keep only one case fonder that is related to most of the sequened subjects in the pedigree. We also want to remove families with no case.
}
\value{
The modified dataPed2 file after the trimming.
}
\references{
Qiao, D. Lange, C., Laird, N.M., Won, S.,  Hersh, C.P., et al. (2017). \emph{Gene-based segregation method for identifying rare variants for family-based sequencing studies.} Genet Epidemiol 41(4):309-319. DOI:10.1002/gepi.22037.

}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GESE}}, \code{\link{trim_oneLineage}}}
\examples{
data(pednew)
data(mapInfo)
data(dataRaw)
data(database)
subjects <- dataRaw[,c(1:2, 6)]
cat("Trimming the families...\n")
cat("Trimming step 1: keep only one lineage \n")
trim <- trim_oneLineage(seqSub=subjects, pednew)
subjects2 <- trim_unrelated(trim$seqSubjUpdate, trim$pedInfoUpdate)

	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pedigree}
