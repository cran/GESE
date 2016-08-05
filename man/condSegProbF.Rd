\name{condSegProbF}
\alias{condSegProbF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes conditional segregation probability for the family
}
\description{
Computes the conditional probability that a variant is segregating in the family conditional on that the variant is present in one of the founders in the family.
}
\usage{
condSegProbF(pedTemp, subjInfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedTemp}{
The data frame that includes the complete pedigree structure for the family}
  \item{subjInfo}{
A data frame that contains the subject phenotype information for the sequenced subjects. it should include the columns FID, IID, and PHENOTYPE.}
}

\value{
returns the conditional segregating probability of a variant in the family
}
\references{
Qiao, D. Lange, C., Laird, N.M., Won, S., Hobbs, B., et al. 2016. {Gene-based segregation method for identifying rare variants for family-based sequencing studies.}
}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GESE}}
}
\examples{
data(pednew)
data(mapInfo)
data(dataRaw)
data(database)
library(kinship2)
pedigrees = kinship2::pedigree(pednew$IID, pednew$faID, pednew$moID,pednew$sex,famid=pednew$FID)
subjects= dataRaw[,c(1,2,6)]
condSegProbF(pedigrees['93'], subjects)
condSegProbF(pedigrees['412'], subjects)
results2 <- GESE(pednew, database, dataRaw, mapInfo, threshold=1e-2)
results2$condSegProb
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{probability}
