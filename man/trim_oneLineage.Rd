\name{trim_oneLineage}
\alias{trim_oneLineage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trims the pedigree structure to include one lineage only.
}
\description{
Trims the families to include only one lineage.
}
\usage{
trim_oneLineage(seqSub, pednew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqSub}{
  A data frame that should include three columns FID (family ID), IID (individual ID), and PHENOTYPE (affection status) for the sequenced subjects in the data. One example is the 1st, 2nd and 6th columns from the plink raw format.
}
  \item{pednew}{
A data frame includes the complete pedigree structure information for all sequenced families in the dataset. The required column names of this data frame include: FID (family ID), IID (individual ID, must be of class character), faID (father ID, NA if unavailable), moID (mother ID, NA if unavailable), and sex.

}
}
\details{
For each subject, only the maternal or the paternal family is included, since the rare variant should be present in only the related subjects. The lineage with the maximal set of sequenced cases will be used as the final pedigree.
}
\value{
\item{pedInfoUpdate}{the complete pedigrees with only the paternal or maternal lineage }
\item{seqSubjUpdate}{The sequenced subjects that are in the selected lineage are returned for the rest of the analysis.}
}
\references{
Qiao, D. Lange, C., Laird, N.M., Won, S., Hobbs, B., et al. 2016. {Gene-based segregation method for identifying rare variants for family-based sequencing studies.}

}
\author{
Dandi Qiao
}
\note{
This function can be used for other analysis of family-based data processing. For example, the pre-processing step for PVAAST analysis.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GESE}}, \code{\link{trim_unrelated}}}
\examples{
data(pednew)
data(mapInfo)
data(dataRaw)
data(database)
subjects <- dataRaw[,c(1:2, 6)]
cat("Trimming the families...\n")
cat("Trimming step 1: keep only one lineage \n")
trim <- trim_oneLineage(seqSub=subjects, pednew)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{pedigree}
