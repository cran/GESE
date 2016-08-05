findIntermediateFounder <-
function(mrcf, founder, famInfo)
{
  ### find all the founders between mrcf and the real founder
  parents <- famInfo[famInfo$id == mrcf$id, c("momid", "dadid")]
  if(founder %in% parents)
  {
    return(as.character(parents[!parents %in% founder]))
  }
  pRelated <- sapply(parents, isRelated, founder, famInfo)
  moreFounder <- parents[!pRelated]
  mrcf <- famInfo[famInfo$id %in% parents[pRelated],]
  moreFounder <- c(moreFounder, findIntermediateFounder(mrcf, founder, famInfo))
  return(as.character(moreFounder))
  
}
