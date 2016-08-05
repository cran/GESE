getFounder <-
function(pedTemp, subjInfo)
{ famInfo <- data.frame(pedTemp)
  if(nrow(subjInfo)==0)
  {	stop("Error: no subject in the subjInfo matrix!\n")
  }
  
  subjInFam <- subjInfo[subjInfo$FID==pedTemp$famid[1],]
  famInfo$dv <- kindepth(famInfo$id, famInfo$dadid, famInfo$momid )
  caseFam <- famInfo[famInfo$id %in% subjInFam$IID[subjInFam$PHENOTYPE==1],]
  controlFam <- famInfo[famInfo$id %in% subjInFam$IID[subjInFam$PHENOTYPE==0],]
  if(nrow(caseFam)==1&nrow(controlFam)==0)
  {
    return(caseFam$id)
  }
  inLevel <- max(caseFam$dv)
  caseCommonFounder <- findMostRecentCommonFounder(caseFam, famInfo)
  if(nrow(controlFam)==0)
  {
    return(caseCommonFounder)
  }
  
  ### if no most recent common founder for all cases, return NA
  if(all(is.na(caseCommonFounder)))
  {
    return(NA)
  }
  offspringsMRCF <- caseCommonFounder
  thisLevelS <- caseCommonFounder
  findAllO <- FALSE
  while(!findAllO)
  { newOffspring <- famInfo[famInfo$dadid %in%thisLevelS  |famInfo$momid %in% thisLevelS , ]$id
    if(length(newOffspring)==0)
    {
      findAllO <- T
    }else
    {
      thisLevelS  <- newOffspring
      offspringsMRCF <- c(offspringsMRCF, thisLevelS)
    }
    
  }
  leftControl <- controlFam[!controlFam$id %in% offspringsMRCF,]
  if(nrow(leftControl)==0)
  {
    return(caseCommonFounder)
  }
  
  foundersFinal <- c()
  for(i in 1:length(caseCommonFounder))
  {
    leftControlRelated <- sapply(leftControl$id, isRelated, caseCommonFounder[i], famInfo)
    leftControl2 <- leftControl[leftControlRelated,] ## only care about controls that are related to the mrcf of the cases
    if(nrow(leftControl2)>0)
    { 	mrcf <- famInfo[famInfo$id == caseCommonFounder[i],]
        subjStudy <- rbind(leftControl2, mrcf)
        ultimateFounder <- findMostRecentCommonFounderControl(leftControl2, mrcf, famInfo)
        intermediateFounder <- findIntermediateFounder(mrcf, ultimateFounder[1], famInfo)
        foundersFinal <- unique(c(foundersFinal, intermediateFounder , ultimateFounder))
      
    }else
    {
      foundersFinal <- unique(c(foundersFinal, caseCommonFounder[i]))
    }
    
    
  }
  return(unlist(foundersFinal))
  
}
