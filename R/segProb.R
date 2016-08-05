segProb <-
function(pedTemp, subjInfo, variantMAF2, indicatorSeg, condSegProb)
{ 
  ### condSegProb will never be NA here
  if(is.na(condSegProb))
  { stop("ERROR: conditional segregation probability is NA!\n")
  }
  founders <- getFounder(pedTemp, subjInfo)
  numFounders <- length(founders) 
  DEFAULT_MAF= 0.0000001
  vmaf <- ifelse((is.na(variantMAF2$MAF)|variantMAF2$MAF==0), DEFAULT_MAF, variantMAF2$MAF)
  varInFam <- 1-(1-as.numeric(vmaf))^(2*numFounders)
  pedProb <- prod((1.0-condSegProb*varInFam))
  if(indicatorSeg==0)
  { ## no variant is segregating in the gene for the family
    return (pedProb)
  }else
  {
    return (1-pedProb)
  }
  
}
