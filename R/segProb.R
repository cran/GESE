segProb <-
function(pedTemp, subjInfo, variantMAF2, indicatorSeg, condSegProb, dbSize)
{ 
  ### condSegProb will never be NA here
  if(is.na(condSegProb))
  { stop("ERROR: conditional segregation probability is NA!\n")
  }
  founders <- getFounder(pedTemp, subjInfo)
  #numFounders <- length(founders) 
  DEFAULT_MAF= 1-(0.95)^(1/(2*dbSize))
  vmaf <- ifelse((is.na(variantMAF2$MAF)|variantMAF2$MAF==0), DEFAULT_MAF, variantMAF2$MAF)
  varInFam <- 1-(1-as.numeric(vmaf))^(2)
  pedProb <- prod((1.0-condSegProb*varInFam))
  if(indicatorSeg==0)
  { ## no variant is segregating in the gene for the family
    return (pedProb)
  }else
  {
    return (1-pedProb)
  }
  
}
