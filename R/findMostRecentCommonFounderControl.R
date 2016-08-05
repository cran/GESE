findMostRecentCommonFounderControl <-
function(controlFam, case1, famInfo)
{   found <- FALSE
    subjLook <- rbind(controlFam, case1)
    parentsTemp <- subjLook$id
    mindv <- min(subjLook$dv)
    mindvsub <- subjLook[subjLook$dv==mindv, ]
    numRepeat = nrow(controlFam)+1
    if(nrow(mindvsub)==2)
    { if( any(subjLook$dadid %in% mindvsub$id[mindvsub$sex=="male"] & subjLook$momid %in% mindvsub$id[mindvsub$sex=="female"] ))
      {
         numRepeat = numRepeat-1
       }
    }
    
    while(!found)
    {
      parentsTemp <- (c(parentsTemp, subjLook$dadid, subjLook$momid))
      parentsTemp <- parentsTemp[parentsTemp!='0']
      
        if(max(table(parentsTemp)[table(parentsTemp)!=0])==numRepeat)
      {
        found <- TRUE
        temp <- table(parentsTemp)
        commonFounders <- names(temp)[temp==max(temp)]
      
      }else
      { if(all(subjLook$dadid=="0") & all(subjLook$momid=="0"))
        {
         return (NA)
       }
      
      temp1 <- c(subjLook$dadid, subjLook$momid)
      temp1 <- temp1[temp1!='0']
      index1 <- sapply(temp1, function(x) which(famInfo$id ==x))
      subjLook <- famInfo[index1,]
      
      }
    }
    return (commonFounders)
    
}
