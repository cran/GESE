findMostRecentCommonFounder <-
function(caseFam, famInfo)
{   found <- FALSE
    if(nrow(caseFam)==1  & caseFam$dadid[1]!="0")
    {	parentsTemp <- c()
    }else
    {
    	parentsTemp <- caseFam$id
    }
    subjLook <- caseFam
    while(!found)
    {
      parentsTemp <- (c(parentsTemp, subjLook$dadid, subjLook$momid))
      parentsTemp <- parentsTemp[parentsTemp!='0']
      if(max(table(parentsTemp)[table(parentsTemp)!=0])==nrow(caseFam))
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
    return (as.character(commonFounders))
    
}
