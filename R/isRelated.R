isRelated <-
function(control, caseCommonFounder, famInfo)
{
	if(length(caseCommonFounder)==0)
	{	return (FALSE)
	}else if(all(is.na(caseCommonFounder)))
	{	return(FALSE)
	}
 	 ress <- c()
 	 for(i in 1:length(caseCommonFounder))
 	 {
 	   subjStudy <- famInfo[famInfo$id %in% c(control,caseCommonFounder[i]),]
 	   ress <- c(ress, findMostRecentCommonFounder(subjStudy, famInfo))
  	}
  	return(!all(is.na(ress)))
  
  
}
