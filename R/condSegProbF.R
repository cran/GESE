condSegProbF <-
function(pedTemp, subjInfo)
{   
  famInfo <- data.frame(pedTemp)
  ## founders that are related to the sequenced subjects
  subjInFam <- subjInfo[subjInfo$FID==pedTemp$famid[1],]
  ### depth of each subject in the pedigree
  famInfo$dv <- kindepth(famInfo$id, famInfo$dadid, famInfo$momid, align=TRUE )
  famInfo <- merge(famInfo, subjInFam, by.x="id", by.y="IID", sort=FALSE, all.x=TRUE)

  if(any(!subjInFam$PHENOTYPE %in% c(1,0, NA)))
  {	 stop("Error: subject phenotype can only be 0, 1 or NA \n")
  }
  caseFam <- famInfo[famInfo$id %in% subjInFam$IID[subjInFam$PHENOTYPE==1],]
  controlFam <- famInfo[famInfo$id %in% subjInFam$IID[subjInFam$PHENOTYPE==0],]
  if(nrow(caseFam)==1&nrow(controlFam)==0)
  {
    return(1)
  }
  
  commonFounders <- getFounder(pedTemp, subjInfo)  
  probA <- 0
  if(all(is.na(commonFounders)))
  {  founderCases <- famInfo[famInfo$PHENOTYPE==1 &!is.na(famInfo$PHENOTYPE) & famInfo$dadid==0 & famInfo$momid==0,]
 	 nonfounderCases <- famInfo[famInfo$PHENOTYPE==1 &!is.na(famInfo$PHENOTYPE) & famInfo$dadid!=0 & famInfo$momid!=0,]
 	 subjInfoTemp <- subjInfo[!subjInfo$IID %in% founderCases$id,]  
	 commonFoundersTemp <- getFounder(pedTemp, subjInfoTemp)
  
     founderCaseRelated <- sapply(founderCases$id, isRelated, commonFoundersTemp, famInfo)
     if(length(founderCaseRelated)>0)
     {
		 founderCaseRelated <- founderCases[founderCaseRelated,]
		 if(nrow(founderCaseRelated)>0)
		 {
			 for(i in 1:nrow(founderCaseRelated))
			 {	### remove all the other founder case
				 caseFamTemp <- rbind(nonfounderCases, founderCaseRelated[i,])
				 famInfoTemp <- famInfo
				famInfoTemp[famInfoTemp$id %in% founderCases$id & famInfoTemp$id != founderCaseRelated$id[i],]$PHENOTYPE  = NA
				max_dv <- max(caseFamTemp$dv, controlFam$dv)
				 #edit 10/31/2016
				 #probA <- probA + (1/nrow(founderCaseRelated))*getTranProb_dv(caseFamTemp, controlFam, famInfoTemp, founderCaseRelated$id[i], max_dv)
				 probA <- probA + getTranProb_dv(caseFamTemp, controlFam, famInfoTemp, founderCaseRelated$id[i], max_dv)

			}
		 }else
		 {	## remove those founder case, only use commonFounders
			 caseFamTemp <- nonfounderCases
			 famInfoTemp <- famInfo
			 famInfoTemp[famInfoTemp$id %in% founderCases$id,]$PHENOTYPE=NA
			 max_dv <- max(caseFamTemp$dv, controlFam$dv)
			 for(i in 1:length(commonFoundersTemp))
			 {	commonFounder <- commonFoundersTemp[i]
				 max_dv <- max(caseFamTemp$dv, controlFam$dv)
				 #edit 10/31/2016
				 #probA <- probA+(1/length(commonFoundersTemp))*getTranProb_dv(caseFamTemp, controlFam, famInfoTemp, commonFounder, max_dv)
				 probA <- probA + getTranProb_dv(caseFamTemp, controlFam, famInfoTemp, commonFounder, max_dv)
			 
			 }
		 }
  	 }else
  	 {	return (NA)
  	 }
     
  }else
  {	 for(i in 1:length(commonFounders))
  	 { 
   	 	commonFounder <- commonFounders[i]
   	 	max_dv <- max(caseFam$dv, controlFam$dv)
   	 	##edit 10/31/2016
   	 	#probA <- probA+(1/length(commonFounders))*getTranProb_dv(caseFam, controlFam, famInfo, commonFounder, max_dv)
    	probA <- probA + getTranProb_dv(caseFam, controlFam, famInfo, commonFounder, max_dv)
    
 	 }
  
  
  }
  
  return(probA)
  
}
