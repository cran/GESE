trim_unrelated <- function(seqSub, pednew2)
{	pedigrees = kinship2::pedigree(id=pednew2$IID, dadid=pednew2$faID, momid = pednew2$moID, sex=pednew2$sex, famid=pednew2$FID)
	fam <- as.character(unique(seqSub$FID))
	numFam <- length(fam)

	for(i in  1:numFam)
	{	cat("Family ", fam[i], "\n")
	    pedTemp <- pedigrees[as.character(fam[i])]
	    famInfo <- data.frame(pedTemp)
	    subjFam <- seqSub[seqSub$FID == fam[i],]
	    if(nrow(subjFam)==1)
	    {	cat("--do not need to trim\n")
	    	next
	    }
	    famInfo$dv <- kindepth(famInfo$id, famInfo$dadid, famInfo$momid, align=TRUE )
	    subjFam2 <- merge(subjFam, famInfo, by.x="IID", by.y="id", sort=FALSE)
		
		commonFounders <- getFounder(pedTemp, subjFam2)
		founderCases <- subjFam2[subjFam2$dadid==0&subjFam2$momid==0&subjFam2$PHENOTYPE==1,]
 	    nonfounderCases <- subjFam2[subjFam2$PHENOTYPE==1 &!is.na(subjFam2$PHENOTYPE) & subjFam2$dadid!=0 & subjFam2$momid!=0,]
 	    ### non founder controls
 	 	controls <- subjFam2[subjFam2$PHENOTYPE==0 &!is.na(subjFam2$PHENOTYPE)& subjFam2$dadid!=0 & subjFam2$momid!=0,] 
       ### if there are two founder cases
  	   if(anyNA(commonFounders) | nrow(founderCases)>1 )
  		{   subjInfoTemp <- subjFam2[!subjFam2$IID %in% founderCases$IID,] 
  		    ### if no nonfounder case left, keep only one founder
  		    if(sum(subjInfoTemp$PHENOTYPE==1)==0 & nrow(founderCases)>1)
  		    {		#### if no nonfounder cases left and there are non-founder controls, select the founder related to the controls
  		    		if(nrow(controls)>0)
  		    		{	founderCaseRelatedControl <- sapply(founderCases$IID, isRelated,  controls, famInfo)
  		    		    founderRelateContrl <- founderCases[founderCaseRelatedControl,]
  		    		    if(nrow(founderRelateContrl)>0)
  		    		    {	leaveFounder <- founderRelateContrl[1,]
  		    		    }else
  		    		    {	leaveFounder <-founderCases[order(founderCases$dv),][1,]
  		    		    }
  		    		}else
  		    		{
  		    			### if no non-founder controls, select the founder at the highest level
  		    
  		           	 leaveFounder <- founderCases[order(founderCases$dv),][1,]
  		           	 }
  					### remove all the other founder case from dataPed
  	    			seqSub <- seqSub[!seqSub$IID %in% founderCases$IID | seqSub$IID == leaveFounder$IID,]
  	     
  		    }else
  		    {
  		    ##  there are non founder cases
 	 		commonFoundersTemp <- getFounder(pedTemp, subjInfoTemp)
  			founderCaseRelated <- sapply(founderCases$IID, isRelated, commonFoundersTemp, famInfo)
    	    ### if any founder case is unrelated to the non-founders common ancestors, remove these founder cases
    	    ### then only work on founder cases that are related to the non-founders
     		founderCaseRelated <- founderCases[founderCaseRelated,]
     		if(nrow(founderCaseRelated)>0)
     		{   ### if any founder is related to the common ancestors of the non-founder cases
     		    ### choose only the founder with highest dv
     		    leaveFounder <- founderCaseRelated[order(founderCaseRelated$dv),][1,]
  				### remove all the other founder case from dataPed
  	    		seqSub <- seqSub[!seqSub$IID %in% founderCases$IID | seqSub$IID == leaveFounder$IID,]
  	     	}else
			{	### no founder is related to the common ancestors of the non-founders, then remove all these founders
				seqSub <- seqSub[!seqSub$IID %in% founderCases$IID,]
			}
			}
		}else
		{	cat("--do not need to trim\n")
		}	
	}
	#### remove families without a case
	caseFam = unique(seqSub[seqSub$PHENOTYPE==1,]$FID)
	famNoCase = seqSub[!seqSub$FID %in% caseFam,]
	seqSub2 <- seqSub[!seqSub$FID %in% famNoCase$FID,]
	cat("Removed families with no case subjects:",  unique(famNoCase$FID), "\n")
	return(seqSub2)
}