trim_oneLineage <-
function(seqSub, pednew)
{ #colnames of pednew)[1:5] include c("FID", "IID", "faID", "moID", "sex")
  #colnames of seqSub)[1:3] include c("FID" "IID", "PHENOTYPE")
  pednew <- pednew[pednew$FID %in% unique(seqSub$FID),]
  pednew$AFFECT <- ifelse(pednew$IID %in% seqSub$IID[seqSub$PHENOTYPE==1&!is.na(seqSub$PHENOTYPE)], 1, ifelse(pednew$IID %in% seqSub$IID[seqSub$PHENOTYPE==0&!is.na(seqSub$PHENOTYPE)], 0, NA))
  fams = unique(seqSub$FID)
  nped <- length(fams)
  pedigreeInfo <-c()
  seqSubRecord <- c()
  for(i in 1:nped)
  { cat("Family ", fams[i], "\n")
    seqSubFam = seqSub[seqSub$FID == fams[i], ]
    pedFam = pednew[pednew$FID==fams[i],]
    if(nrow(seqSubFam)==1)
    {	pedigreeInfoTemp = pedFam[pedFam$IID == seqSubFam$IID,]
    	pedigreeInfoTemp$faID = NA
        pedigreeInfoTemp$moID = NA
    	pedigreeInfo <-rbind(pedigreeInfo, pedigreeInfoTemp)
        seqSubRecord <- rbind(seqSubRecord, seqSubFam)
		next
    } 

    ### look at the couples that are both founders and establish their lineage, if all sequenced subjects are in either lineaage, just report that lineage, and ignore other lineage
    foundersFam = pedFam[is.na(pedFam$faID) & is.na(pedFam$moID), ]
    couplesSet = pedFam[,c("faID", "moID")]
    cSetFounder = couplesSet[couplesSet$faID %in% foundersFam$IID & couplesSet$moID %in% foundersFam$IID,]
    colnames(cSetFounder) = c("ID1", "ID2")
    temp = cSetFounder[,c(2,1)]
    colnames(temp) = c("ID1", "ID2")
    
    cset = unique(rbind(cSetFounder, temp))
    parentRecord = list()
    ### find extended parent set
    for(j in 1:nrow(cset))
    { parentTemp = as.character(cset[j,1:2])
      temp = cset[cset$ID1 %in% parentTemp, ]$ID2
      addParent = temp[!temp%in%parentTemp]
      parentTemp = unique(c(parentTemp, addParent))
      while(length(addParent)>0)
      {   
          temp = cset[cset$ID1 %in% addParent, ]$ID2
          addParent = temp[!temp%in%parentTemp]
          parentTemp = unique(c(parentTemp, addParent))
        
      }
      parentRecord[[j]] = sort(parentTemp)
      
    }
    parentRecord2 = unique(parentRecord)
   
    nlineage = length(parentRecord2)
    lineageV = list()
    for(nl in 1:nlineage)
    {   
        parents = parentRecord2[[nl]]
        children = pedFam[!pedFam$IID %in% parents & (pedFam$moID %in% parents | pedFam$faID %in% parents),]
        ## find partners of these children
        partners = c(pedFam[pedFam$moID %in% children$IID,]$faID, pedFam[pedFam$faID %in% children$IID,]$moID)
        addIn  <- unique(c(children$IID, partners))
        parents = unique(c(parents, addIn))
        while(length(addIn)!=0)
        { 
          children = pedFam[!pedFam$IID %in% parents & (pedFam$moID %in% children$IID | pedFam$faID %in% children$IID),]
          ## find partners of these children
          partners = c(pedFam[pedFam$moID %in% children$IID,]$faID, pedFam[pedFam$faID %in% children$IID,]$moID)
          addIn  <- unique(c(children$IID, partners))
          parents = unique(c(parents, addIn))
          
        }
        ## add in the partners of 
        lineageV[[nl]] = parents
        
    }
    nSegLineage = sapply(lineageV, function(x) sum(seqSubFam$IID %in% x))
    
    ### if not all sequenced subjects are in one lineage, only record the lineage with the most sequenced subjects, and remove the other sequenced subjects from the pedigree
    if(any(nSegLineage == nrow(seqSubFam)))
    {   pedigreeInfoTemp = pedFam[pedFam$IID %in% lineageV[[which(nSegLineage==nrow(seqSubFam))[1]]],]
        index0 = which(((!pedigreeInfoTemp$faID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$faID))|((!pedigreeInfoTemp$moID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$moID)))
        if(length(index0)>0)
        {
           pedigreeInfoTemp[index0,]$faID = NA
           pedigreeInfoTemp[index0,]$moID = NA
        }
       
        pedigreeInfo <-rbind(pedigreeInfo, pedigreeInfoTemp)
        seqSubRecord <- rbind(seqSubRecord, seqSubFam)
        #cat("All sequenced subjects are in the same lineage \n")
      
    }else
    {   ### check if  all sequenced cases are in one lineage,
    	casesID = seqSubFam[seqSubFam$PHENOTYPE == 1 &!is.na(seqSubFam$PHENOTYPE),]
    	nSegLineageCase = sapply(lineageV, function(x) sum(casesID$IID %in% x))
    	if(any(nSegLineageCase == nrow(casesID)))
    	{	pedigreeInfoTemp = pedFam[pedFam$IID %in% lineageV[[which(nSegLineageCase==nrow(casesID))[1]]],]
        	index0 = which(((!pedigreeInfoTemp$faID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$faID))|((!pedigreeInfoTemp$moID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$moID)))
        	if(length(index0)>0)
        	{
           		pedigreeInfoTemp[index0,]$faID = NA
           		pedigreeInfoTemp[index0,]$moID = NA
        	}
        	pedigreeInfo <-rbind(pedigreeInfo, pedigreeInfoTemp)
        	seqSubRecord <- rbind(seqSubRecord, seqSubFam[seqSubFam$IID %in% pedigreeInfoTemp$IID,])
       		cat("All sequenced cases are in the same lineage \n")
    	    cat("Exclude ", sum(!seqSubFam$IID %in% pedigreeInfoTemp$IID), " control subjects! \n")
    	}else
    	{
    		### if not all sequenced cases are in one lineage, only record the lineage with the most sequenced cases, and remove the other sequenced subjects from the pedigree
       	 	temp = order(nSegLineageCase, decreasing=TRUE)
        	largestLineage = lineageV[[temp[1]]]
        	pedigreeInfoTemp = pedFam[pedFam$IID %in% largestLineage,]
        	index0 = which(((!pedigreeInfoTemp$faID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$faID))|((!pedigreeInfoTemp$moID %in% pedigreeInfoTemp$IID) & !is.na(pedigreeInfoTemp$moID)))
         	if(length(index0)>0)
        	{
         	  pedigreeInfoTemp[index0,]$faID = NA
         	  pedigreeInfoTemp[index0,]$moID = NA
        	}
       
        	pedigreeInfo <-rbind(pedigreeInfo,  pedigreeInfoTemp)
       		seqSubRecord <- rbind(seqSubRecord, seqSubFam[seqSubFam$IID %in% largestLineage,])
        	#cat("Number of sequenced sujects in each lineage: ", nSegLineage, "\n")
        	cat("Exclude ", sum(!seqSubFam$IID %in% largestLineage), " subjects! \n")
        }
      
    }
    
    ### are there two founder cases in one pedigree
    foundersFamNew = pedigreeInfoTemp[is.na(pedigreeInfoTemp$faID) & is.na(pedigreeInfoTemp$moID), ]
    cases = seqSubRecord[seqSubRecord$AFFECT==1,]
    if(sum(cases$IID %in% foundersFamNew) > 1)
    {
      cat(sum(cases$IID %in% foundersFamNew)," founder cases in the same pedigree!!\n")
    }
    
    
    

  }
  return (list(pedInfoUpdate = pedigreeInfo,
  seqSubjUpdate = seqSubRecord))
}
