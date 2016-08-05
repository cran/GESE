getTranProb_dv <-
function(caseFam, controlFam, famInfo, commonFounder, max_dv)
{
  if(nrow(controlFam)>0)
  {
  	controlRelated <- sapply(controlFam$id, isRelated, commonFounder, famInfo)
  	controlFam <- controlFam[controlRelated,] 
  }

  studySubj <- rbind(caseFam, controlFam)
  studySubj_dv <- studySubj[studySubj$dv==max_dv,]
  while(nrow(studySubj_dv)==0)
  {
    max_dv=max_dv-1
    studySubj_dv <- studySubj[studySubj$dv==max_dv,]
  }
  dadInfo <- sort(table(studySubj_dv$dadid),decreasing=TRUE)
  momInfo <- sort(table(studySubj_dv$momid),decreasing=TRUE)
  ### only care about parents that are related to the founder
  dadRelate <- sapply(names(dadInfo), isRelated, commonFounder, famInfo)
  momRelate <- sapply(names(momInfo), isRelated, commonFounder, famInfo)
  dadInfo2 <- dadInfo[dadRelate]
  momInfo2 <- momInfo[momRelate]
  info2 <- c(dadInfo2, momInfo2)
  info2 <- info2[names(info2)!="0"]
 
  probBase <- 1
  comFounder_dv <- famInfo[famInfo$id==commonFounder,]$dv
  if(max_dv==comFounder_dv+1)
  {
    if(!all(names(info2) %in% commonFounder))
    {
      stop("ERROR: dv level 1, but not all parents are common founders!!\n")
      return(NA)
    }
    if(commonFounder %in% famInfo$id[famInfo$PHENOTYPE==0])
    {
      probBase <- 0
    }else
    {
      probBase <- (1/2)^info2[names(info2)==commonFounder]
    }
    
    return(probBase)
  }else  if(max_dv==comFounder_dv)
  {		if(studySubj_dv$id != commonFounder)
  		{	stop("ERROR: cases at the same level as common founder is not the common founder!")
  		}
  		return (1)
  		
  }
   
  allcontrols <- rep(NA, length(info2))
  for(i in 1:length(info2))
  {
    studySubjOS <- studySubj[studySubj$momid%in%names(info2)[i]|studySubj$dadid%in%names(info2)[i],]
    allcontrols[i] <- all(studySubjOS$id %in% controlFam$id)
  }
  
  probOneLevel <- 1 ## the probability of the offspring's status given variant introduced by the common founder
  newCase <- newControl <- c()
  for(i in 1:length(info2))
  {
    if(allcontrols[i]==FALSE)
    { 
      
      if(names(info2)[i] %in% controlFam$id)
      {
        cat("the related parent(to common founder) of at least one case in family ", famInfo$FID[1]," is a non-carrier!\n")
        probOneLevel <- NA
        
      }else
      {
        probOneLevel <- probOneLevel*((1/2)^(info2[i]))
        newCase <- c(newCase, names(info2)[i] )
      }
      
    }
  }

  probFinal <- probOneLevel
  subjConsider <- names(info2)[allcontrols]
  if(length(subjConsider)==0)
  { 
    newCaseFam <- rbind(famInfo[famInfo$id %in% newCase&!famInfo$id %in% caseFam$id,],caseFam[caseFam$dv<max_dv,])
    newControlFam <- controlFam[controlFam$dv<max_dv,]
    probFinal <- probFinal*getTranProb_dv(newCaseFam, newControlFam, famInfo, commonFounder, max_dv-1)
    return(probFinal)
  }else
  {
    

    if(any(subjConsider %in% controlFam$id))
    { newCaseFam <- rbind(famInfo[famInfo$id %in% newCase&!famInfo$id %in% caseFam$id,],caseFam[caseFam$dv<max_dv,])
      newControlFam <- controlFam[!(controlFam$momid %in% subjConsider[subjConsider %in% controlFam$id]| controlFam$dadid %in% subjConsider[subjConsider %in% controlFam$id]) & !controlFam$dadid %in% newCase & !controlFam$momid %in% newCase ,]
      probFinal <- probFinal*1*getTranProb_dv(newCaseFam, newControlFam, famInfo, commonFounder, max_dv)
      return (probFinal)
    }

    if(any(subjConsider %in% caseFam$id))
    { newCaseFam <- rbind(famInfo[famInfo$id %in% newCase&!famInfo$id %in% caseFam$id,],caseFam[caseFam$dv<max_dv,])
      newControlFam <- controlFam[!(controlFam$momid %in% subjConsider[subjConsider %in% caseFam$id]| controlFam$dadid %in% subjConsider[subjConsider %in% caseFam$id]) & !controlFam$dadid %in% newCase & !controlFam$momid %in% newCase,]
      numMei <- sum(info2[subjConsider[subjConsider %in% caseFam$id]])
      probFinal <- probFinal*(1/2)^numMei*getTranProb_dv(newCaseFam, newControlFam, famInfo, commonFounder, max_dv)
      return (probFinal)
    }
    

    newCaseFam <- rbind(famInfo[famInfo$id %in% newCase&!famInfo$id %in% caseFam$id  ,],caseFam[caseFam$dv<max_dv,])
    newControlFam <-  rbind(controlFam[!(controlFam$dadid%in%subjConsider[1]|controlFam$momid%in%subjConsider[1]) & !controlFam$dadid %in% newCase & !controlFam$momid %in% newCase,], famInfo[famInfo$id==subjConsider[1],])
    
    newCase2 <- c(newCase, subjConsider[1])
    newCaseFam2 <- rbind(famInfo[famInfo$id %in% newCase2 &!famInfo$id %in% caseFam$id,], caseFam[caseFam$dv<max_dv,])
    newControlFam2 <- controlFam[!(controlFam$dadid%in%subjConsider[1]|controlFam$momid%in%subjConsider[1])& !controlFam$dadid %in% newCase & !controlFam$momid %in% newCase,]
    
    probFinal <- probFinal*(((1/2)^info2[names(info2)==subjConsider[1]])*getTranProb_dv(newCaseFam2, newControlFam2, famInfo, commonFounder, max_dv)+getTranProb_dv(newCaseFam, newControlFam, famInfo, commonFounder, max_dv))
    
    return (probFinal)
  }
}
