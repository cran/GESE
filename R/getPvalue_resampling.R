getPvalue_resampling <- function(pedigrees, subjects, variantInformation, dbSize, segregatingFam, threshold, familyWeight)
{ ## @@ pedigrees: pedigree list of pedigree objects for all the families
  ## @@ subjects: subject information contains FID, IID, PHENOTYPE
  ### 0 indicates control, 1 indicates cases
  ## @@ variantInformation: variantInformation for all the variants to consider for each gene
  ### columns SNP, GENE, MAF
  ### @@ segregatingFam: segregation information 
  ### columns are families, rows are gene names

  ### step 1: calculate the conditional segregating probability for each family
  cat("Step 1: calculate the conditional segregating probability for each family ...\n")
  fam <- as.character(unique(subjects$FID))
  numFam <- length(fam)
  condSegFam <- data.frame(fam=fam, condP=-1)
  for(i in 1:numFam)
  { 
    condSegFam$condP[i] <- (condSegProbF(pedigrees[as.character(fam[i])], subjects))
  }
  
  ### step 2: calculate the segregating probability for each family for each gene
  ### considering all the variants
  cat("Step 2: calculate the segregating probability for each family for each gene ...\n")
  segProbGene <- segregatingFam
  segProbGene[,-1] <- NA
  for(i in 1:nrow(segProbGene))
  { cat(segProbGene$GENE[i], "\n")
    variantMAF <- variantInformation[variantInformation$GENE==segProbGene$GENE[i], ]    
   
    for( j in 1:numFam)
    { condSegProb <- condSegFam[condSegFam$fam==fam[j], 2]
      if(!is.na(condSegProb))
      {
        pedTemp <- pedigrees[as.character(fam[j])]
        indicatorSeg <- 1
        segProbGene[i, j+1] <- segProb(pedTemp, subjects, variantMAF, indicatorSeg, condSegProb, dbSize)
        
      }
      
    }
  }
  
  ### step 3: calculate the p-value
  obsP <- data.frame(GENE=segregatingFam$GENE, obsProb=-1, obsWeightStat=NA)
  pvalue <- data.frame(GENE=segregatingFam$GENE, pvalue=-1, nsim=-1, pvalue_weighted=NA)
  pvalue_nfam <- data.frame(GENE=segregatingFam$GENE, nseg=-1)
  cat("Step 3: calculate the p-value for each gene ... \n")

  for(i in  1:nrow(segProbGene))
  {   
	  cat(segProbGene$GENE[i], "\n")
	  variantMAF <- variantInformation[variantInformation$GENE==segProbGene$GENE[i], ]
	  
	  if(!anyNA(familyWeight)) ### all weights needs to be specified, no missing value allowed
  		{	if(ncol(familyWeight) == 2)
  			{	matchFam  = match( colnames(segProbGene)[-1], familyWeight$FID)
  				weightFam <- familyWeight[matchFam,]
  				colnames(weightFam)[2] =  "weight"
  			}else ### different family weights for each gene 
  			{	matchGene = match(segProbGene$GENE[i], colnames(familyWeight))
  				if(anyNA(matchGene))
  				{	weightFam <- NA
  				}else
  				{
  					matchFam = match(colnames(segProbGene)[-1], familyWeight$FID)
  					weightFam <- familyWeight[matchFam, matchGene]
  					weightFam <- data.frame(FID=colnames(segProbGene)[-1], weight = weightFam)
  				}
  			}
  		}else
  		{  weightFam = NA
  		}

      ### observed prob
      probb <- sapply(2:(numFam+1), function(x) ifelse(segregatingFam[i, x]==0, 1-segProbGene[i, x], segProbGene[i, x] ))
      obsP[i, 2] <- prod(as.numeric(probb), na.rm=TRUE)
      if(!anyNA(weightFam))
      {	 obsP[i, 3] <- (prod(as.numeric(probb)^as.numeric(weightFam$weight), na.rm=TRUE))
      }
      
      thisGene <- segProbGene[i, 2:(numFam+1)]
      nseg <- sum(as.numeric(segregatingFam[i, -1]))
 	  pvalue_nfam[i,2] <- nseg
      if(nseg == 0)
	  {	 pvalue[i,2] <- 1
	  	 pvalue[i,3] <- NA
	  	 if (!anyNA(weightFam))
	  	 { pvalue[i,4] <- 1
	  	 }
	  }else
	  {  if(nrow(variantMAF)>0)
	  	 {
	  	 	dd <- computeP_resampling(thisGene, obsP[i,2], obsP[i,3], weightFam, threshold)
	 		pvalue[i, 2] <- dd$pvalue
	 	 	pvalue[i, 3] <- dd$nsim
	 	 	pvalue[i, 4] <- dd$pvalue_weighted
	 	 }else
	 	 {	pvalue[i, 2] <- 0
	 	 	pvalue[i, 3] <- NA
	 	 	pvalue[i, 4] <- 0
	 	 }
	 	 
	  }

  }
  geneDatabase <- unique(variantInformation$GENE)
  geneDatabase2 <- geneDatabase[!geneDatabase %in% pvalue$GENE]
  if(anyNA(familyWeight))
  {
  	results <- cbind(pvalue, obsP$obsProb, pvalue_nfam[,2])
  	colnames(results) <- c("GENE", "pvalue","numSim", "pvalue_weighted", "obs_prob",  "N_seg")
  	results2 <- results[order(results$pvalue), c("GENE", "obs_prob", "pvalue", "numSim", "N_seg")]
  	if (length(geneDatabase2)>0)
  	{	results3 <- data.frame(GENE=sort(geneDatabase2),  obs_prob=NA, pvalue=1, numSim=0, N_seg=NA)
  		results4 = rbind(results2, results3)
  	}else
  	{	results4 = results2
  	}
  	return (list(results = results4, condSegProb =condSegFam, segProbGene=segProbGene)) 
  }else
  {	results <- cbind(pvalue, obsP$obsProb, obsP$obsWeightStat, pvalue_nfam[,2])
  	colnames(results) <- c("GENE", "pvalue","numSim", "pvalue_weighted",  "obs_prob", "obs_weight_stat",  "N_seg")
  	results2 <- results[order(results$pvalue), c("GENE", "obs_prob", "pvalue", "obs_weight_stat", "pvalue_weighted", "numSim", "N_seg")]
  	if (length(geneDatabase2)>0)
  	{	results3 <- data.frame(GENE=sort(geneDatabase2),  obs_prob=NA, pvalue=1, obs_weight_stat=NA, pvalue_weighted=1, numSim=0, N_seg=NA)
  		results4 = rbind(results2, results3)
  	}else
  	{	results4 = results2
  	}

  	return (list(results = results4, condSegProb =condSegFam, segProbGene=segProbGene)) 

  } 
}
