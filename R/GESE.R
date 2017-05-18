GESE <- function(pednew, variantInformation, dbSize, dataPed, mapInfo, threshold=1e-7, onlySeg=FALSE, familyWeight=NA)
{	#@@ pednew:  data frame containing complete pedigree infomraiont for the subjects, the family IDs should be character variables. The required columns are FID, IID, faID, moID, and sex.
	#@@ variantInformation: this data frame should include the information for all the variants satisfying the same filtering criteria in the chosen reference genome. The columns should include at least, but not limited to: unique SNP ID (colname SNP), GENE name (GENE), minor allele frequency in the reference database (MAF).
	#@@ dataPed	plink raw file format, recoded to including only polymorphic variants and sequenced subjects. Also important to make sure the recoding is with respect to the minor allele in the population. The affection status of this file will be used as phenotype. The column names of this data frame should be: FID, IID, PAT, MAT, SEX, PHENOTYPE, (variant ids, do not need to match the SNP ids in the mapInfo parameter, but should be in the same order as the mapInfo variable).
	#@@ mapInfo: this data frame should include the information for the variants in the sequenced dat. The columns include at least, but not limited to: unique SNP ID (colname SNP), corresponding GENE name (GENE).
	#@@ threshold: the precision needed for the final p-value
	#@@ onlySeg: a logical value, whether only segregation information is needed. If FALSE, p-value for the tests will be calculated too
	#@@ familyWeight: a data frame gives the weight for the families. If it is NA, no weighting scheme is used. Otherwise, its dimenstion should be (number of families)x(number of genes with weight). The first column should be family name (column name FID). If the weights for the families are the same for all the genes, the second column should just be weight (columns name weight), otherwise the second column and above should be the gene names (columns names are GENE names).
	 
	#### check files formats
	if(!is.data.frame(dataPed))
	{	stop("the raw file of the study is not supplied correctly!\n")
	}
	if(ncol(dataPed)< 7|!all(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE") %in% colnames(dataPed)))
	{	stop("The required columns of dataPed include: FID, IID, PAT, MAT, SEX, PHENOTYPE, SNP ids,...\n")
	}
	nvariants = ncol(dataPed)-6
	fam <- as.character(unique(dataPed$FID))
	numFam <- length(fam)
			
	if(!is.data.frame(pednew))
	{	stop("Complete pedigree information should be supplied!\n")
	}
	if(!all(c("FID", "IID", "faID", "moID", "sex") %in% colnames(pednew)))
	{	stop("The columns of pednew need to include: FID, IID, faID, moID, sex. \n")
	}
	if(!is.data.frame(variantInformation))
	{	stop("VariantInformation is missing!\n")
	}	
	if(!all(c("SNP", "GENE", "MAF") %in% colnames(variantInformation)))
	{	stop("The columns of variantInformation needs to include: SNP, GENE, MAF.\n")
	}
	
	if(!is.data.frame(mapInfo))
	{	stop("The mapInfo parameter is missing!\n")
	}
	if(nvariants!=nrow(mapInfo))
	{	stop("The number of variants in the raw file is not the same as in mapInfo!\n")
	}	
	if(!all(c("SNP", "GENE") %in% colnames(mapInfo)))
	{	stop("The columns of mapInfo need to include: SNP, GENE.\n")
	}
	nGenes = length(unique(mapInfo$GENE))
	if(is.data.frame(familyWeight))
	{	if(ncol(familyWeight)==2 & !all(c("FID", "weight") %in% colnames(familyWeight)) & nGenes!=1)
		{		stop("If there are two colulmns in familyWeight, the column names are required to be FID, weight.\n")
		}
		if( colnames(familyWeight)[1] != "FID")
		{	stop("The first column of familyWeight should be family IDs with column name FID.\n")
		}
		if (any(is.na(familyWeight)))
		{	stop("No missing value is allowed in familyWeight!\n")
		}
	

		
	}
	
	if(!all(unique(dataPed$PHENOTYPE) %in% c(NA, 0, 1)))
	{		dataPed$PHENOTYPE = ifelse(dataPed$PHENOTYPE %in% c(0, -9), NA, 	dataPed$PHENOTYPE-1)
			if(!all(unique(dataPed$PHENOTYPE) %in% c(NA, 1, 0)))
			{	stop("Phenotype conversion to NA, 1 and 0 is not correct!\n")
			}
	}
	subjects <- dataPed[,c(1:2, 6)]
	cat("Trimming the families...\n")
	cat("Trimming step 1: keep only one lineage \n")
	trim <- trim_oneLineage(seqSub=subjects, pednew)
	pednew2 <- trim$pedInfoUpdate
	dataPed2 = dataPed[dataPed$IID %in% trim$seqSubjUpdate$IID,]
	
	cat("Trimming step 2: remove unrelated subjects \n")
	subjects <- dataPed2[,c(1:2, 6)]
	subjects2 <- trim_unrelated(seqSub=subjects, pednew2)
	dataPed3 = dataPed2[dataPed2$IID %in% subjects2$IID,]
	cat("Subjects ", unique(dataPed$IID[!dataPed$IID %in% dataPed3$IID]), " were trimmed out.\n")
	
	cat("Imputing missing genotypes to be 0 ...\n")
	naValues = which(is.na(dataPed3), arr.ind=TRUE)
	imputeIndex = naValues[naValues[,2]>6,]
	dataPed3[imputeIndex] = 0
	
	cat("Calculating segregating information ...\n")
	pedigrees = kinship2::pedigree(id=pednew2$IID, dadid=pednew2$faID, momid = pednew2$moID, sex=pednew2$sex, famid=pednew2$FID)
	subjects <- dataPed3[,c(1:2,6)]
	fam <- as.character(unique(subjects$FID))
	numFam <- length(fam)
	segInfoVar <- data.frame(mapInfo[,c("GENE", "SNP")], t(rep(NA, numFam)))
	for(i in 1:numFam)
 	{		cat("Family ", fam[i], "\n")
  	    	famTemp <- dataPed3[dataPed3$FID==fam[i], ]
  	  		## whether any variant segregate
  	  		cases <- famTemp[famTemp$PHENOTYPE==1,]
  	  		if(nrow(cases)>0)
  	  		{
 	  			controls <- famTemp[famTemp$PHENOTYPE==0,]
 	  			controlNoncarrier <- NA
 	  			if (nrow(controls)==0)
 	  			{	controlNoncarrier <- TRUE
 	  			}else
 	  			{	if(ncol(famTemp)==7)
 	  				{	controlNoncarrier <-( sum(controls[,7]==0, na.rm=TRUE) == nrow(controls))
 	  				}else
 	  				{	controlNoncarrier <- sapply(controls[,7:(ncol(famTemp))],  function(x) sum(x==0, na.rm=TRUE))==nrow(controls)
 	  				}
 	  			}
 	  			if(ncol(famTemp)==7)
 	  			{	segInfoVar[,(i+2)] <- (sum(cases[,7:(ncol(famTemp))]>0,na.rm=TRUE)==nrow(cases) & controlNoncarrier )
 	  			}else
 	  			{	segInfoVar[,(i+2)] <- (apply(cases[,7:(ncol(famTemp))], 2, function(x) sum(x>0, na.rm=TRUE))==nrow(cases) & controlNoncarrier)
 	  		
 	  			}
 	  			
 	  		}else
 	  		{	segInfoVar[,(i+2)] <- NA
 	  		}
    }
    colnames(segInfoVar)[-c(1,2)] <- fam
    segInfo <- aggregate(segInfoVar[,-c(1,2)], by=list(segInfoVar$GENE), FUN="any", simplify=TRUE, na.rm=TRUE)
    colnames(segInfo)[1] <- "GENE"
    segInfo2 <- segInfo[order(segInfo$GENE),!is.na(segInfo[1,])]
    
    ress = cbind(segInfo, apply(segInfo[,-1], 1, sum, na.rm=TRUE))
    colnames(ress)[ncol(ress)] <- "numSegFam"
	ress2 = ress[order(as.integer(ress$numSegFam), decreasing=TRUE),]
	
	ress3 = cbind(segInfoVar, apply(segInfoVar[,-c(1,2)], 1, sum, na.rm=TRUE))
    colnames(ress3)[ncol(ress3)] <- "numSegFam"
	ress4 = ress3[order(as.integer(ress3$numSegFam), decreasing=TRUE),]
	
	### 11/03/2016
	### updating variant information
	### 05/17/2017: add the check 
	if(any(!mapInfo$SNP %in% variantInformation$SNP))
	{	notInDatabase = mapInfo[!mapInfo$SNP %in% variantInformation$SNP,]
		notInDatabase$MAF = 0
		variantInformation2 = rbind(variantInformation, notInDatabase)
	}else
	{	variantInformation2 = variantInformation
	}
	

    if(!onlySeg)
    {
  		cat("Calculating p-value...\n")
		results <- getPvalue_resampling(pedigrees, subjects, variantInformation2, dbSize, segInfo2, threshold, familyWeight)
		results$segregation <- ress2
		results$varSeg <- ress4
		return (results)
	}else
	{	return (list(segregation=ress2, varSeg = ress4))
	}
}
