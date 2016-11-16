#### This file is used to compute the gene-based and variant-based segregation information using recessive and compound heterozygous model

getSegInfo <- function(pednew, dataPed, mapInfo, mode="recessive")
{		
	#@@ pednew: data frame containing complete pedigree infomraiont for the subjects, the family IDs should be character variables. The required columns are FID, IID, faID, moID, and sex.
	#@@ variantInformation: this data frame should include the information for all the variants satisfying the same filtering criteria in the chosen reference genome. The columns should include at least, but not limited to: unique SNP ID (colname SNP), GENE name (GENE), minor allele frequency in the reference database (MAF).
	#@@ dataPed	plink raw file format, recoded to including only polymorphic variants and sequenced subjects. Also important to make sure the recoding is with respect to the minor allele in the population. The affection status of this file will be used as phenotype. The column names of this data frame should be: FID, IID, PAT, MAT, SEX, PHENOTYPE, (variant ids, do not need to match the SNP ids in the mapInfo parameter, but should be in the same order as the mapInfo variable).
	#@@ mapInfo: this data frame should include the information for the variants in the sequenced dat. The columns include at least, but not limited to: unique SNP ID (colname SNP), corresponding GENE name (GENE).
	#@@ threshold: the precision needed for the final p-value
	#@@ tol: the tolerance value set for numerical euqality comparison
	#@@ onlySeg: a logical value, whether only segregation information is needed. If FALSE, p-value for the tests will be calculated too
	#@@ resampling: a logical value, whether resampling is used to compute the p-values.
	#@@ nsegLimit: an integer value. For any gene with the number of families segregated less than this value, the segregating probability will be computed using recursive function implemented in computeP.so
	#@@ familyWeight: a data frame gives the weight for the families. If it is NA, no weighting scheme is used. Otherwise, its dimenstion should be (number of families)x(number of genes+1). The first column should be family name (column name FID). If the weights for the families are the same for all the genes, the second column should just be weight (columns name weight), otherwise the second column and above should be the gene names (columns names are GENE names).
	 
	#### check files formats
	if(!is.data.frame(dataPed))
	{	stop("the raw file of the study is not supplied correctly!\n")
	}
	if(ncol(dataPed)< 7|!all(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE") %in% colnames(dataPed)))
	{	stop("The required columns of dataPed include: FID, IID, PAT, MAT, SEX, PHENOTYPE, (SNP ids),...\n")
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
	allGenes = sort(unique(mapInfo$GENE))
	
	if (mode == "dominant")
	{
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
						{	controlNoncarrier <- apply(controls[,7:(ncol(famTemp))], 2,  function(x) sum(x==0, na.rm=TRUE))==nrow(controls)
						}
					}
					if(ncol(famTemp)==7)
					{	segInfoVar[,(i+2)] <- (sum(cases[,7:(ncol(famTemp))]>0,na.rm=TRUE)==nrow(cases) & controlNoncarrier )
					}else
					{	segInfoVar[,(i+2)] <- (apply(cases[,7:(ncol(famTemp))], 2, function(x) sum(x>0, na.rm=TRUE))==nrow(cases) & controlNoncarrier)
			
					}
					#segInfo[j,(i+1)] <- segOrNot
				}else
				{	segInfoVar[,(i+2)] <- NA
				}
		}
		colnames(segInfoVar)[-c(1,2)] <- fam
		segInfo <- aggregate(segInfoVar[,-c(1,2)], by=list(segInfoVar$GENE), FUN="any", simplify=TRUE, na.rm=TRUE)
		colnames(segInfo)[1] <- "GENE"
		ress = cbind(segInfo, apply(segInfo[,-1], 1, sum, na.rm=TRUE))
		colnames(ress)[ncol(ress)] <- "numSegFam"
		ress2 = ress[order(as.integer(ress$numSegFam), decreasing=TRUE),]
	
		ress3 = cbind(segInfoVar, apply(segInfoVar[,-c(1,2)], 1, sum, na.rm=TRUE))
		colnames(ress3)[ncol(ress3)] <- "numSegFam"
		ress4 = ress3[order(as.integer(ress3$numSegFam), decreasing=TRUE),]
		return (list(geneSeg = ress2, varSeg = ress4))
		
	}else if (mode == "recessive")
	{
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
					{	if(ncol(famTemp)==7)## if there is only one variant
						{	controlNoncarrier <-( sum(controls[,7]<=1, na.rm=TRUE) == nrow(controls))
						}else
						{	controlNoncarrier <- apply(controls[,7:(ncol(famTemp))], 2, function(x) sum(x<=1, na.rm=TRUE))==nrow(controls)
						}
					}
					if(ncol(famTemp)==7)
					{	segInfoVar[,(i+2)] <- (sum(cases[,7:(ncol(famTemp))]==2,na.rm=TRUE)==nrow(cases) & controlNoncarrier )
					}else
					{	segInfoVar[,(i+2)] <- (apply(cases[,7:(ncol(famTemp))], 2, function(x) sum(x==2, na.rm=TRUE))==nrow(cases) & controlNoncarrier)
			
					}
					#segInfo[j,(i+1)] <- segOrNot
				}else
				{	segInfoVar[,(i+2)] <- NA
				}
		}
		colnames(segInfoVar)[-c(1,2)] <- fam
		segInfo <- aggregate(segInfoVar[,-c(1,2)], by=list(segInfoVar$GENE), FUN="any", simplify=TRUE, na.rm=TRUE)
		colnames(segInfo)[1] <- "GENE"
		ress = cbind(segInfo, apply(segInfo[,-1], 1, sum, na.rm=TRUE))
		colnames(ress)[ncol(ress)] <- "numSegFam"
		ress2 = ress[order(as.integer(ress$numSegFam), decreasing=TRUE),]
	
		ress3 = cbind(segInfoVar, apply(segInfoVar[,-c(1,2)], 1, sum, na.rm=TRUE))
		colnames(ress3)[ncol(ress3)] <- "numSegFam"
		ress4 = ress3[order(as.integer(ress3$numSegFam), decreasing=TRUE),]
		return (list(geneSeg = ress2, varSeg = ress4))
	}else if (mode == "CH") ## compound heterozygous
	{		### consider for each family, only the variants present in all cases
			segInfoVar <- data.frame(mapInfo[,c("GENE", "SNP")], t(rep(NA, numFam)))
			segALL <- c()
			colnames(segInfoVar)[-c(1,2)] <- fam
			for( i in 1:numFam)
			{	cat("Family ", fam[i], "\n")
				famTemp <- dataPed3[dataPed3$FID==fam[i], ]
				## whether any variant segregate
				cases <- famTemp[famTemp$PHENOTYPE==1,]
				if(nrow(cases)>0)
				{   ## segInfoVar computes whether the variant is present in all cases in the family i
					if(ncol(famTemp)==7)
					{	segInfoVar[,(i+2)] <- (sum(cases[,7:(ncol(famTemp))]>0,na.rm=TRUE)==nrow(cases))
					}else
					{	segInfoVar[,(i+2)] <- (apply(cases[,7:(ncol(famTemp))], 2, function(x) sum(x>0, na.rm=TRUE))==nrow(cases))
			
					}
					
				}else
				{	segInfoVar[,(i+2)] <- NA
				}
				varTable <- data.frame(mapInfo[,c("GENE", "SNP")], t(famTemp[,-c(1:6)]))
				varTable2 = varTable[varTable$SNP %in% segInfoVar[segInfoVar[,(i+2)]==TRUE & !is.na(segInfoVar[,(i+2)]),]$SNP,]
				if (nrow(varTable2)<=1)
				{	cat("There is no pair of variants present in this family.\n")
					if (length(segALL)==0)
					{	segALL <- data.frame(ID=NA, GENE.x=NA, GENE.y=NA)
						segALL[,eval(fam[i])] <- FALSE	
					}else
					{ segALL[,eval(fam[i])] <- FALSE
					}
					next
				}
				varTable2$joinCol <- "ALL"
				varTable2 = varTable2[order(varTable2$SNP),]
				varTableCross <- merge(varTable2, varTable2, by="joinCol", sort=FALSE)
				varTableCross2 <- varTableCross[varTableCross$SNP.x<varTableCross$SNP.y,]

				varTableCross2$ID <- paste(varTableCross2$SNP.x, varTableCross2$SNP.y, sep="-")
				
				segInfoVar2 <- data.frame(ID=varTableCross2$ID, GENE.x=varTableCross2$GENE.x, GENE.y =varTableCross2$GENE.y,  FAM=NA)	
 				colnames(segInfoVar2)[-c(1,2,3)] <- fam[i]
				controlsTemp <- famTemp[famTemp$PHENOTYPE==0,]
				if (nrow(controlsTemp)==0)
				{	segInfoVar2[, 4] <- TRUE
				}else
				{
					for(j in 1:nrow(famTemp))
					{	if(famTemp$IID[j] %in% controlsTemp$IID)
						{
							subjName = famTemp[j,]$IID
							varTableCross2[,eval(subjName)] <- (varTableCross2[,3+j]>0 & varTableCross2[,ncol(varTable2)+2+j]>0)
						}
					}
					controls <- varTableCross2[,colnames(varTableCross2) %in% c("ID", controlsTemp$IID)]
					if (ncol(controls)==2)
					{	segInfoVar2[,(4)] <- sapply(controls[,-1], function(x) x==FALSE)
							
					}else
					{	segInfoVar2[,(4)]<- apply(controls[,-1], 1, function(x) sum(x==FALSE, na.rm=TRUE))== nrow(controlsTemp)
					}
				}
				if(length(segALL)==0)
				{	segALL <- segInfoVar2
				}else
				{	segALL <- merge(segALL, segInfoVar2, by=c("ID", "GENE.x", "GENE.y"), sort=F, all.x=TRUE, all.y=TRUE)
				}
		

			}
			segALL = segALL[!is.na(segALL$ID),]
			naValues = which(is.na(segALL), arr.ind=TRUE)
			segALL[naValues] = FALSE
			segALL$numSegFam <- apply(segALL[,-c(1,2,3)], 1, sum, na.rm=TRUE)
			segALL <- segALL[order(segALL$numSegFam, decreasing=TRUE),]
			
			### within one gene PER GENE info
			segALLsub <- segALL[segALL$GENE.x==segALL$GENE.y, -c(ncol(segALL))]
			if(nrow(segALLsub)>0)
			{
				segInfoGene <- aggregate(segALLsub[,-c(1,2,3)], by=list(segALLsub$GENE.x), FUN="any", simplify=TRUE, na.rm=TRUE)
				colnames(segInfoGene)[1] <- "GENE"
				ress = cbind(segInfoGene, apply(segInfoGene[,-1], 1, sum, na.rm=TRUE))
				colnames(ress)[ncol(ress)] <- "numSegFam"
				ress2 = ress[order(as.integer(ress$numSegFam), decreasing=TRUE),]
			}else
			{	ress2 = NA
			}
			
			### Per gene pair info
			segALL2 = segALL[, -c(ncol(segALL))]
			if(nrow(segALL2)>0)
			{
				segInfoGene2 <- aggregate(segALL2[,-c(1,2,3)], by=list(segALL2$GENE.x, segALL2$GENE.y), FUN="any", simplify=TRUE, na.rm=TRUE)
				ress3 = cbind(segInfoGene2, apply(segInfoGene2[,-c(1,2)], 1, sum, na.rm=TRUE))
				colnames(ress3)[ncol(ress3)] <- "numSegFam"
				ress4 = ress3[order(as.integer(ress3$numSegFam), decreasing=TRUE),]
			}else
			{	ress4 = NA
			}
			return (list(geneSeg = ress2, genePairSeg = ress4, varSeg = segALL))

	}else
	{	stop("This mode of inheritance is not implemented yet.")
	}

	
	

}