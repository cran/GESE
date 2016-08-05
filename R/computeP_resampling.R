computeP_resampling <-
function(prob, obsP, obsWeightStat, familyWeight, threshold=1e-8)
{	probSeg= prob[!is.na(prob)]
	if(!anyNA(familyWeight))
	{	familyWeightGene = familyWeight[!is.na(prob),]
	}else
	{	familyWeightGene = NA
	}
	
	if (threshold > 1e-5)
	{	nsim = 1/threshold
	}else
	{	nsim = 100000
	}	
	
	pv = oneSetSim(nsim, probSeg,obsP, obsWeightStat, familyWeightGene)
	limit = 1/threshold
	totalsim = pv$nsim
	pvalue = pv$pvalue
	pvalue_weighted = pv$pvalue_weighted
	
	if (!is.na(pvalue_weighted))
	{
	
		while((pvalue < 1/totalsim*10 | pvalue_weighted < 1/totalsim*10) & totalsim < limit)
		{	cat("P-value (", pvalue, ") or weighted p-value (", pvalue_weighted, " is less than ", 1/totalsim*10, ", running more simulations...\n")
			moresim = min(totalsim*9, limit-totalsim)
			ntimes = moresim/pv$nsim
			for(simset in 1:ntimes)
			{	newpv = oneSetSim(nsim, prob=probSeg, obsP=obsP, obsWeightStat=obsWeightStat, familyWeightGene=familyWeightGene)
				pvalue = (pvalue*totalsim+ newpv$pvalue*newpv$nsim)/(totalsim+newpv$nsim)
				pvalue_weighted = (pvalue_weighted*totalsim+ newpv$pvalue_weighted*newpv$nsim)/(totalsim+newpv$nsim)
				totalsim = totalsim + newpv$nsim
				cat("-")
			
			}		
			cat("\n")
		}
	}else
	{	while((pvalue < 1/totalsim*10) & totalsim < limit)
		{	cat("P-value (", pvalue, ")  is less than ", 1/totalsim*10, ", running more simulations...\n")
			moresim = min(totalsim*9, limit-totalsim)
			ntimes = moresim/pv$nsim
			for(simset in 1:ntimes)
			{	newpv = oneSetSim(prob=probSeg, obsP=obsP)
				pvalue = (pvalue*totalsim+ newpv$pvalue*newpv$nsim)/(totalsim+newpv$nsim)
				totalsim = totalsim + newpv$nsim
				cat("-")
			
			}	
			cat("\n")	
		}

	}

	return (list(pvalue=pvalue, nsim=totalsim, pvalue_weighted = pvalue_weighted))

}
