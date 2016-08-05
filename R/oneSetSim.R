oneSetSim <- function(nsim = 100000, prob, obsP, obsWeightStat=NA, familyWeightGene=NA)
{	
	temp = matrix(rbinom(n=length(prob)*nsim,size=1, prob=prob), nsim, length(prob), byrow=TRUE)
	pvalue = sum(apply(temp,1, getProb, prob) <= obsP)/nsim
	pvalue_weighted = NA
	cat("Running ", nsim, " simulations, obtain p-value: ", pvalue, "\n")
	if(!anyNA(familyWeightGene))
	{	pvalue_weighted = sum(apply(temp,1, getProb_weight, prob, familyWeightGene) <= obsWeightStat)/nsim
		#cat("Running ", nsim, " simulations, obtain weighted p-value: ", pvalue_weighted,"\n")
	
	}
	return (list(nsim=nsim, pvalue=pvalue, pvalue_weighted=pvalue_weighted))
	
}

