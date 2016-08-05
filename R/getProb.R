getProb <-
function(Seg, prob)
{	
	return (prod(prob^Seg*(1-prob)^(1-Seg)))
}
