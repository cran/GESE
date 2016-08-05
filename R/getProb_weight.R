getProb_weight <-
function(Seg, prob, familyWeightGene)
{      
    return (prod((prob^Seg*(1-prob)^(1-Seg))^familyWeightGene$weight))
}
