# An allele cause causes and odds ratio increase in the chance that
# the person carrying the allele is a case.
# p(case|allele) = odds.ratio * p(control|allele)
# p(control|allele) = maf in controls
# maf = (maf_ctrl * controls + maf_case * cases) / cases+controls

# p(alelle|case) = odds.ratio * p(allele|control)

# p(case|alelle) = odds.ratio * p(control|allele) # may be solvable
# ->
# p(allele|case) = (odds.ratio * p(control|allele)) * _p(allele) / _p(case)

# First, solve p(case|alelle) = odds.ratio * p(control|allele)
# p(case|allele) / p(control|allele) = odds.ratiol



'++'<-function(a,b,...)  { c(a,b,...)} 

# Rsquared sampling
sampleRawGenotype <- function(odds.ratio = 1, cases = 100, controls = 100, maf = 0.5) {
  # Generate alleles
  num.variant.alleles <- round(maf*2*(cases+controls))
  
  allele.pool <- c(rep(T,num.variant.alleles), rep(F, (2*(cases+controls))-num.variant.alleles))
  stopifnot(length(allele.pool) == 2*(cases+controls))
  # sampling way (not precise!)
  #allele.pool <- runif(2*(cases+controls)) < maf

  sample.alleles <- function(scale.factor, cases, controls, maf) {
    sample_probabilities <- ifelse(allele.pool==F, 1, scale.factor)
    print(unique(sample_probabilities))
    # ensure sum(probs) = 1
    sample_probabilities <- sample_probabilities / sum(sample_probabilities)
    print(unique(sample_probabilities))
    
    which.cases.alleles <- sample(1:length(allele.pool), size=cases*2, replace=FALSE, prob = sample_probabilities)
    which.controls.alleles <- which(!(1:length(allele.pool) %in% which.cases.alleles))

    list(
      cases = allele.pool[which.cases.alleles],
      controls = allele.pool[which.controls.alleles]
    )
  }
  
  lower.scale.factor = 0
  upper.scale.factor = odds.ratio * 2
  scale.factor = odds.ratio
  epsilon = 0.01 # One percent different is acceptable
  iteration = 1
  best.scale.factor = scale.factor
  sampled = NA 
  sample.odds.ratio = NA
  
  repeat {
    sampled = sample.alleles(scale.factor,cases,controls,maf)

    sample.odds.ratio = (sum(sampled$cases) / cases) / (sum(sampled$controls) / controls)
     
    sample.odds.ratio.diff <- abs(sample.odds.ratio - odds.ratio)

    print(paste("iteration =",iteration, 
                "sample.odds=",sample.odds.ratio, 
                "lower=",lower.scale.factor,
                "upper=", upper.scale.factor,
                "scale=",scale.factor, 
                "diff=", sample.odds.ratio.diff))
    
    if (sample.odds.ratio.diff < epsilon)
        break
    else if (sample.odds.ratio < odds.ratio)
      lower.scale.factor = scale.factor
    else if (sample.odds.ratio > odds.ratio)
      upper.scale.factor <- scale.factor
    
    scale.factor <- lower.scale.factor + (upper.scale.factor - lower.scale.factor)/2
    
    if (lower.scale.factor == upper.scale.factor) {
      break
      if(sample.odds.ratio.diff > epsilon) {
        lower.scale.factor<-lower.scale.factor - 1
        upper.scale.factor<-upper.scale.factor + 1
      }
    }
    
    iteration <- iteration + 1
    if (iteration > 10)
      break
  }

  # pairs diploid-alleles assorted according hardy-weinberg equibrilium principles iff alleles are in random order
  alleles2genotypes <- function(alleles)
    alleles[seq(1,length(alleles)-1,2)] + alleles[seq(2,length(alleles),2)]
  
  list(
    cases = cases,
    controls = controls,
    maf = maf,
    cases.maf = sum(sampled$cases) / (cases*2),
    controls.maf = sum(sampled$controls) / (controls*2),
    true.odds.ratio = odds.ratio,
    sample.odds.ratio = (sum(sampled$cases) / cases) / (sum(sampled$controls) / controls),
    genotypes = list(
      controls = alleles2genotypes(sampled$controls),
      cases = alleles2genotypes(sampled$cases)))
  
}

imputedGenotype <- function(genotype, info) {
  repeat {
    dosages <- jitter(geno,amount=0)
    
    lm
  }
}

normalize.dosages <- function(dosages) 2 * (dosages + abs(min(dosages))) / ((max(dosages)+abs(min(dosages))) - (min(dosages)+abs(min(dosages))))

snp.table <- function(sampled.snp)
  with(sampled.snp,
    data.frame(
      case = c(rep(T,cases), rep(F,controls)),
      genotype = c(genotypes$cases, genotypes$controls)))


sampleMarker <- function(odds.ratio = 1, cases = 100, controls = 100, rsquared = 1) {
}

