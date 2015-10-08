library("aod")

# vim: set tabstop=2 shiftwidth=2 expandtab
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



# Rsquared sampling
sample.genotypes <- function(odds.ratio = 1, cases = 100, controls = 100, maf = 0.5) {
  # Generate alleles
  num.variant.alleles <- round(maf*2*(cases+controls))
  
  allele.pool <- sample(c(rep(T,num.variant.alleles), rep(F, (2*(cases+controls))-num.variant.alleles)), 2*(cases+controls))
  stopifnot(length(allele.pool) == 2*(cases+controls))
  # sampling way (not precise!)
  #allele.pool <- runif(2*(cases+controls)) < maf

  sample.alleles <- function(scale.factor, cases, controls, maf) {
    sample_probabilities <- ifelse(allele.pool==F, 2, scale.factor)
    # ensure sum(probs) = 1
    sample_probabilities <- sample_probabilities / sum(sample_probabilities)
    
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

    #print(paste("iteration =",iteration,"sample.odds=",sample.odds.ratio, "lower=",lower.scale.factor,"upper=", upper.scale.factor,"scale=",scale.factor, "diff=", sample.odds.ratio.diff))
    
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

scale.to.dosages <- function(unscaled) {
  minima <- min(unscaled) 
  maxima <- max(unscaled)
  2 * (unscaled + abs(minima)) / (maxima+abs(minima) - (minima+abs(minima)))
}

# The usual R jitter function adds noise from a uniform distributio
jitter.dosages <- function (dosages, amount = 0) 
  scale.to.dosages(dosages + (amount*rnorm(length(dosages))))


# Sample dosages from genotypes using fast binary search for r-squared 
# It will result in  approximately correct r-squared, but have a tendency to "look" to nice..
sample.dosages <- function(genotypes, info,epsilon=0.01) {
  jitter.amount = 0
  dosages <- genotypes
  
  repeat {
    rsq <- summary(lm(genotypes ~ dosages))$r.squared
    diff <- abs(rsq - info)
	
    #print(paste("info=", info,  "rsq=", rsq, "diff=", diff,"jitter.amount=", jitter.amount))

    if(diff < epsilon)
      break
    else if (rsq > info)
      jitter.amount <- jitter.amount + (diff / 2)
    else
      jitter.amount <- jitter.amount - (diff / 2)

    dosages <- jitter(genotypes,amount=jitter.amount)
  }
  list(
    # Make sure dosages are normalized to [0-2] interval:
   dosages = dosages,
#2 * (dosages + abs(min(dosages))) / ((max(dosages)+abs(min(dosages))) - (min(dosages)+abs(min(dosages)))),
    r.squared = rsq
  )
}

sample.dosages2 <- function(genotypes, info,epsilon=0.01) {
  dosages <- scale.to.dosages(runif(length(genotypes)))
  rsq = 0
  
  repeat {
    rsq <- summary(lm(genotypes ~ dosages))$r.squared
     diff <- abs(rsq - info)
    
    if(diff < epsilon)
      break
    if(rsq > info+epsilon)
      break
    

    # pick random "dosages" and update them to be closer to true genotypes
    # the amount of genotypes to update per iteration is proportional to 
    # to the level of info 
    how.many <- sum(runif(length(genotypes)) < info)
    picked <- sample(1:length(genotypes), how.many)
    dosages[picked] <- (dosages[picked] + genotypes[picked]) / 2

    print(paste("info=", info,  "rsq=", rsq, "diff=", diff))

  }
  list(dosages = dosages,r.squared = rsq)
}

# Test of sample.dosages2
test.sample.dosages2 <- function() {
  r.squared.values <- sapply(1:250, function(x) {
    info <- runif(1)
    maf <- runif(1)/2
    or <- runif(1)*4
    print(paste('or=',or,"maf=",maf))
    snp = sample.genotypes(or, 100, 100, maf)
    (info - sample.dosages2(c(snp$genotypes$cases, snp$genotypes$controls),info)$r.squared)
  })
  print(r.squared.values)
  hist(r.squared.values,breaks=20)
  qqnorm(r.squared.values)
}

snp.table <- function(sampled.snp)
  with(sampled.snp,
    data.frame(
      case = c(rep(T,cases), rep(F,controls)),
      genotype = c(genotypes$cases, genotypes$controls)))

dosages.table <- function(snp,info) {
  dosages <- sample.dosages2(c(snp$genotypes$cases, snp$genotypes$controls),info)
  data.frame(
    case = c(rep(T,snp$cases), rep(F,snp$controls)),
    genotype = dosages$dosages)
}

sample.snp <- function(maf, info, odds.ratio, n_cases, n_controls) {
  snp <- sample.genotypes(odds.ratio,n_cases,n_controls,maf)
  dosages <- sample.dosages(c(snp$genotypes$cases, snp$genotypes$controls),info)
  print(dosages$r.squared)
  d <- data.frame(
    case = c(rep(T,snp$cases), rep(F,snp$controls)),
    genotype = c(snp$genotypes$cases, snp$genotypes$controls),
    dosage = dosages$dosages)
}

test.snp <- function(d) {
  fit<-glm(case ~ dosage, data=d, family=binomial()) 
  test<-wald.test(b = coef(fit), Sigma=vcov(fit), Terms=1)
  test$result$chi2[3]
}


# Takes a dataframe with columns 'info' and 'maf'
calculate.power <- function(imputed, n_cases, n_controls, odds.ratio, max.maf=0.5, max.info=1, number.of.tests=NA, significance=0.05, multiple.test.adjust="bonferoni") {
  
  if (is.na(number.of.tests)) {
    cat("Assumming number.of.tests = nrow(imputed): ")
    cat("Note that multiple test correction will be inaccurate if nrow(imputed) is different from actual number imputed markers.")
  }
  
  maf.idx = which(colnames(imputed)=="maf")
  info.idx = which(colnames(imputed)=="info")
  
  p.values <- sapply(1:nrow(imputed), function(i) {
    snp = sample.snp(imputed[i,maf.idx], imputed[i,info.idx], odds.ratio, n_cases, n_controls)
    test.snp(snp)
  })
  
  # Return the proportion of signifant tests
  sum(p.adjust(p.values) < significance) / length(p.values)
}

# sample a data.frame of imputation qualities - (maf,info) pairs for used in power calculations
# This function samples  (maf,info) pairs with probability proportional to observing that 
# maf+info pair in the imputed data
sample.info.maf <- function(info.files, sampling.iterations) {
  info.file.names <- Sys.glob(info.files)
  
  # We record the frequency of observed info+maf imputations
  counts <- table(factor(levels=1:100/100), factor(levels=1:50/100)) 
  for(file in info.file.names) {
    print(paste("processing file:", file))
    info <- read.table(file,h=T)
    counts <- counts + table(
      factor(round(info$info,digits=2), levels=1:100/100),
      factor(round(ifelse(info$exp_freq_a1 > 0.5, 1-info$exp_freq_a1,info$exp_freq_a1),digits=2), levels=1:50/100))
  }
  frequencies <- counts / 5000
  
  sampled.imputation.quals <- data.frame(maf = rep(NA,sampling.iterations), info = rep(NA,sampling.iterations))
  iteration <- 0

  
  repeat {
    iteration <- iteration + 1
    if (iteration > sampling.iterations)
      break
    index <- sample(1:5000, 1,prob=frequencies)
    info_pct   <- ifelse(index %% 100 == 0, 100, index %% 100) # This would be so much simpler with 0-based indexing..
    maf_pct  <- (index %/% 100) + 1
    if (maf_pct == 51) maf_pct <- 50
    sampled.imputation.quals[iteration,] <- c(maf_pct/100, info_pct/100)
  }
  sampled.imputation.quals
}

impute2.pwr <- function(info.files, sampling.iterations, ...) {
  calculate.power(sampled.imputation.quals,...)
}

tst.calc.pwr <- function() {
  tests = 250
  cal
  # Let's asssume the LuCamp scenario of 1000/1000 cases/controls
  n_cases = 1000
  n_controls = 1000
  
  
  proposed.maf <- runif(tests)/2
  imputed.values <- data.frame(
    maf = ifelse(proposed.maf >= 0.01, proposed.maf, 0.01),
    info = runif(tests)
  )
  
  calculate.power(imputed.values,n_cases, n_controls, 1.35)
}

tst.calc.pwr2 <- function() {
  # Let's asssume the LuCamp scenario of 1000/1000 cases/controls
  n_cases = 1000
  n_controls = 1000
  imputed.table <- sample.info.maf(info.files="info-files/chr21.info", sampling.iterations = 100)
  powers <- sapply(seq(1,2,0.1), function(or) calculate.power(imputed.table,n_cases, n_controls, or))
  plot(seq(1,2,0.1), powers)
}

testit <- function() {
  tst.calc.pwr()
}