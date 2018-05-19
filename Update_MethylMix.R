########################################################
# change r-squared threshold for MethylMix and correct Beta
# estimation step
########################################################

MethylMix_ModelGeneExpression_ModifiedRsquared <- function(METcancer, GEcancer, CovariateData = NULL) {      
  
  # overlapping samples
  OverlapSamples = intersect(colnames(METcancer), colnames(GEcancer))
  cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
  GEcancer = GEcancer[, OverlapSamples, drop = FALSE]
  METcancer = METcancer[, OverlapSamples, drop = FALSE]
  if (!is.null(CovariateData)) CovariateData = CovariateData[OverlapSamples, , drop = FALSE]
  
  Rsquares = matrix(0, nrow = length(rownames(METcancer)), ncol = 1)    
  Genes = rownames(METcancer)  
  PvalueThreshold = 0.001
  RsquareThreshold = .05
  
  cat("Correlating methylation data with gene expression...\n")
  i <- NULL # to avoid "no visible binding for global variable" in R CMD check
  Rsquares = foreach::foreach(i = 1 : length(rownames(METcancer)), .combine = 'c') %dopar% {
    Rsq = 0
    tmpGene = unlist(strsplit(Genes[i], '---'))[1]
    pos = which(rownames(GEcancer) == tmpGene)
    if (length(pos) > 0) {
      if (!is.null(CovariateData)) {
        res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ] + factor(CovariateData))
        res.summary = summary(res)
        an = anova(res)
        if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold ) { 
          # Note: to also request  methylation effect to be bigger than tissue add above: & an$`Pr(>F)`[1] < an$`Pr(>F)`[2], but tissue effect is always bigger so now I remove it
          Rsq = res.summary$r.squared
        }
      } else {           
        res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ])
        res.summary = summary(res)
        if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold) {
          Rsq= res.summary$r.squared 
        }
      }
    }
    Rsq
  }
  Rsquares = matrix(Rsquares, ncol = 1)
  
  # Rsquare threshold
  FunctionalGenes = Genes[Rsquares > RsquareThreshold]
  cat("\nFound", length(FunctionalGenes), "transcriptionally predictive genes.\n")
  return(FunctionalGenes)
}

betaEst_3 <-function (Y, w, weights, lowerBound = -Inf, upperBound = Inf, correctionFailed = FALSE) {
  y=Y
  yobs = !is.na(y)
  if (sum(yobs) <= 1) 
    return(c(1, 1))
  y = y[yobs]
  w = w[yobs]
  weights = weights[yobs]
  N <- sum(weights * w)
  p <- sum(weights * w * y)/N
  v <- sum(weights * w * y * y)/N - p * p
  logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 1))
  if (sum(yobs) == 2) 
    return(exp(logab))
  opt <- try(optim(logab, RPMM::betaObjf, ydata = y, wdata = w, weights = weights, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")){
    return(c(1, 1))
  }
  
  a = exp(opt$par[1])
  b = exp(opt$par[2])
  
  if(!correctionFailed & ( (a < 1 & b < 1) | (a == 1 & b == 1) ) ){
    cat("\tMaking Correction...\n")
    coef = .50
    low = .01
    
    min = min(y[w != 0 & y > lowerBound])
    lowerBound = min
    pos = which(y <= lowerBound & w != 0 & weights > low)
    
    max = max(y[w != 0 & y < upperBound])
    upperBound = max
    pos = append(pos, which(y >= upperBound & w != 0 & weights > low))
    
    weights[pos] <- weights[pos] * coef
    
    if((lowerBound > max(y[w != 0]) | upperBound < min(y[w != 0]))){
      cat("\tCorrection Failed...\n")
      correctionFailed = TRUE
      weights <- rep(1:length(w))
      lowerBound = -Inf
      upperBound = Inf
    }
    return(betaEst_2(y, w, weights, lowerBound, upperBound, correctionFailed))
  }
  c(a, b)
}

MethylMix_Update <- function(){
  detach("package:MethylMixBeta", unload = TRUE) 
  
  environment(betaEst_3) <- environment(MethylMixBeta:::betaEst_2)
  attributes(betaEst_3) <- attributes(MethylMixBeta:::betaEst_2)
  assignInNamespace("betaEst_2", value = betaEst_3, ns = "MethylMixBeta")
  
  environment(MethylMix_ModelGeneExpression_ModifiedRsquared) <- environment(MethylMixBeta:::MethylMix_ModelGeneExpression)
  attributes(MethylMix_ModelGeneExpression_ModifiedRsquared) <- attributes(MethylMixBeta:::MethylMix_ModelGeneExpression)
  assignInNamespace("MethylMix_ModelGeneExpression", value = MethylMix_ModelGeneExpression_ModifiedRsquared, ns = "MethylMixBeta")
  
  library(MethylMixBeta)
}

MethylMix_Revert <- function(){
  detach("package:MethylMixBeta", unload = TRUE)
  library(MethylMixBeta)
}
