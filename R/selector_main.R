
#install these packages first, if necessary
library(psych)
library(CompQuadForm)
library(reshape2)
library(lavaan)
library(ggplot2)
library(goftest)
library(MASS)

selector <- function(f, b.reps=1000, seed=1){
  #continuous case, we assume ML
  if(f@Options$estimator != "ML"){
    cat("Object f must be fitted with estimator= ML\n")
    return()
  }
    
  orig.sample <- lavInspect(f, what="data")
  startvalues <- parTable(f)
  df <- fitmeasures(f, "df")
  
  porigs <- test_pvalues(f)
  
  outputstring <- paste("\n Candidate p-values are: \n",paste(names(porigs), round(porigs,3), sep = ":", collapse = ",  "))
  cat(paste(outputstring, "\n"))
  #####
  # transform sample a la Bollen-Stine
  Sigma.hat <- lavaan:::computeSigmaHat(lavmodel = f@Model)
  sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[1]])
  S.inv.sqrt <- lav_matrix_symmetric_sqrt(f@SampleStats@icov[[1]])
  sample.transformed <- data.frame(as.matrix(orig.sample) %*% S.inv.sqrt %*% sigma.sqrt)
  colnames(sample.transformed) <- colnames(orig.sample)
  #####

  ## run bootstrap in serial
  boot.df <- run_bootstrap(X=seed, b.reps=b.reps, sample.transformed=sample.transformed, startvalues=startvalues)
  
  colnames(boot.df)<- names(porigs)

  distUniform= as.vector(sapply(boot.df,function(x) ad.test(x, null="punif")$statistic ))

  selected <- which.min(distUniform)
  selectedname <- colnames(boot.df)[selected]
  pselected <- porigs[selectedname]

  cat("\n The selector selects: ", selectedname, " with p-value", pselected, "\n")

  rownames(boot.df) = NULL
  melted = melt(boot.df, id.vars = NULL)
  melted$statistic <- melted$variable

  #plot of the p-values
  ggplot(melted, aes(value))+geom_histogram(binwidth=0.02)+facet_wrap(~statistic)

}



